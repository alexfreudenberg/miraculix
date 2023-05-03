/***************************************************************************************************
 * Copyright (c) 2022 - 2023 Alexander Freudenberg
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **************************************************************************************************/

#include <chrono>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <inttypes.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <unistd.h>

// CUTLASS includes
#include "cutlass/arch/mma.h"
#include "cutlass/arch/wmma.h"
#include "cutlass/array.h"
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/layout/matrix.h"
#include "cutlass/numeric_types.h"
#include "cutlass/util/host_tensor.h"
#include "cutlass/util/reference/device/gemm.h"
#include "cutlass/util/reference/host/gemm.h"
#include "cutlass/util/reference/host/tensor_compare.h"
#include "cutlass/util/reference/host/tensor_copy.h"
#include "cutlass/util/reference/host/tensor_fill.h"
#include "cutlass/util/reference/host/tensor_reduce.h"
#include "cutlass/util/tensor_view_io.h"

#include "matvec_u8x4f64.h"
#include "reference_u2d.h"

#ifdef STANDALONE
int main()
{
    const Reference ref = reference_none;
    using ColumnMajor = cutlass::layout::ColumnMajor;
    using RowMajor = cutlass::layout::RowMajor;

    using ElementA_ = uint8_t;
    using ElementB_ = cutlass::u4f64_t;
    ;
    using ElementC_ = double;
    using ElementAccumulator_ = double;

    const int kStages = 2;

    using LayoutA_ = LayoutA__;
    using LayoutB_ = LayoutB__;
    using LayoutC_ = RowMajor;
    using OperatorClass_ = cutlass::arch::OpClassSimt;
    using ArchTag_ = cutlass::arch::Sm61;
    using ThreadblockShape_ = cutlass::gemm::GemmShape<32, 32, 16>; // Don't increase this -- too high memory requirements
    using WarpShape_ = cutlass::gemm::GemmShape<16, 16, 16>;
    using InstructionShape_ = cutlass::gemm::GemmShape<1, 1, 1>;
    const int AlignmentA = 1;
    const int AlignmentB = 1;
    const bool SplitKSerial = false;

    // Problem size
    const long length_m = 1024 * 50;
    const long length_n = 16;
    const long length_k = 1024 * 1500;

    // Create a tuple of problem size for matrix multiplication
    // Terminology:
    //      Full: Full problem set in unpacked u8 x F64 format
    //      Packed: Full problem set in packed u2 x F&4 format
    //      Small: Smaller submatrix to be recycled for the full problem set to avoid long simulation times
    cutlass::gemm::GemmCoord problem_size_full(length_m, length_n, length_k);
    cutlass::gemm::GemmCoord problem_size_small(256, length_n, length_k);
    cutlass::gemm::GemmCoord problem_size_packed(length_m, length_n, length_k / 4);

    ElementC_ alpha = ElementC_(1);
    ElementC_ beta = ElementC_(0);

    printf("Allocating tensors\n");

    // Initialize tensors using CUTLASS helper functions
    cutlass::HostTensor<ElementA_, LayoutA_> tensor_a_small(
        problem_size_small.mk());
    cutlass::HostTensor<ElementA_, LayoutA_> tensor_a_full((ref == reference_none) ? problem_size_small.mk() : problem_size_full.mk());
    cutlass::HostTensor<uint8_t, LayoutA_> tensor_a_packed(
        problem_size_packed.mk());
    cutlass::HostTensor<double, LayoutB_> tensor_b_full(
        problem_size_full.kn());
    cutlass::HostTensor<ElementB_, LayoutB_> tensor_b_packed(
        problem_size_packed.kn());
    cutlass::HostTensor<ElementC_, LayoutC_> tensor_c(
        problem_size_full.mn());
    cutlass::HostTensor<ElementC_, LayoutC_> tensor_d(
        problem_size_full.mn());
    cutlass::HostTensor<ElementC_, LayoutC_> tensor_ref_d(
        problem_size_full.mn());

    printf("Simulating tensor content\n");
    // Simulate tensors using CUTLASS helper functions
    cutlass::reference::host::TensorFillRandomUniform(
        tensor_a_small.host_view(), 1, ElementA_(3), ElementA_(0),
        -1);
    cutlass::reference::host::TensorFillRandomUniform(
        tensor_b_full.host_view(), 1, double(4), double(0),
        -1);
    // Fallback constant fill
    // cutlass::reference::host::TensorFill(
    //     tensor_a_small.host_view(), uint8_t(2));
    // cutlass::reference::host::TensorFill(
    //     tensor_b_full.host_view(), double(2.0));

    printf("Packing data 1\n");

    std::memcpy(tensor_a_full.host_data(), tensor_a_small.host_data(), sizeof(uint8_t) * problem_size_small.k() * problem_size_small.m());
    std::memcpy(tensor_b_packed.host_data(), tensor_b_full.host_data(), sizeof(double) * problem_size_full.k() * problem_size_full.n());

    printf("Packing data 2\n");

#pragma omp parallel num_threads(64)
    for (long i = 0; i < long(problem_size_packed.k()) * long(problem_size_small.m()); i++) {
        tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4);
        tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4 + 1) << 2;
        tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4 + 2) << 4;
        tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4 + 3) << 6;
    }

    printf("Packing data 3\n");
    uint8_t *a_ptr_full = tensor_a_full.host_data(),
            *a_ptr_packed = tensor_a_packed.host_data();

// Copying data from the first cycle to the rest of the tensor - creating a cyclical vector
#pragma omp parallel num_threads(64)
    for (long i = 1; i < problem_size_full.m() / problem_size_small.m(); i++) {
        if (ref != reference_none)
            std::memcpy(a_ptr_full + i * problem_size_full.k() * problem_size_small.m(), a_ptr_full, sizeof(uint8_t) * problem_size_full.k() * problem_size_small.m());

        std::memcpy(a_ptr_packed + i * problem_size_packed.k() * problem_size_small.m(), a_ptr_packed, sizeof(uint8_t) * problem_size_packed.k() * problem_size_small.m());
    }

    cutlass::reference::host::TensorFill(
        tensor_c.host_view());
    cutlass::reference::host::TensorFill(
        tensor_d.host_view());
    cutlass::reference::host::TensorFill(
        tensor_ref_d.host_view());

    // Copy data from host to GPU
    printf("Copying data to the device\n");
    tensor_a_packed.sync_device();
    tensor_b_packed.sync_device();
    tensor_c.sync_device();
    tensor_d.sync_device();
    tensor_ref_d.sync_device();

    // Declare Gemm problem
    using CutlassGemm =
        typename cutlass::gemm::device::Gemm<
            ElementA_, LayoutA_, ElementB_, LayoutB_, ElementC_, LayoutC_,
            ElementAccumulator_, OperatorClass_, ArchTag_, ThreadblockShape_,
            WarpShape_, InstructionShape_, cutlass::epilogue::thread::LinearCombination<ElementC_, 1, ElementAccumulator_, ElementC_>, typename cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
            kStages, AlignmentA, AlignmentB, SplitKSerial, cutlass::gemm::device::DefaultGemmConfiguration<OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_, ElementAccumulator_>::Operator, true>;

    // Define a CUTLASS GEMM type
    CutlassGemm gemm_operator;

    // Define CUTLASS GEMM arguments
    typename CutlassGemm::Arguments arguments {
        problem_size_packed, // <- problem size of matrix multiplication
        tensor_a_packed.device_ref(), // <- reference to matrix A on device
        tensor_b_packed.device_ref(), // <- reference to matrix B on device
        tensor_c.device_ref(), // <- reference to matrix C on device
        tensor_d.device_ref(), // <- reference to matrix D on device
        { alpha, beta } // <- tuple of alpha and beta
    };

    // Calculate CUTLASS GEMM workspace size
    size_t workspace_size = CutlassGemm::get_workspace_size(arguments);

    // Allocate workspace memory
    cutlass::device_memory::allocation<double> workspace(workspace_size);

    // Instantiate CUTLASS kernel depending on templates
    CutlassGemm gemm_op;

    // Test if problem can be implemented
    cutlass::Status status = gemm_op.can_implement(arguments);
    if (status != cutlass::Status::kSuccess)
        printf("Can't implement");

    // Initialize CUTLASS kernel with arguments and workspace pointer
    status = gemm_op.initialize(arguments, workspace.get());
    if (status != cutlass::Status::kSuccess)
        printf("Initialize");

    // Start calculation
    printf("Start calculation\n");
    int niter = 10; // Number of iterations of GEMM op
    auto start = std::chrono::high_resolution_clock::now();
    // Launch initialized CUTLASS kernel
    for (int i = 0; i < niter; i++) {
        status = gemm_op(); // Actual gemm op
        cudaDeviceSynchronize(); // Required to sync streams on the device
    }
    if (status != cutlass::Status::kSuccess)
        printf("Operation error %d\n", status);

    auto stop = std::chrono::high_resolution_clock::now();

    // Calculate calculatiom time
    std::chrono::duration<double> duration = (stop - start);
    printf("Duration perf %.5lf\n", duration.count() / niter);

    // Copy back to host
    tensor_d.sync_host();

    // Start reference calculation
    printf("Start reference calculation\n");

    tensor_a_full.sync_device();
    tensor_b_full.sync_device();
    cutlass::reference::host::TensorFill(tensor_c.host_view(), 0.0);
    tensor_c.sync_device();

    start = std::chrono::high_resolution_clock::now();
    if (ref == reference_device) {
        cutlass::reference::device::compute_gemm<
            uint8_t, LayoutA_, double, LayoutB_, ElementC_, LayoutC_,
            ElementAccumulator_, ElementAccumulator_>(
            problem_size_full, alpha, tensor_a_full.device_ref(), tensor_b_full.device_ref(), beta,
            tensor_c.device_ref(), tensor_ref_d.device_ref(), 0.0);
        tensor_ref_d.sync_host();
    } else if (ref == reference_mkl) {
        Reference_Gemm<double> ref_gemm(length_m, length_n, length_k);
        ref_gemm(tensor_ref_d.host_data(), tensor_a_full.host_data(), tensor_b_full.host_data());
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = (stop - start);
    printf("Duration ref %.5lf\n", duration.count());

    for (int i = 0; i < min((long)10, length_n * length_m); i++) {
        printf("(%.3lf, %.3lf) ", tensor_d.host_data(i), tensor_ref_d.host_data(i));
    }

    printf("\n");
    // Calculate difference between reference calc and main calc
    double diff = cutlass::reference::host::TensorSumSqDiff(
        tensor_ref_d.host_view(), tensor_d.host_view());

    // Print difference between reference calc and main calc
    printf("Dimensions m=%d, n=%d, k=%d\n", problem_size_full.m(), problem_size_full.n(), problem_size_full.k());
    printf("Difference %.lf, ref_d norm %10.lf, d norm %10.1lf\n", diff,
        cutlass::reference::host::TensorNorm(tensor_ref_d.host_view()),
        cutlass::reference::host::TensorNorm(tensor_d.host_view()));
    printf("\n");

    return 0;
}
#endif