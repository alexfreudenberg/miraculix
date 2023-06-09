/*
 Alexander Freudenberg, alexander.freudenberg@stads.de

 Copyright (C) 2020-2023 Alexander Freudenberg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// Include order is important!!!
// Namespace conflicts if order is changed

#include <inttypes.h>
#include <thrust/device_vector.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

// CUTLASS includes
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/system/cuda/memory.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/discard_iterator.h>

#include "cuda_utils.h"
#include "mmagpuIntern.h"

static void gpuCrossprodIntern(unsigned int *CGM, size_t snps,
                               size_t individuals, double *ans,
                               size_t TileSize) {
  bool verbose = get_print_level() >= 0;
  const char *omp_num_threads = getenv("OMP_NUM_THREADS");
  int num_threads = 4;
  if (omp_num_threads != NULL) {
    num_threads = atoi(omp_num_threads);
  }
  if (verbose) {
    char print_message[] = (omp_num_threads != NULL)
                               ? ("OMP_NUM_THREADS is set")
                               : ("OMP_NUM_THREADS is not set");
    printf("%s, using %d threads.", print_message, num_threads);
  }
    //
    const size_t BytesPerRow =
        (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / 32) * 32;
    const size_t IntsPerRow = 1 + (BytesPerRow - 1) / sizeof(unsigned int);

    // sanity checks
    // limit Tilesize to individuals
    TileSize = TileSize > individuals ? individuals : TileSize;

    // calculates total memory requirements
    size_t req_mem = 2 * BytesPerRow * TileSize * CodesPerByte +
                    TileSize * TileSize * sizeof(unsigned int);
    if (req_mem > free_mem)
        ERR("Not enough memory available.");

    // Input data
    cutlass::uint4b_t *d_Z_block1;
    cutlass::uint4b_t *d_Z_block2;
    // Buffer for output
    int32_t *d_M;
    // Buffer for copying back results from device to host
    int32_t *h_M;

    const int size_of_input =
        BytesPerRow * TileSize * CodesPerByte / MEMORY_FACTOR;
    const int size_of_output = sizeof(int32_t) * TileSize * TileSize;
    // Initialization of buffers: We calculate n_streams of tile matrix
    // multiplications in parallel and allocate the corresponding amount of memory
    cudaMalloc((void **)&d_Z_block1, n_streams * size_of_input);
    cudaMalloc((void **)&d_Z_block2, n_streams * size_of_input);
    cudaMalloc((void **)&d_M, n_streams * size_of_output);
    cudaMallocHost((void **)&h_M, n_streams * size_of_output);
    err_check("Memory allocation: ");

    // initialization of cutlass gemm operators
    using ColumnMajor = cutlass::layout::ColumnMajor;
    using RowMajor = cutlass::layout::RowMajor;
    using TensorOp = cutlass::arch::OpClassTensorOp;
    using ElementA_ = cutlass::uint4b_t;
    using LayoutA_ = RowMajor;
    using ElementB_ = cutlass::uint4b_t;
    using LayoutB_ = ColumnMajor;
    using ElementC_ = int32_t;
    using LayoutC_ = RowMajor;
    using ElementAccumulator_ = ElementC_;
    using OperatorClass_ = TensorOp;
    using ArchTag_ = cutlass::arch::Sm75;
    using ThreadblockShape_ =
        typename cutlass::gemm::device::DefaultGemmConfiguration<
            OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
            ElementAccumulator_>::ThreadblockShape;
    using WarpShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::WarpShape;
    using InstructionShape_ =
        typename cutlass::gemm::device::DefaultGemmConfiguration<
            OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
            ElementAccumulator_>::InstructionShape;
    using EpilogueOutputOp_ =
        typename cutlass::gemm::device::DefaultGemmConfiguration<
            OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
            ElementAccumulator_>::EpilogueOutputOp;
    using ThreadblockSwizzle_ =
        typename cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>;
    const int Stages = cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::kStages;
    const int AlignmentA = cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::kAlignmentA;
    const int AlignmentB = cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::kAlignmentB;
    const bool SplitKSerial = false;
    using Operator_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::Operator;
    const bool IsBetaZero = false;

    using CutlassGemm = cutlass::gemm::device::Gemm<
        ElementA_, LayoutA_, ElementB_, LayoutB_, ElementC_, LayoutC_,
        ElementAccumulator_, OperatorClass_, ArchTag_, ThreadblockShape_,
        WarpShape_, InstructionShape_, EpilogueOutputOp_, ThreadblockSwizzle_,
        Stages, AlignmentA, AlignmentB, SplitKSerial, cutlass::arch::CustomOp,
        IsBetaZero>;

    // Define a CUTLASS GEMM type
    CutlassGemm gemm_operator;

    // Main loop
    // Calculates matrix multiplications in parallel: Each thread in this loop sends
    // its data to a different stream on the device. The threads calculate
    // concurrently and send the output back to main memory. Memory copies are
    // asynchronous to take full advantage of the memory bandwidth.
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(                                          \
            n_streams < 1 + (individuals - 1) / TileSize                           \
                    ? n_streams                                                    \
                    : 1 + (individuals - 1) / TileSize) schedule(dynamic)
    #endif
    for (int64_t i = 0; i < individuals; i += TileSize) {
        int threadidx = omp_get_thread_num();
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        cudaStreamSynchronize(stream);

        if (err != cudaSuccess)
        ERR("Stream couldn't be created");

        // Pointer to the first element of current rows
        unsigned int *x = (CGM + i * IntsPerRow);
        cutlass::uint4b_t *d_tile1 =
            d_Z_block1 + threadidx * TileSize * BytesPerRow;
        cutlass::uint4b_t *d_tile2 =
            d_Z_block2 + threadidx * TileSize * BytesPerRow;

        // Number of rows in matrix
        size_t const rows_left = individuals - i;
        // Size x of current tile
        size_t const x_tile_size = TileSize < rows_left ? TileSize : rows_left;

        cudaMemcpyAsync(d_tile1, x, x_tile_size * BytesPerRow,
                        cudaMemcpyHostToDevice, stream);

        cudaStreamSynchronize(stream);
        err_check("Copy 1:");

        // Inner loop
        for (int64_t j = i; j < individuals; j += TileSize) {

        // Same as above with y
        size_t const columns_left = individuals - j;
        size_t const y_tile_size =
            TileSize < columns_left ? TileSize : columns_left;
        unsigned int *y = (CGM + j * IntsPerRow);

        cudaMemcpyAsync(d_tile2, y, y_tile_size * BytesPerRow,
                        cudaMemcpyHostToDevice, stream);
        err_check("Copy 2:");
        cudaStreamSynchronize(stream);

        // initialize gemm arguments
        CutlassGemm::Arguments args(
            {int(x_tile_size), int(y_tile_size),
            int(BytesPerRow * CodesPerByte / COMPRESSION_GPU)},
            {d_tile1, int(BytesPerRow * CodesPerByte / COMPRESSION_GPU)},
            {d_tile2, int(BytesPerRow * CodesPerByte / COMPRESSION_GPU)},
            {d_M + threadidx * TileSize * TileSize, int(y_tile_size)},
            {d_M + threadidx * TileSize * TileSize, int(y_tile_size)}, {1, 0});
        cudaStreamSynchronize(stream);

        // compute Multiplication
        cutlass::Status status;
    #pragma omp critical
        status = gemm_operator(args, nullptr, stream);

        cudaStreamSynchronize(stream);
        err_check("Calculation:");

        // Copy results back to host
        cudaMemcpyAsync(h_M + threadidx * TileSize * TileSize,
                        d_M + threadidx * TileSize * TileSize,
                        TileSize * TileSize * sizeof(int32_t),
                        cudaMemcpyDeviceToHost, stream);
        err_check("Copying back:");

        cudaStreamSynchronize(stream);

        if (*(h_M + threadidx * TileSize * TileSize) == 0) {
            printf("Computation failed at thread %d, (%d,%d)\n", threadidx, i, j);
            print_kernel<<<1, 1>>>((int32_t *)d_M +
                                threadidx * TileSize * TileSize);
            j -= TileSize;
            continue;
        }
        err_check("Copy back:");

    // Loop over tile and store values in output matrix
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(n_streams) schedule(static)
    #endif
        for (int64_t di = 0; di < x_tile_size; ++di) {
            for (int64_t dj = 0; dj < y_tile_size; ++dj) {
            // Get result
            const auto Mij =
                *(h_M + threadidx * TileSize * TileSize + dj + di * y_tile_size);
            // Create pointers to the output matrix (because it is symmetric we
            // use two pointers)
            double *ans0 = ans + (i + di), *ans1 = ans + (i + di) * individuals;
            // Store result in ouput matrix
            ans0[(j + dj) * individuals] = ans1[j + dj] = (double)Mij;
            }
        }
        }
        cudaStreamDestroy(stream);
    }

    // Free memory
    cudaFree(d_Z_block1);
    cudaFree(d_Z_block2);
    cudaFree(d_M);
    cudaFreeHost(h_M);
}

static void crossprodIntern(Uint *CM, Uint snps, Uint individuals,
        double *ans) {
// tilse_size needs to stay the same: for smaller values we experience undocumented calculation failures on the device
const size_t tilesize = 2048; 

// Initialize host pointers and copy input data cuda managed memory
Uint* h_CM; 
const size_t BytesPerIndiv = UnitsPerIndiv256(snps) * BytesPerUnit;  cudaMallocHost((void**)&h_CM, individuals * BytesPerIndiv);
MEMCOPY(h_CM, CM, individuals * BytesPerIndiv);

gpuCrossprodIntern(h_CM, snps, individuals, ans, tilesize);
cudaFreeHost(h_CM);

}


extern "C" {
void crossprod_mmagpu(Uint *CGM, Uint snps, Uint individuals, double *ans) {
crossprodIntern(CGM, snps, individuals, ans);
}
}
