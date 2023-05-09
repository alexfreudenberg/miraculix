/*
 Authors 
 Alexander Freudenberg, alexander.freudenberg@stads.de

 Copyright (C) 2022-2023 Alexander Freudenberg

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

#include "dgemm_compressed_cuda.h"
#include "reference_u2d.h"

enum Reference { reference_device, reference_mkl, reference_none };

int main() {
  const Reference ref = reference_mkl;
  int internal_err;

  using ElementA_ = uint8_t;
  using ElementB_ = cutlass::u4f64_t;
  ;
  using ElementC_ = double;
  using ElementAccumulator_ = double;

  using LayoutA_ = cutlass::layout::RowMajor;
  using LayoutB_ = cutlass::layout::ColumnMajor;
  using LayoutC_ = cutlass::layout::RowMajor;

  int device_count;
  cudaGetDeviceCount(&device_count);

  // Problem size 
  //TODO: make these args
  int length_m = 1024 * 10;
  int length_n = 16;
  int length_k = 1024 * 5;

  // Create a tuple of problem size for matrix multiplication
  // Terminology:
  //      Full: Full problem set in unpacked u8 x F64 format
  //      Packed: Full problem set in packed u2 x F&4 format
  //      Small: Smaller submatrix to be recycled for the full problem set to
  //      avoid long simulation times
  cutlass::gemm::GemmCoord problem_size_full(length_m, length_n, length_k);
  cutlass::gemm::GemmCoord problem_size_small(256, length_n, length_k);
  cutlass::gemm::GemmCoord problem_size_packed(length_m, length_n,
                                               length_k / 4);
  cutlass::gemm::GemmCoord problem_size_t_packed(length_m / 4, length_n,
                                                 length_k);
  ElementC_ alpha = ElementC_(1);
  ElementC_ beta = ElementC_(0);

  printf("Allocating tensors\n");

  // Initialize tensors using CUTLASS helper functions
  cutlass::HostTensor<ElementA_, LayoutA_> tensor_a_small(
      problem_size_small.mk());
  cutlass::HostTensor<ElementA_, LayoutA_> tensor_a_full(
      (ref == reference_none) ? problem_size_small.mk()
                              : problem_size_full.mk());
  cutlass::HostTensor<uint8_t, LayoutA_> tensor_a_packed(
      problem_size_packed.mk());
  cutlass::HostTensor<uint8_t, LayoutA_> tensor_a_t_packed(
      problem_size_t_packed.mk());
  cutlass::HostTensor<double, LayoutB_> tensor_b_full(problem_size_full.kn());
  cutlass::HostTensor<ElementB_, LayoutB_> tensor_b_packed(
      problem_size_packed.kn());
  cutlass::HostTensor<ElementC_, LayoutC_> tensor_c(problem_size_full.mn());
  cutlass::HostTensor<ElementC_, LayoutC_> tensor_d(problem_size_full.mn());
  cutlass::HostTensor<ElementC_, LayoutC_> tensor_d_t(problem_size_full.kn());
  cutlass::HostTensor<ElementC_, LayoutC_> tensor_ref_d(problem_size_full.mn());

  printf("Simulating tensor content\n");
  // Simulate tensors using CUTLASS helper functions
  cutlass::reference::host::TensorFillRandomUniform(
      tensor_a_small.host_view(), 1, ElementA_(3), ElementA_(0), -1);
  cutlass::reference::host::TensorFillRandomUniform(
      tensor_b_full.host_view(), 1, double(4), double(0), -1);
  // Fallback constant fill
  // cutlass::reference::host::TensorFill(
  //     tensor_a_small.host_view(), uint8_t(2));
  // cutlass::reference::host::TensorFill(
  //     tensor_b_full.host_view(), double(2.0));

  printf("Packing data 1\n");

  std::memcpy(tensor_a_full.host_data(), tensor_a_small.host_data(),
              sizeof(uint8_t) * problem_size_small.k() *
                  problem_size_small.m());
  std::memcpy(tensor_b_packed.host_data(), tensor_b_full.host_data(),
              sizeof(double) * problem_size_full.k() * problem_size_full.n());

  printf("Packing data 2\n");

#pragma omp parallel num_threads(omp_get_max_threads())
  for (long i = 0;
       i < long(problem_size_packed.k()) * long(problem_size_small.m()); i++) {
    tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4);
    tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4 + 1) << 2;
    tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4 + 2) << 4;
    tensor_a_packed.host_data(i) |= tensor_a_full.host_data(i * 4 + 3) << 6;
  }

  printf("Packing data 3\n");
  uint8_t *a_ptr_full     = tensor_a_full.host_data(),
          *a_ptr_packed   = tensor_a_packed.host_data(),
          *a_t_ptr_packed = tensor_a_t_packed.host_data();


// Copying data from the first cycle to the rest of the tensor - creating a
// cyclical vector
#pragma omp parallel num_threads(omp_get_max_threads())
  for (long i = 1; i < problem_size_full.m() / problem_size_small.m(); i++) {
    if (ref != reference_none)
      std::memcpy(
          a_ptr_full + i * problem_size_full.k() * problem_size_small.m(),
          a_ptr_full,
          sizeof(uint8_t) * problem_size_full.k() * problem_size_small.m());

    std::memcpy(
        a_ptr_packed + i * problem_size_packed.k() * problem_size_small.m(),
        a_ptr_packed,
        sizeof(uint8_t) * problem_size_packed.k() * problem_size_small.m());
  }
  
  printf("Transposing packed matrix\n");
// #pragma omp parallel num_threads(64)
//   for(long i = 0; i < problem_size_full.m(); i++){
//     for(long j = 0; j < problem_size_full.k(); j++){
//       long new_row = j * problem_size_full.m()/4;
//       long old_row = i * problem_size_packed.k();
//       a_t_ptr_packed[new_row + i] |= ((a_ptr_packed[old_row + j/4 ] >> (2 * (j%4) )) & 0x3 ) << (2 * (i%4));
//     }
//   }

  cutlass::reference::host::TensorFill(tensor_c.host_view());
  cutlass::reference::host::TensorFill(tensor_d.host_view());
  cutlass::reference::host::TensorFill(tensor_ref_d.host_view());

  // Copy data from host to GPU
  printf("Copying data to the device\n");
  void *GPU_obj;

  printf("Dimensions: m %d, k %d, n %d\n", length_m, length_k, length_n);
  internal_err = plink2gpu(
      (char *)tensor_a_packed.host_data(), (char *)tensor_a_t_packed.host_data(),
      length_m, length_k,  tensor_b_full.host_data(), length_n, &GPU_obj);
  if (internal_err != 0) {
    printf("Plink2GPU failed with error code %d\n", internal_err);
    return 1;
  }

  tensor_c.sync_device();
  tensor_d.sync_device();
  tensor_d_t.sync_device();
  tensor_ref_d.sync_device();

 
  // Start calculation on GPU
  printf("Start calculation\n");
  int err = dgemm_compressed_gpu(
      true,
      GPU_obj, 
      length_n,
      tensor_b_full.host_data(), 
      length_k,
      tensor_d.host_data(),
      length_n);
  if (err != 0)
    return 1;
  // int err = dgemm_compressed_gpu(
  //     true,
  //     GPU_obj, 
  //     length_n,
  //     tensor_b_full.host_data(), 
  //     length_k,
  //     tensor_d_t.host_data(),
  //     length_n);
  

  // Start reference calculation
  printf("Syncing reference tensors\n");
  tensor_a_full.sync_device();
  tensor_b_full.sync_device();
  cutlass::reference::host::TensorFill(tensor_c.host_view(), 0.0);
  tensor_c.sync_device();

  printf("Start reference calculation\n");
  if (ref == reference_device) {
    cutlass::reference::device::compute_gemm<
        uint8_t, LayoutA_, double, LayoutB_, ElementC_, LayoutC_,
        ElementAccumulator_, ElementAccumulator_>(
        problem_size_full, alpha, tensor_a_full.device_ref(),
        tensor_b_full.device_ref(), beta, tensor_c.device_ref(),
        tensor_ref_d.device_ref(), 0.0);
    tensor_ref_d.sync_host();
  } else if (ref == reference_mkl) {
    Reference_Gemm<double> ref_gemm(length_m, length_n, length_k);
    ref_gemm(tensor_ref_d.host_data(), tensor_a_full.host_data(),
             tensor_b_full.host_data());
  }

  printf("\n");
  // Calculate difference between reference calc and main calc
  double diff = cutlass::reference::host::TensorSumSqDiff(
      tensor_ref_d.host_view(), tensor_d.host_view());

  // Print difference between reference calc and main calc
  printf("Dimensions m=%d, n=%d, k=%d\n", problem_size_full.m(),
         problem_size_full.n(), problem_size_full.k());
  printf("Difference %.e, ref_d norm %10.1lf, d norm %10.1lf\n", diff,
         cutlass::reference::host::TensorNorm(tensor_ref_d.host_view()),
         cutlass::reference::host::TensorNorm(tensor_d.host_view()));
  if(diff > pow(10,-5)){
    for(int i = 0; i < 10; i++){
      printf("(%.3lf, %.3lf), ", *(tensor_ref_d.host_data()+i),*(tensor_d.host_data()+i));
    }
  }
  printf("\n");

  freegpu(&GPU_obj);
  return 0;
}