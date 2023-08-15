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
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <thrust/fill.h>
#include <thrust/device_vector.h>

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>

#include "dgemm_compressed_cuda.h"
#include "cuda_utils.h"


// 
// plink2gpu function
// 
int plink2gpu(char *genotype, char *genotype_transposed, int snps,
              int indiv, double *f, int n, void **GPU_obj) {
  /*
    Moves SNP matrix, its transposed and the according allele frequencies to the
    device and stores the pointer to this data in a separate object, a pointer
    to which is then returned.

    Parameters:
      genotype: Pointer to SNP matrix in plink format of size ceil(indiv/4) * snps
    +3
      genotype_transpoed: Pointer to transposed SNP matrix in plink format
      f: Pointer to vector of allele frequencies
      snps: Pointer to number of snps
      indiv: Pointer to number of individuals
      GPU_obj: void pointer to a pointer in which the GPU object is stored
  */

  // Print compile info 
  print_compile_info("dgemm_compressed");
  
  //
  // Initialize CUDA variables
  //
  cudaError_t err;
  uint8_t *d_genotype, *d_genotype_transposed;
  double *d_f, *d_unit, *d_C, *d_D;
  double *d_B;

  int n_datasets = int(genotype != NULL) + int(genotype_transposed != NULL);
  long n_bytes_per_snp =
      (indiv - 1) / 4 + 1; // number of columns of Z if individuals
                           // are zero-padded to be a multiple of 4
  long n_bytes_per_indiv =
      (snps - 1) / 4 +
      1; // number of columns of Z^T if SNPs are zero-padded to be a multiple of 4
  long size_buffer = 4 * ((long(max(snps, indiv)) - 1) / 4 + 1) * long(n);
  // Maximal size of the matrices B and C on the device
  // Matrices are forced to have a number of rows which is a multiple of 4 by
  // zero-padding This allows us to deal with SNP matrices with unaligned
  // dimensions which are themselves zero-padded

  debug_info("Dimensions: (%d,%d), size in bytes: %ld", snps, indiv, n_bytes_per_snp * long(snps) + long(indiv) * n_bytes_per_indiv);

  // Check if CUDA installation is correct
  if (checkCuda() != 0) {
    return 1;
  }
  // Switch to correct device
  int device = switchDevice();
  if (device == -1) {
    return 1;
  }
  // Check if a genotype dataset was supplied
  if (n_datasets == 0) {
    checkError(__func__, __LINE__, cudaErrorInvalidHostPointer);
    return 1;
  }

  // Check if enough memory is available
  size_t required_mem = 3 * size_buffer * sizeof(double);
  if (genotype != NULL)
    required_mem += n_bytes_per_snp * long(snps);
  if (genotype_transposed != NULL)
    required_mem += n_bytes_per_indiv * long(indiv);

  if (checkDevMemory(required_mem) != 0) {
    return 1;
  }

  //
  // Allocate device memory
  //
  if (genotype != NULL) {
    err = cudaMalloc((void **)&d_genotype, n_bytes_per_snp * long(snps));
    if (checkError(__func__, __LINE__, err) != 0)
      return 1;
  }
  if (genotype_transposed != NULL) {
    err = cudaMalloc((void **)&d_genotype_transposed,
                     n_bytes_per_indiv * long(indiv));
    if (checkError(__func__, __LINE__, err) != 0)
      return 1;
  }
  err = cudaMalloc((void **)&d_f, sizeof(double) * snps);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaMalloc((void **)&d_unit, sizeof(double) * indiv);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaMalloc((void **)&d_B, sizeof(double) * size_buffer);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaMalloc((void **)&d_C, sizeof(double) * size_buffer);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaMalloc((void **)&d_D, sizeof(double) * size_buffer);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;

  //
  // Copy data to device
  //
  if (genotype != NULL) {
    err = cudaMemcpy(d_genotype, genotype, long(n_bytes_per_snp) * long(snps),
                     cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0)
      return (1);
  }
  if (genotype_transposed != NULL) {
    err = cudaMemcpy(d_genotype_transposed, genotype_transposed,
                     long(n_bytes_per_indiv) * long(indiv),
                     cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0)
      return 1;
  }
  err = cudaMemcpy(d_f, f, sizeof(double) * long(snps), cudaMemcpyHostToDevice);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;

  // Fill d_unit with 1.0s
  thrust::device_ptr<double> d_unit_thrust(d_unit);
  thrust::fill(d_unit_thrust, d_unit_thrust + indiv, 1.0);
  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;

  //
  // Initialize GPU_gemm_storage object
  //
  struct GPU_gemm_storage *GPU_storage_obj =
      (struct GPU_gemm_storage *)malloc(sizeof(struct GPU_gemm_storage));

  GPU_storage_obj->d_genotype            = d_genotype;
  GPU_storage_obj->d_genotype_transposed = d_genotype_transposed;
  GPU_storage_obj->d_f                   = d_f;
  GPU_storage_obj->d_unit                = d_unit;
  GPU_storage_obj->d_B                   = d_B;
  GPU_storage_obj->d_C                   = d_C;
  GPU_storage_obj->d_D                   = d_D;
  GPU_storage_obj->size_buffer           = size_buffer;
  GPU_storage_obj->snps                  = snps;
  GPU_storage_obj->indiv                 = indiv;
  GPU_storage_obj->device                = device;

  // Set pointer to initialized object
  *GPU_obj = (void *)GPU_storage_obj;

  return 0;
}


// 
// freegpu function
// 
int freegpu(void **GPU_obj){
  cudaError_t err;

  if (checkCuda() != 0) {
    return 1;
  }
  if (*GPU_obj == NULL)
    return checkError(__func__, __LINE__, cudaErrorInvalidHostPointer);

  // Free device memory and derefence storage object
  struct GPU_gemm_storage *GPU_storage_obj = (struct GPU_gemm_storage *) (*GPU_obj);

  if (switchDevice(GPU_storage_obj->device) == -1) {
    return 1;
  }

  if (GPU_storage_obj->d_genotype != NULL) {
    err = cudaFree(GPU_storage_obj->d_genotype);
    if (checkError(__func__, __LINE__, err) != 0)
      return 1;
  }
  if (GPU_storage_obj->d_genotype_transposed != NULL) {
    err = cudaFree(GPU_storage_obj->d_genotype_transposed);
    if (checkError(__func__, __LINE__, err) != 0)
      return 1;
  }
  err = cudaFree(GPU_storage_obj->d_f);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaFree(GPU_storage_obj->d_unit);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaFree(GPU_storage_obj->d_B);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaFree(GPU_storage_obj->d_C);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaFree(GPU_storage_obj->d_D);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  free(GPU_storage_obj);

  GPU_obj = NULL;
  return 0;
}

// 
// dgemm_compressed_gpu function
//
int dgemm_compressed_gpu(bool transA, void *GPU_obj, int n, double *B, int ldb,
                         int centered, int normalized, double *C, int ldc) {
  /*
  Performs one of the operations
      C <- alpha * op(M - 2 * 1_indiv f^T) * op(B) + beta * C
  on the GPU, where
      op(X) = X or op(X) = X^T,
      alpha and beta are scalars,
      M is a compressed genotype matrix of dimensions indiv times snps stored in row-major,
      f is a vector of allele frequencies of length snps,
      op(B) is a matrix of double precision and number of rows equal to number of columns of op(M) 
      C is a matrix of double precision

  Parameters:
      transa: Specifies the form of op(M) used in the matrix multiplication. If
  transa = true, then op(M) = M^T. If transa = false, then op(M) =
  M.
      m: Specifies the number of rows of op(M - 2* 1_k*f^T) and C
      n: Specifies the number of columns of op(B) and C
      k: Specifies the number of columns of op(M - 2* 1_k*f^T) and rows of
  op(B) k1: Specifies the number of columns of op(M) in compressed format
      alpha: Not supported currently, only alpha=1 is allowed. Specifies the
  scalar alpha Storage: Struct which stores device pointers to both M and its
  transposed as well as a device pointer to the vector of allele frequencies f
      lda: Not supported currently
      B: Specifies the matrix B
      ldb: Not supported currently
      beta: Specifies the scalar beta. When beta is
      equal to zero, then C need not be set on input.
      C: Specifies the matrix C
      ldc: Not supported currently


  A further boost in performance can be achieved if the whole PCG is transfered
  to the GPU to avoid data movement.
  */

  // 
  // Initialization
  //

  // Initialize helper variables
  cudaError_t err;
  cublasStatus_t cublas_status;
  cublasHandle_t cublas_handle;

  if (GPU_obj == NULL)
    return checkError(__func__, __LINE__, cudaErrorInvalidHostPointer);

  struct GPU_gemm_storage *GPU_storage_obj = (struct GPU_gemm_storage *) GPU_obj;

  // Initialize device pointer for M
  // cutlass function assumes row major for M and PLINK bed uses SNP major, hence
  // pointer genotype is needed for transA = 't' and genotype_transposed if transA =
  // 'N'
  uint8_t          *d_M         = transA ? GPU_storage_obj->d_genotype : GPU_storage_obj->d_genotype_transposed;
  cutlass::u4f64_t *d_B         = reinterpret_cast<cutlass::u4f64_t *>(GPU_storage_obj->d_B);
  double           *d_f         = GPU_storage_obj->d_f,
                   *d_unit      = GPU_storage_obj->d_unit;
  double           *d_C         = GPU_storage_obj->d_C;
  double           *d_D         = GPU_storage_obj->d_D;
  double           *d_workspace = NULL;

  long m           = transA ? GPU_storage_obj->snps : GPU_storage_obj->indiv;
  long k           = transA ? GPU_storage_obj->indiv : GPU_storage_obj->snps;
  long size_buffer = GPU_storage_obj->size_buffer;
  long k1          = (k - 1) / 4 + 1;

  const  double alpha = 1.0,
         alpha_n2     = -2.0,
         beta         = 0.0;

  // Check CUDA installation
  if(checkCuda() != 0){
    return 1;
  }
  // Switch to device of object
  if (switchDevice(GPU_storage_obj->device) == -1) {
    return 1;
  }

  // Check correctness of pointers
  if (d_M == NULL) {
    printf(
        "Storage object does not hold data for transpose operation %s\n",
        transA ? "true" : "false");
    checkError(__func__, __LINE__, cudaErrorInvalidDevicePointer);
    return 1;
  }

  debug_info("\tEntering GPU multiplication\n");
  debug_info("Pointer: d_M %d, Dimensions: m %ld, k %ld, k1 %ld, n %ld", d_M, m, k, k1, n);

  // Create cuBLAS handle
  cublas_status = cublasCreate(&cublas_handle);
  if (checkError(__func__, __LINE__, cublas_status) != 0)
    return 1;

  // Zero-fill matrices C and B to avoid spurious results
  err = cudaMemset(d_B, 0, sizeof(double) * size_buffer);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaMemset(d_C, 0, sizeof(double) * err);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  err = cudaMemset(d_D, 0, sizeof(double) * size_buffer);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;

  // Copy data to device
  debug_info("Memcpy dstpitch %d, srcpitch %d, width %d, height %d",
             sizeof(double) * 4 * ((k - 1) / 4 + 1), sizeof(double) * k,
             sizeof(double) * k, sizeof(double) * n);
  err = cudaMemcpy2D(d_B, sizeof(double) * 4 * ((k - 1) / 4 + 1), B,
                     sizeof(double) * k, sizeof(double) * k, n,
                     cudaMemcpyHostToDevice);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;

  //
  // SNP matrix multiplication
  // The following section multiplies the SNP matrix with a vector of doubles in
  // cutlass
  //

  // Create a problem size struct for matrix multiplication
  cutlass::gemm::GemmCoord problem_size_packed(m, n, k1);

  // Declare Gemm problem
  using CutlassGemm = typename cutlass::gemm::device::Gemm<
      uint8_t, cutlass::layout::RowMajor, cutlass::u4f64_t,
      cutlass::layout::ColumnMajor, double, cutlass::layout::RowMajor, double,
      cutlass::arch::OpClassSimt,
      cutlass::arch::Sm61, // TODO: Check template hierachy if this is required
      cutlass::gemm::GemmShape<32, 32, 16>, // Do not change
      cutlass::gemm::GemmShape<16, 16, 16>, // Might be tuned but must be
                                            // smaller than ThreadblockShape
      cutlass::gemm::GemmShape<1, 1, 1>,    // Do not change
      cutlass::epilogue::thread::LinearCombination<double, 1, double, double>,
      typename cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>, 2,
      1, 1, false,
      cutlass::arch::OpMultiplyAdd // Operator
      >;

  debug_info("Size: %ld", k1 * n);
  // Define CUTLASS GEMM arguments
  typename CutlassGemm::Arguments arguments{
      problem_size_packed, // problem size of matrix multiplication
      {d_M, k1},           // reference to matrix M on device
      {d_B, k1},           // reference to matrix B on device
      {d_C, n},            // reference to matrix C on device
      {d_C, n},            // reference to matrix D on device
      {alpha, beta}        // tuple of alpha and beta
  };

  // Calculate CUTLASS GEMM workspace size
  size_t workspace_size = CutlassGemm::get_workspace_size(arguments);

  // Allocate workspace memory
  cudaMalloc((void **)&d_workspace, max(workspace_size, sizeof(double) * n));

  // Instantiate CUTLASS kernel depending on templates
  CutlassGemm gemm_op;

  // Test if problem can be implemented
  cutlass::Status status = gemm_op.can_implement(arguments);
  if (status != cutlass::Status::kSuccess)
    printf("Can't implement\n");

  // Initialize CUTLASS kernel with arguments and workspace pointer
  status = gemm_op.initialize(arguments, d_workspace);
  if (status != cutlass::Status::kSuccess)
    printf("Error in initialization\n");

  // Launch initialized CUTLASS kernel
  status = gemm_op(); // Actual gemm op
  cudaDeviceSynchronize();
  if (status != cutlass::Status::kSuccess)
    printf("Operation error %d\n", (int) status);

  // Catch all accumulated errors from previous cuda launches
  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;


  //
  // Transpose
  // Switch to column major
  //
  cublas_status = cublasDgeam(cublas_handle, //
                              CUBLAS_OP_T,   // C needs to be transposed
                              CUBLAS_OP_N,   // No-op on B
                              m, // Number of rows of C after transposing
                              n, // Number of columns of C after transposing
                              &alpha, // alpha
                              d_C,    // matrix A
                              n,      // lda
                              &beta,  // beta
                              d_D,    // matrix B
                              m,      // ldb
                              d_D,    // matrix C
                              m);     // ldb
  if (checkError(__func__, __LINE__, cublas_status) != 0)
    return (1);

  //
  // Genotype centering
  // The following section performs genotype centering by substracting op(2 *
  // 1_indiv * f^T ) B from C
  //

  // dgemv performs C <- alpha op(A) x + beta y
  // We calculate f^T B if transa = true or 1_k^T B if transa = false
  // cuBLAS only supports op(A) x, so we calculate B^T f (B^T 1_k resp.), which
  // returns the same as the result is a vector
  // B is of dimension (snps, n) if
  // transa = false and of dimension (indiv,n) if transa = true cuBLAS assumes
  // column-major and d_B is stored in column-major, hence trans = CUBLAS_OP_T
  switch(centered){
    case 0: break;
    case 1: 
    {
      debug_info("Centering: B^f");
      cublas_status = cublasDgemv(cublas_handle, // handle
                                  CUBLAS_OP_T,   // trans
                                  k,             // number of rows of A
                                  n,             // number of cols of A
                                  &alpha_n2,     // alpha
                                  reinterpret_cast<double *>(d_B), // matrix A
                                  ((k - 1) / 4 + 1) * 4,           // lda
                                  transA ? d_unit : d_f,           // vector x
                                  1,                               // incx
                                  &beta,                           // beta
                                  d_workspace,                     // vector y
                                  1);                              // incy
      if (checkError(__func__, __LINE__, cublas_status) != 0)
        return 1;

      debug_info("Centering: C");

      // Now every column i is scaled with alpha_i 1_k if transa = true of alpha_i
      // f if transa = false, where alpha_i = B_i^T f (alpha_i = B_i^T 1_k resp)
      cublas_status =
          cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_DEVICE);
      if (checkError(__func__, __LINE__, cublas_status) != 0)
        return 1;
      for (int i = 0; i < n; i++) {
        cublas_status = cublasDaxpy(cublas_handle,         // handle
                                    m,                     // number of rows
                                    d_workspace + i,       // alpha
                                    transA ? d_f : d_unit, // x
                                    1,                     // incx
                                    d_D + i * m,           //  y
                                    1);                    // incy
        if (checkError(__func__, __LINE__, cublas_status) != 0)
          return 1;
      }
      break;
    }
    default: checkError(__func__, __LINE__, cudaErrorInvalidValue); return 1;  
   }

  //
  // Wrap-up
  //
  // Copy results back to the host
  debug_info("Copy back");
  err = cudaMemcpy(C, d_D, sizeof(double) * m * n, cudaMemcpyDeviceToHost);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;

  debug_info("Free pointers");

  cudaDeviceSynchronize();
  cublas_status = cublasDestroy(cublas_handle);
  if (checkError(__func__, __LINE__, cublas_status) != 0)
    return 1;
    
  cudaFree(d_workspace);
  debug_info("Return");
 
  return 0;
}
