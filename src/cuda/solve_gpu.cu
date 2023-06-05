
/*
    Authors
    Alexander Freudenberg
    Copyright (C) 2020 -- 2023 Alexander Freudenberg

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.


    This file provides interfaces to the solving functionality in the cuSOLVER
   and cuSPARSE libraries by NVIDIA. It aims at providing convenience wrappers
   to these functions to alleviate the burden of CUDA memory management and
   other low-level details.

    Functionality:
    -   Dense Solve: Solves an equation system defined by a dense matrix once
   through the Cholesky decomposition and optionally compute the log-determinant
   of the matrix.
    -   Sparse Solve: Prepares a sparse matrix to be solved multiple times for
   use in, e.g., iterative solver algorithms.

    cuSPARSE and cuSOLVER are released as part of the CUDA toolkit under the
   copyright Copyright (C) 2012 -- 2023, NVIDIA Corporation & Affiliates.
*/

#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cusolverSp.h>
#include <cusparse.h>

#include <chrono>
#include <cublasLt.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <vector>

#include <thrust/gather.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>


#include "solve_gpu.h"
#include "cuda_utils.h"


/*
 *
 * C++
 *
 */
#define THREADS_PER_BLOCK 1024

__global__ void logdet_kernel(double* d_matrix, long* d_size, double* d_logdet);

int dense_solve(double *A, // Pointer to the input matrix A in row-major order
                unsigned int input_size, // The dimension of the square matrix A
                double *B, // Pointer to the input matrix B in row-major order
                unsigned int rhs_cols, // The number of columns in the
                                       // right-hand side matrix B
                double *X, // Pointer to the result matrix X in row-major order
                double *logdet,   // Pointer to store the log-determinant of A
                int oversubscribe // Controls whether to use managed memory for
                                  // device oversubscription (1: true, 0: false)
) {

  // Initialize process variables
  unsigned long size = (unsigned long)input_size;
  size_t bufferSize_device = 0, bufferSize_host = 0;
  cudaDataType dataTypeA = CUDA_R_64F;
  int *info = NULL;
  int h_info = 0;

  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  cusolverDnHandle_t handle = NULL;
  cusolverDnParams_t params = NULL;
  cudaStream_t stream = NULL;
  cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
  cudaError_t err = cudaSuccess;

  double *buffer_device = NULL, *buffer_host = NULL, *d_matrix = NULL,
         *d_B = NULL, *d_logdet = NULL;
  long *d_size = NULL;

  // Start time measurement
  clock_t start, start_potrf, start_potrs, start_logdet;
  start = clock();

  // Switch to requested CUDA device
  int device = switchDevice();

  // Initialize handle and stream, calculate buffer size needed for cholesky
  cusolverDnCreate(&handle);
  cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  cusolverDnSetStream(handle, stream);
  cusolverDnCreateParams(&params);

  //
  // Allocate memory on device
  //
  cudaMalloc(&info, sizeof(int));

  if (oversubscribe == 1) {
    int managed_available = 0;
    cudaDeviceGetAttribute(&managed_available, cudaDevAttrManagedMemory,
                           device);
    if (managed_available != 1) {
      printf("This GPU does not suppport managed memory.\n");
    }
    err = cudaMallocManaged((void **)&d_matrix, sizeof(double) * size * size);
  } else if (oversubscribe == 0) {
    err = cudaMalloc((void **)&d_matrix, sizeof(double) * size * size);
  } else {
    checkError(__func__, __LINE__, cudaErrorInvalidValue);
    return 1;
  }
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  err = cudaMalloc((void **)&d_B, sizeof(double) * size * rhs_cols);
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  cudaMemset(info, 0, sizeof(int));

  //
  // Copy data to device
  //
  err = cudaMemcpy(d_matrix, A, sizeof(double) * size * size,
                   cudaMemcpyHostToDevice);
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);
  err = cudaMemcpy(d_B, B, sizeof(double) * size * rhs_cols,
                   cudaMemcpyHostToDevice);
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  // Check for uncaught errors
  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  start_potrf = clock();
  debug_info("Time for initializing: %.3fs", difftime(start_potrf, start));

  //
  // Calculate buffer size
  //
  status = cusolverDnXpotrf_bufferSize(handle, params, uplo, size, dataTypeA,
                                       d_matrix, size, dataTypeA,
                                       &bufferSize_device, &bufferSize_host);
  cudaDeviceSynchronize();
  if (checkError(__func__, __LINE__, status) != 0)
    return (1);

  err = cudaMalloc(&buffer_device, sizeof(double) * bufferSize_device);
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);
  err = cudaMallocHost(&buffer_host, sizeof(double) * bufferSize_host);
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  //
  // Calculate cholesky factorization and store it in d_matrix
  //
  status = cusolverDnXpotrf(handle, params, uplo, size, dataTypeA, d_matrix,
                            size, dataTypeA, buffer_device, bufferSize_device,
                            buffer_host, bufferSize_host, info);
  cudaDeviceSynchronize(); // Synchronize is necessary, otherwise error code
                           // "info" returns nonsense
  if (checkError(__func__, __LINE__, status) != 0)
    return (1);
  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  // Check for errors
  cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  if (0 != h_info) {
    if (h_info > 0)
      printf("Error: Cholesky factorization failed at minor %d \n", h_info);
    if (h_info < 0)
      printf("Error: Wrong parameter in cholesky factorization at %d entry\n",
             h_info);
    err = cudaDeviceReset();
    if (err != cudaSuccess)
      printf("Device reset not successful");
    return (1);
  }

  start_potrs = clock();
  debug_info("Time for potrf: %.3fs", difftime(start_potrs, start_potrf));

  //
  // Calculate x = A\b
  //
  status = cusolverDnXpotrs(handle, params, uplo, size, rhs_cols, dataTypeA,
                            d_matrix, size, dataTypeA, d_B, size, info);
  cudaDeviceSynchronize();
  if (checkError(__func__, __LINE__, status) != 0)
    return (1);
  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  start_logdet = clock();
  debug_info("Time for potrs: %.3fs", difftime(start_logdet, start_potrs));

  //
  // Calculate log determinant of A through its Cholesky root
  //
  if (logdet != NULL) {
    cudaMalloc((void **)&d_logdet, sizeof(double));
    cudaMalloc((void **)&d_size, sizeof(long));
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0)
      return (1);

    cudaMemcpy(d_size, &size, sizeof(long), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, status) != 0)
      return (1);

    logdet_kernel<<<(size - 1) / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK>>>(
        d_matrix, d_size, d_logdet);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, status) != 0)
      return (1);

    err = cudaMemcpy(logdet, d_logdet, sizeof(double), cudaMemcpyDeviceToHost);
    if (checkError(__func__, __LINE__, status) != 0)
      return (1);
    cudaFree(d_size);
    cudaFree(d_logdet);
  }

  debug_info("Time for logdet calculation: %.3fs",
             difftime(clock(), start_logdet));
  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0)
    return (1);

  //
  // Copy  solution from device to vector on host
  //
  err = cudaMemcpy(X, d_B, sizeof(double) * size * rhs_cols,
                   cudaMemcpyDeviceToHost);
  if (checkError(__func__, __LINE__, status) != 0)
    return (1);

  debug_info("Total time: %.3fs", difftime(clock(), start));

  // Free allocated memory
  cudaFree(info);
  cudaFree(buffer_device);
  cudaFree(buffer_host);
  cudaFree(d_matrix);
  cudaFree(d_B);
  cusolverDnDestroy(handle);
  cudaStreamDestroy(stream);
  return 0;
};

void sparse_solve_init(
    double *V,      // Vector of matrix values (COO format)
    int *I,         // Vector of row indices (COO format)
    int *J,         // Vector of column indices (COO format)
    long nnz,       // Number of nonzero values (length of V)
    long m,         // Number of rows and columns of matrix
    long max_ncol,  // Maximum number of columns of RHS in equation systems
    void **GPU_obj, // Pointer in which GPU object for iterative solver will be stored
    int *status     // Holds error code
) {

    // Print compile info
    print_compile_info("cuSPARSE triangular solve interface");


    //
    // Initialize CUDA variables
    //
    cusparseHandle_t handle;
    cusparseMatDescr_t descrL, descrLt;
    bsrsm2Info_t info_csc, info_csr;
    cusparseStatus_t sp_status;
    cudaError_t err;

    // Declare device pointers
    double *d_X               = NULL,
           *d_V               = NULL,
           *d_B               = NULL,
           *d_cscVal          = NULL,
           *d_csrVal          = NULL;
    int    *d_I               = NULL,
           *d_J               = NULL,
           *d_cscColPtr       = NULL,
           *d_csrRowPtr       = NULL,
           *d_cscRowInd       = NULL,
           *d_csrColInd       = NULL;
    void   *d_pBuffer_csc     = NULL,
           *d_pBuffer_csr     = NULL,
           *d_pBuffer_CSC2CSR = NULL;

    int    pBufferSizeInBytes_csc = 0,
           pBufferSizeInBytes_csr = 0,
           structural_zero        = 0;
    size_t pBufferSizeCSC2CSR     = 0;

    // Check CUDA installation
    if (checkCuda() != 0) {
        *status = 1;
        return;
    }

    size_t required_mem = (2 * m * max_ncol + nnz) * sizeof(double) +
                          sizeof(int) * (2 * nnz + (m + 1));
    
    if(checkDevMemory(required_mem) != 0){
        *status = 1;
        return;
    }
    // Allocate memory for device objects
    cudaMalloc((void**)&d_X, sizeof(double) * m * max_ncol);
    cudaMalloc((void**)&d_B, sizeof(double) * m * max_ncol);
    cudaMalloc((void**)&d_I, sizeof(int) * nnz);
    cudaMalloc((void**)&d_J, sizeof(int) * nnz);
    cudaMalloc((void**)&d_V, sizeof(double) * nnz);

    cudaMalloc((void **)&d_cscColPtr, sizeof(int) * (m + 1));
    cudaMalloc((void **)&d_csrRowPtr, sizeof(int) * (m + 1));
    cudaMalloc((void **)&d_csrColInd, sizeof(int) * nnz);
    cudaMalloc((void**)&d_csrVal, sizeof(double) * nnz);

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Copy I J V data to device
    cudaMemcpy(d_I, I, sizeof(int) * nnz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_J, J, sizeof(int) * nnz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_V, V, sizeof(double) * nnz, cudaMemcpyHostToDevice);

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Set up auxiliary stuff
    cusparseCreate(&handle);
    cusparseCreateMatDescr(&descrLt);
    cusparseCreateMatDescr(&descrL);
    cusparseSetMatDiagType(descrLt, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatFillMode(descrLt, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
    cusparseSetMatIndexBase(descrLt, CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatType(descrLt, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);

    //
    // Sort COO data by column -- this is strictly required as CSR routine fails
    // otherwise
    //
    debug_info("Calculating buffer");
    sp_status = cusparseXcoosort_bufferSizeExt(handle, m, m, nnz, d_I, d_J,
                                               &pBufferSizeCSC2CSR);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    debug_info("Allocating");
    cudaMalloc((void**)&d_pBuffer_CSC2CSR, pBufferSizeCSC2CSR);
    
    debug_info("Sequencing");
    thrust::sequence(thrust::device, d_csrColInd, d_csrColInd + nnz);

    debug_info("Sorting");
    sp_status = cusparseXcoosortByColumn(handle, m, m, nnz, d_I, d_J, d_csrColInd, d_pBuffer_CSC2CSR);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    cudaDeviceSynchronize();

    debug_info("Sorting Values");
    thrust::gather(thrust::device, d_csrColInd, d_csrColInd + nnz, d_V, d_csrVal);
    cudaDeviceSynchronize();

    cudaMemcpy(d_V, d_csrVal, sizeof(double) * nnz, cudaMemcpyDeviceToDevice);
    cudaFree(d_pBuffer_CSC2CSR);
    cudaMemset(d_csrColInd, 0, sizeof(int) * nnz);
    cudaMemset(d_csrVal, 0, sizeof(double) * nnz);

    // Initialize CSC data
    sp_status = cusparseXcoo2csr(handle, d_J, nnz, m, d_cscColPtr,
                                 CUSPARSE_INDEX_BASE_ONE);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    d_cscRowInd = d_I;
    d_cscVal    = d_V;
    cudaFree(d_J);

    // Initialize CSR data -- construct from CSC by transposing

    sp_status = cusparseCsr2cscEx2_bufferSize(
        handle, m, m, nnz, d_cscVal, d_cscColPtr, d_cscRowInd, d_csrVal,
        d_csrRowPtr, d_csrColInd, CUDA_R_64F, CUSPARSE_ACTION_NUMERIC,
        CUSPARSE_INDEX_BASE_ONE, CUSPARSE_CSR2CSC_ALG1,
        &pBufferSizeCSC2CSR); // The cusparseCsr2CscAlg_t alg argument is
                              // undocumented in the official docs -- cf
                              // cusparse.h
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    // Allocate buffer memory for CSC2CSR procedure
    required_mem += pBufferSizeCSC2CSR;
    if (checkDevMemory(required_mem) != 0) {
        *status = 1;
        return;
    }
    cudaMalloc((void**)&d_pBuffer_CSC2CSR, pBufferSizeCSC2CSR);
    

#ifdef DEBUG
    int *h_csrRowPtr = NULL, *h_csrColInd = NULL;
    cudaMallocHost((void **)&h_csrRowPtr, sizeof(int) * (m + 1));
    cudaMallocHost((void **)&h_csrColInd, max(sizeof(double) * nnz, pBufferSizeCSC2CSR));
    err = cudaMemcpy(h_csrRowPtr, d_cscColPtr, sizeof(int) * (m + 1),
               cudaMemcpyDeviceToHost);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    err = cudaMemcpy(h_csrColInd, d_cscRowInd, sizeof(int) * nnz,
               cudaMemcpyDeviceToHost);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    debug_info("m %d, nnz %d", m, nnz);
    for (int i = 0; i < min(m + 1, (long) 20); i++)
        printf("%d, ", h_csrRowPtr[i]);
    debug_info(" - cscColPtr\n");
    for (int i = 0; i < min(nnz, (long) 20); i++)
        printf("%d, ", h_csrColInd[i]);
    debug_info(" - RowInd\n");
#endif


    debug_info("m %d nnz %d, Buffer %zu", m, nnz, pBufferSizeCSC2CSR);
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    sp_status = cusparseCsr2cscEx2(
        handle, m, m, nnz, (void *)d_cscVal, d_cscColPtr, d_cscRowInd,
        (void *)d_csrVal, d_csrRowPtr, d_csrColInd, CUDA_R_64F,
        CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ONE,
        CUSPARSE_CSR2CSC_ALG1, d_pBuffer_CSC2CSR);
    cudaDeviceSynchronize();
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        cudaDeviceReset();
        return;
    }
    cudaFree(d_pBuffer_CSC2CSR);

    //
    // Init phase for triangular solve
    //
    sp_status = cusparseCreateBsrsm2Info(&info_csc);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    sp_status = cusparseCreateBsrsm2Info(&info_csr);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    // Get required buffer size of triangular solve
    sp_status = cusparseDbsrsm2_bufferSize(
        handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, m, max_ncol, nnz, descrLt, d_cscVal,
        d_cscColPtr, d_cscRowInd, 1, info_csc, &pBufferSizeInBytes_csc);
    cudaDeviceSynchronize();

    debug_info("Buffer size CSC %d", pBufferSizeInBytes_csc);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    sp_status = cusparseDbsrsm2_bufferSize(
        handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, m, max_ncol, nnz, descrL, d_csrVal,
        d_csrRowPtr, d_csrColInd, 1, info_csr, &pBufferSizeInBytes_csr);
    cudaDeviceSynchronize();

    debug_info("Buffer size CSR %d", pBufferSizeInBytes_csr);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    // Allocate buffer memory for triangular solve
    required_mem += 2*pBufferSizeInBytes_csc;
    if (checkDevMemory(required_mem) != 0) {
        *status = 1;
        return;
    }
    cudaMalloc((void**)&d_pBuffer_csc, pBufferSizeInBytes_csc);
    cudaMalloc((void**)&d_pBuffer_csr, pBufferSizeInBytes_csr);

    // Perform analysis phase of triangular solve 
    sp_status = cusparseDbsrsm2_analysis(
        handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, m, max_ncol, nnz, descrLt, d_cscVal,
        d_cscColPtr, d_cscRowInd, 1, info_csc, CUSPARSE_SOLVE_POLICY_NO_LEVEL,
        d_pBuffer_csc);
    cudaDeviceSynchronize();

    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    // Check for solvability in Cholesky root L
    sp_status = cusparseXbsrsm2_zeroPivot(handle, info_csc, &structural_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == sp_status) {
        printf("Structural zero in Cholesky root CSC: L(%d,%d) is missing\n",
               structural_zero, structural_zero);
        *status = 1;
        return;
    }

    sp_status = cusparseDbsrsm2_analysis(
        handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, m, max_ncol, nnz, descrL, d_csrVal,
        d_csrRowPtr, d_csrColInd, 1, info_csr, CUSPARSE_SOLVE_POLICY_NO_LEVEL,
        d_pBuffer_csr);
    cudaDeviceSynchronize();

    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    // Check for solvability in Cholesky root L 
    sp_status = cusparseXbsrsm2_zeroPivot(handle, info_csr, &structural_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == sp_status) {
        printf("Structural zero in Cholesky root CSR: L(%d,%d) is missing\n",
               structural_zero, structural_zero);
        *status = 1;
        return;
    }

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    //
    // Initialize GPU_sparse_storage object
    //
    struct GPU_sparse_storage *GPU_storage_obj =
        (struct GPU_sparse_storage *)malloc(sizeof(struct GPU_sparse_storage));
    GPU_storage_obj->d_cscColPtr   = d_cscColPtr;
    GPU_storage_obj->d_cscRowInd   = d_cscRowInd;
    GPU_storage_obj->d_cscVal      = d_cscVal;
    GPU_storage_obj->d_csrRowPtr   = d_csrRowPtr;
    GPU_storage_obj->d_csrColInd   = d_csrColInd;
    GPU_storage_obj->d_csrVal      = d_csrVal;
    GPU_storage_obj->nnz           = nnz;
    GPU_storage_obj->m             = m;
    GPU_storage_obj->max_ncol      = max_ncol;
    GPU_storage_obj->d_X           = d_X;
    GPU_storage_obj->d_B           = d_B;
    GPU_storage_obj->d_pBuffer_csc = d_pBuffer_csc;
    GPU_storage_obj->d_pBuffer_csr = d_pBuffer_csr;
    GPU_storage_obj->info_csc      = info_csc;
    GPU_storage_obj->info_csr      = info_csr;

    debug_info("Pointer pBuffer_csc %d", d_pBuffer_csc);

    // Set pointer to initialized object
    *GPU_obj = (void *)GPU_storage_obj;

    sp_status = cusparseDestroy(handle);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
};
void sparse_solve_destroy(void **GPU_obj, // Pointer to GPU object
                          int *status     // Holds error code
) {
    cudaError_t err;
    bsrsm2Info_t info_csc, info_csr;

    // Check CUDA installation
    if (checkCuda() != 0) {
        *status = 1;
        return;
    }
    // Check if valid pointer
    if (*GPU_obj == NULL) {
        return 1;
    }

    cudaDeviceSynchronize();

    // Get GPU storage object
    struct GPU_sparse_storage *GPU_storage_obj =
        (struct GPU_sparse_storage *)(*GPU_obj);

    debug_info("Pointer value %d", GPU_storage_obj);
    if((GPU_obj == NULL) || (GPU_storage_obj == NULL)){
        checkError(__func__, __LINE__, cudaErrorIllegalAddress);
        *status = 1;
        return;
    }

    // Declare device pointers
    double *d_X           = GPU_storage_obj->d_X,
           *d_B           = GPU_storage_obj->d_B,
           *d_cscVal      = GPU_storage_obj->d_cscVal,
           *d_csrVal      = GPU_storage_obj->d_csrVal;
    int    *d_cscColPtr   = GPU_storage_obj->d_cscColPtr,
           *d_cscRowInd   = GPU_storage_obj->d_cscRowInd,
           *d_csrRowPtr   = GPU_storage_obj->d_csrRowPtr,
           *d_csrColInd   = GPU_storage_obj->d_csrColInd;
    void   *d_pBuffer_csc = GPU_storage_obj->d_pBuffer_csc,
           *d_pBuffer_csr = GPU_storage_obj->d_pBuffer_csr;

    debug_info("Pointer pBuffer_csc %d", d_pBuffer_csc);

    info_csc = GPU_storage_obj->info_csc;
    info_csr = GPU_storage_obj->info_csr;

    err = cudaFree(d_X);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_B);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_cscVal);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_cscColPtr);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_cscRowInd);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_pBuffer_csc);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_csrVal);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_csrRowPtr);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_csrColInd);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_pBuffer_csr);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    cusparseDestroyBsrsm2Info(info_csc);
    cusparseDestroyBsrsm2Info(info_csr);
    free(GPU_storage_obj);
    GPU_obj = NULL;

    cudaDeviceReset();

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
};

void sparse_solve_compute(void *GPU_obj, // Pointer to GPU object
                          double *B,  // Pointer to RHS matrix of size m x ncol
                          int ncol,   // Number of columns of B and X
                          double *X,  // Solution matrix of size size m x ncol
                          int *status // Holds error code
) {

    //
    // Initialize CUDA variables
    //
    cusparseHandle_t handle;
    cusparseMatDescr_t descrL, descrLt;
    bsrsm2Info_t info_csc, info_csr;
    cusparseStatus_t sp_status;
    cudaError_t err;

    // Check CUDA installation
    if (checkCuda() != 0) {
        *status = 1;
        return;
    }

    // Get GPU storage object
    struct GPU_sparse_storage *GPU_storage_obj =
        (struct GPU_sparse_storage *)GPU_obj;

    // Get problem dimensions
    long m = GPU_storage_obj->m, nnz = GPU_storage_obj->nnz;
    int max_ncol = GPU_storage_obj->max_ncol;

    if (ncol > max_ncol) {
        printf("Sparse solve interface has been initialized with %d columns, "
               "but %d columns are requested by the calculation function.\n",
               max_ncol, ncol);
        *status = 1;
        return;
    }
    debug_info("Start calc");
    debug_info("%d %d %d", m, nnz, max_ncol);
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Declare device pointers
    double *d_X = GPU_storage_obj->d_X, *d_B = GPU_storage_obj->d_B,
           *d_cscVal = GPU_storage_obj->d_cscVal,
           *d_csrVal = GPU_storage_obj->d_csrVal;
    int *d_cscColPtr = GPU_storage_obj->d_cscColPtr,
        *d_cscRowInd = GPU_storage_obj->d_cscRowInd,
        *d_csrRowPtr = GPU_storage_obj->d_csrRowPtr,
        *d_csrColInd = GPU_storage_obj->d_csrColInd;
    void *d_pBuffer_csc = GPU_storage_obj->d_pBuffer_csc,
         *d_pBuffer_csr = GPU_storage_obj->d_pBuffer_csr;

    int numerical_zero = 0;
    const double alpha = 1.0;

    // Get CUDA auxiliary variables
    info_csc = GPU_storage_obj->info_csc;
    info_csr = GPU_storage_obj->info_csr;
    cusparseCreate(&handle);

    cusparseCreateMatDescr(&descrLt);
    cusparseCreateMatDescr(&descrL);
    cusparseSetMatDiagType(descrLt, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatFillMode(descrLt, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
    cusparseSetMatIndexBase(descrLt, CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatType(descrLt, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Reset memory for X and B on device
    err = cudaMemset(d_B, 0.0, sizeof(double) * m * max_ncol);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaMemset(d_X, 0.0, sizeof(double) * m * max_ncol);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Copy data to device
    err = cudaMemcpy(d_B, B, sizeof(double) * m * ncol, cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    //
    // Solving equation system on device - forward substitution
    //
    auto start = clock();

    sp_status = cusparseDbsrsm2_solve(
        handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, m, ncol, nnz, &alpha, descrL,
        d_csrVal, d_csrRowPtr, d_csrColInd, 1, info_csr, d_B, m, d_X, m,
        CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_pBuffer_csr);

    cudaDeviceSynchronize();
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    sp_status = cusparseXbsrsm2_zeroPivot(handle, info_csr, &numerical_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == sp_status) {
        printf("Numerical zero during solving: L(%d,%d) is zero\n",
               numerical_zero, numerical_zero);
        *status = 1;
        return;
    }

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Backward substitution to get result of L L^T X = B
    sp_status = cusparseDbsrsm2_solve(
        handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE,
        CUSPARSE_OPERATION_NON_TRANSPOSE, m, ncol, nnz, &alpha, descrLt,
        d_cscVal, d_cscColPtr, d_cscRowInd, 1, info_csc, d_X, m, d_B, m,
        CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_pBuffer_csc);

    cudaDeviceSynchronize();
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    sp_status = cusparseXbsrsm2_zeroPivot(handle, info_csc, &numerical_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == sp_status) {
        printf("Numerical zero during solving: L(%d,%d) is zero\n",
               numerical_zero, numerical_zero);
        *status = 1;
        return;
    }

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    debug_info("Time: %.3f", (double)(clock() - start) / CLOCKS_PER_SEC);

    // Copy results back to device
    err = cudaMemcpy(X, d_B, sizeof(double) * m * ncol, cudaMemcpyDeviceToHost);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    debug_info("Returning");
    cusparseDestroyMatDescr(descrLt);
    cusparseDestroyMatDescr(descrL);
    cusparseDestroy(handle);
};


__global__ void logdet_kernel(double* d_matrix, long* d_size, double* d_logdet)
{
    /* This CUDA kernel calculates the logdeterminant of a matrix by determining the trace of its cholesky decomposition
    Input:
        d_matrix pointer to matrix
        d_size size of matrix
    Output:
        d_logdet pointer to logdeterminant on device
    */
    __shared__ double logdet_loc;
    __shared__ double submatrix[THREADS_PER_BLOCK];
    logdet_loc = 0.0;
    *d_logdet = 0.0;
    long idx = blockDim.x * blockIdx.x + threadIdx.x,
         thread = threadIdx.x;
    if (idx < *d_size) {
        submatrix[thread] = d_matrix[idx * (*d_size + 1)];
    }
    __syncthreads();
    atomicAdd(&logdet_loc, idx >= *d_size ? 0 : 2.0*(log(submatrix[thread])));

    __syncthreads();
    if (threadIdx.x == 0) {
        atomicAdd(d_logdet, logdet_loc);
    };
};


extern "C" {
/*
 * This section contains a series of 'extern "C"' functions which are wrappers 
 * for the original C++ functions defined below. 
 * 
 * The purpose of these functions is to provide an interface that allows 
 * other languages, like Fortan, R and Julia, to call the C++ functions, as they can only 
 * directly interact with C-style functions.
 * 
 * For a detailed explanation and documentation of each function's behavior, 
 * parameters, and return values, refer to the comments provided in the 
 * original C++ functions. Keep in mind that these wrapper functions maintain 
 * the same functionalities and behaviors as their original C++ counterparts.
 */

int potrs_solve(double *A, unsigned int input_size, double *B,
                unsigned int rhs_cols, double *X, double *logdet,
                int oversubscribe) {
    return potrs_solve(A, input_size, B, rhs_cols, X, logdet, oversubscribe);
};

void sparse2gpu(double *V, int *I, int *J, long nnz, long m, long max_ncol,
                void **GPU_obj, int *status) {
    sparse_solve_init(V, I, J, nnz, m, max_ncol, GPU_obj, status);
};

void dcsrtrsv_solve_gpu(void *GPU_obj, double *B, int ncol, double *X,
                        int *status) {
    sparse_solve_compute(GPU_obj, B, ncol, X, status);
};

void free_sparse_gpu(void **GPU_obj, int *status) {
    sparse_solve_destroy(GPU_obj, status);
};

void potrs_solve_gpu(double *A, unsigned int input_size, double *B,
                     unsigned int rhs_cols, double *X, double *logdet,
                     int oversubscribe, int *status) {
    *status = dense_solve(A, input_size, B, rhs_cols, X, logdet, oversubscribe);
};

}
