
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
    long *I,        // Vector of row indices (COO format)
    long *J,        // Vector of column indices (COO format)
    long nnz,       // Number of nonzero values (length of V)
    long m,         // Number of rows and columns of matrix
    long ncol,      // Maximum number of columns of RHS in equation systems
    int is_lower,   // If matrix is lower triangular or upper triangular
    void **GPU_obj, // Pointer in which GPU object for iterative solver will be
                    // stored
    int *status     // Holds error code
) {

  // Print compile info
  print_compile_info("cuSPARSE triangular solve interface");
  debug_info("Init params: m %d nnz %d ncol %d is_lower %d", m, nnz, ncol,
             is_lower);

  //
  // Initialize CUDA variables
  //
  cusparseHandle_t handle;
  cudaError_t err;
  cusparseStatus_t sp_status;
  cusparseFillMode_t fill_mode =
      is_lower ? CUSPARSE_FILL_MODE_LOWER : CUSPARSE_FILL_MODE_UPPER;
  cusparseDiagType_t diag_type = CUSPARSE_DIAG_TYPE_NON_UNIT;
  cusparseConstSpMatDescr_t *matA = NULL;
  cusparseConstDnMatDescr_t *matB = NULL;
  cusparseDnMatDescr_t *matC = NULL;
  cusparseSpSMDescr_t *spsmDescr_noop = NULL, *spsmDescr_trans = NULL;

  spsmDescr_noop = (cusparseSpSMDescr_t *)malloc(sizeof(cusparseSpSMDescr_t));
  spsmDescr_trans = (cusparseSpSMDescr_t *)malloc(sizeof(cusparseSpSMDescr_t));
  matA = (cusparseConstSpMatDescr_t *)malloc(sizeof(cusparseSpMatDescr_t));
  matB = (cusparseConstDnMatDescr_t *)malloc(sizeof(cusparseDnMatDescr_t));
  matC = (cusparseDnMatDescr_t *)malloc(sizeof(cusparseDnMatDescr_t));

  // Declare device pointers
  double *d_X = NULL, *d_V = NULL, *d_B = NULL;
  long *d_I = NULL, *d_J = NULL;
  void *d_buffer_noop = NULL, *d_buffer_trans = NULL;

  size_t bufferSize_noop = 0, bufferSize_trans = 0;
  const double alpha = 1.0;

  // Check CUDA installation
  if (checkCuda() != 0) {
    *status = 1;
    return;
  }

  size_t required_mem = (2 * m * ncol + nnz) * sizeof(double) +
                        sizeof(long) * (2 * nnz + (m + 1));

  if (checkDevMemory(required_mem) != 0) {
    *status = 1;
    return;
  }
  // Allocate memory for device objects
  cudaMalloc((void **)&d_X, sizeof(double) * m * ncol);
  cudaMalloc((void **)&d_B, sizeof(double) * m * ncol);
  cudaMalloc((void **)&d_I, sizeof(long) * nnz);
  cudaMalloc((void **)&d_J, sizeof(long) * nnz);
  cudaMalloc((void **)&d_V, sizeof(double) * nnz);

  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0) {
    *status = 1;
    return;
  }

  // Copy I J V data to device
  debug_info("Copying");
  cudaMemcpy(d_I, I, sizeof(long) * nnz, cudaMemcpyHostToDevice);
  cudaMemcpy(d_J, J, sizeof(long) * nnz, cudaMemcpyHostToDevice);
  cudaMemcpy(d_V, V, sizeof(double) * nnz, cudaMemcpyHostToDevice);

  err = cudaGetLastError();
  if (checkError(__func__, __LINE__, err) != 0) {
    *status = 1;
    return;
  }

  // Create handle
  cusparseCreate(&handle);

  //
  // Set matrix description holding all the information required for the solve
  // routine
  //
  debug_info("Creating mat descriptions");
  sp_status =
      cusparseCreateConstCoo(matA, // sparse matrix description to be filled
                             m,    // Number of rows
                             m,    // Number of columns
                             nnz,  // Number of non-zeros
                             d_I,  // Vector of rows
                             d_J,  // Vector of columns
                             d_V,  // Vector of values
                             CUSPARSE_INDEX_64I,      // Index data type
                             CUSPARSE_INDEX_BASE_ONE, // Index base
                             CUDA_R_64F);             // Value type
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  // Set matrix fill mode - lower or upper depending on input
  sp_status = cusparseSpMatSetAttribute((cusparseSpMatDescr_t)*matA,
                                        CUSPARSE_SPMAT_FILL_MODE, &fill_mode,
                                        sizeof(cusparseFillMode_t));
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  // Set matrix diag type
  sp_status = cusparseSpMatSetAttribute((cusparseSpMatDescr_t)*matA,
                                        CUSPARSE_SPMAT_DIAG_TYPE, &diag_type,
                                        sizeof(cusparseDiagType_t));
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  // Set matrix B descriptor
  sp_status = cusparseCreateConstDnMat(matB, m, ncol, m, d_B, CUDA_R_64F,
                                       CUSPARSE_ORDER_COL);
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  // Set matrix C descriptor
  sp_status = cusparseCreateDnMat(matC, m, ncol, m, d_X, CUDA_R_64F,
                                  CUSPARSE_ORDER_COL);
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }

  //
  // Set up triangular solve descriptors which include storage information
  //
  sp_status = cusparseSpSM_createDescr(spsmDescr_noop);
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  sp_status = cusparseSpSM_createDescr(spsmDescr_trans);
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }

  // Get required buffer sizes for triangular solve
  sp_status = cusparseSpSM_bufferSize(
      handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
      CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, *matA, *matB, *matC, CUDA_R_64F,
      CUSPARSE_SPSM_ALG_DEFAULT, *spsmDescr_noop, &bufferSize_noop);
  cudaDeviceSynchronize();
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }

  sp_status = cusparseSpSM_bufferSize(
      handle, CUSPARSE_OPERATION_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha, *matA, *matB, *matC, CUDA_R_64F, CUSPARSE_SPSM_ALG_DEFAULT,
      *spsmDescr_trans, &bufferSize_trans);
  cudaDeviceSynchronize();
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }

  if (checkError(__func__, __LINE__, cudaGetLastError()) != 0) {
    *status = 1;
    return;
  }
  //
  // Allocate buffer memory
  //
  required_mem += bufferSize_noop + bufferSize_trans;
  if (checkDevMemory(required_mem) != 0) {
    *status = 1;
    return;
  }
  cudaMalloc((void **)&d_buffer_noop, bufferSize_noop);
  cudaMalloc((void **)&d_buffer_trans, bufferSize_trans);

  if (checkError(__func__, __LINE__, cudaGetLastError()) != 0) {
    *status = 1;
    return;
  }

  //
  // Analysis phase for triangular solve
  //
  sp_status = cusparseSpSM_analysis(
      handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
      CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, *matA, *matB, *matC, CUDA_R_64F,
      CUSPARSE_SPSM_ALG_DEFAULT, *spsmDescr_noop, d_buffer_noop);
  cudaDeviceSynchronize();
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  if (checkError(__func__, __LINE__, cudaGetLastError()) != 0) {
    *status = 1;
    return;
  }

  sp_status = cusparseSpSM_analysis(
      handle, CUSPARSE_OPERATION_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
      &alpha, *matA, *matB, *matC, CUDA_R_64F, CUSPARSE_SPSM_ALG_DEFAULT,
      *spsmDescr_trans, d_buffer_trans);
  cudaDeviceSynchronize();
  if (checkError(__func__, __LINE__, sp_status) != 0) {
    *status = 1;
    return;
  }
  if (checkError(__func__, __LINE__, cudaGetLastError()) != 0) {
    *status = 1;
    return;
  }

  //
  // Initialize GPU_sparse_storage object
  //
  struct GPU_sparse_storage *GPU_storage_obj =
      (struct GPU_sparse_storage *)malloc(sizeof(struct GPU_sparse_storage));
  GPU_storage_obj->d_I = d_I;
  GPU_storage_obj->d_J = d_J;
  GPU_storage_obj->d_V = d_V;
  GPU_storage_obj->nnz = nnz;
  GPU_storage_obj->m = m;
  GPU_storage_obj->ncol = ncol;
  GPU_storage_obj->is_lower = is_lower;
  GPU_storage_obj->d_X = d_X;
  GPU_storage_obj->d_B = d_B;
  GPU_storage_obj->matA = matA;
  GPU_storage_obj->matB = matB;
  GPU_storage_obj->matC = matC;
  GPU_storage_obj->spsmDescr_noop = spsmDescr_noop;
  GPU_storage_obj->spsmDescr_trans = spsmDescr_trans;
  GPU_storage_obj->d_buffer_noop = d_buffer_noop;
  GPU_storage_obj->d_buffer_trans = d_buffer_trans;

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
    cusparseStatus_t sp_status;
    cusparseConstSpMatDescr_t *matA            = NULL;
    cusparseConstDnMatDescr_t *matB            = NULL;
    cusparseDnMatDescr_t      *matC            = NULL;
    cusparseSpSMDescr_t       *spsmDescr_noop  = NULL,
                              *spsmDescr_trans = NULL;

    double *d_X, *d_B, *d_V;
    long *d_I, *d_J;
    void *d_buffer_noop, *d_buffer_trans; 

    // Check CUDA installation
    if (checkCuda() != 0) {
        *status = 1;
        return;
    }
    cudaDeviceSynchronize();

    // Get GPU storage object
    struct GPU_sparse_storage *GPU_storage_obj =
        (struct GPU_sparse_storage *)(*GPU_obj);

    if((GPU_obj == NULL) || (GPU_storage_obj == NULL)){
        checkError(__func__, __LINE__, cudaErrorIllegalAddress);
        *status = 1;
        return;
    }

    // Declare device pointers
    d_X             = GPU_storage_obj->d_X;
    d_B             = GPU_storage_obj->d_B;
    d_V             = GPU_storage_obj->d_V;
    d_I             = GPU_storage_obj->d_I;
    d_J             = GPU_storage_obj->d_J;
    matA            = GPU_storage_obj->matA;
    matB            = GPU_storage_obj->matB;
    matC            = GPU_storage_obj->matC;
    spsmDescr_noop  = GPU_storage_obj->spsmDescr_noop;
    spsmDescr_trans = GPU_storage_obj->spsmDescr_trans;
    d_buffer_noop   = GPU_storage_obj->d_buffer_noop;
    d_buffer_trans  = GPU_storage_obj->d_buffer_trans;

    // Free device memory
    err = cudaFree(d_buffer_noop);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_buffer_trans);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
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
    err = cudaFree(d_V);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_I);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaFree(d_J);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    // Destory cuSPARSE descriptors
    sp_status = cusparseDestroySpMat(*matA);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    sp_status = cusparseDestroyDnMat(*matB);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    sp_status = cusparseDestroyDnMat(*matC);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    sp_status = cusparseSpSM_destroyDescr(*spsmDescr_noop);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    sp_status = cusparseSpSM_destroyDescr(*spsmDescr_trans);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }


    free(GPU_storage_obj);
    *GPU_obj = NULL;

    cudaDeviceReset();

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
};

void sparse_solve_compute(
    void *GPU_obj, // Pointer to GPU object
    char transA,   // If A should be transposed ('T', 't') or not ('N', 'n's)
    double *B,     // Pointer to RHS matrix of size m x ncol
    long ncol,     // Number of columns of B and X
    double *X,     // Solution matrix of size size m x ncol
    int *status    // Holds error code
) {

    //
    // Initialize CUDA variables
    //
    cusparseHandle_t handle;
    cusparseStatus_t sp_status;
    cudaError_t err;

    cusparseConstSpMatDescr_t *matA;
    cusparseConstDnMatDescr_t *matB;
    cusparseDnMatDescr_t *matC;
    cusparseSpSMDescr_t *spsmDescr_noop;
    cusparseSpSMDescr_t *spsmDescr_trans;

    double *d_X, *d_B;
    bool trans;
    const double alpha = 1.0;

    switch (transA) {
    case 'T':
        trans = true;
        break;
    case 't':
        trans = true;
        break;
    case 'N':
        trans = false;
        break;
    case 'n':
        trans = false;
        break;
    default:
        checkError(__func__, __LINE__, cudaErrorInvalidValue);
        debug_info("transA: %c", transA);
        *status = 1;
        return;
    }

    // Check CUDA installation
    if (checkCuda() != 0) {
        *status = 1;
        return;
    }

    // Get GPU storage object
    struct GPU_sparse_storage *GPU_storage_obj =
        (struct GPU_sparse_storage *)GPU_obj;

    // Get problem dimensions
    long m   = GPU_storage_obj->m,
         nnz = GPU_storage_obj->nnz;

    // Validate correct problem dimension
    if (ncol != GPU_storage_obj->ncol) {
        printf("Sparse solve interface has been initialized with %d columns, "
               "but %d columns are requested by the compute function.\n",
               GPU_storage_obj->ncol, ncol);
        *status = 1;
        return;
    }
    debug_info("Compute params: m %d nnz %d ncol %d", m, nnz, ncol);
    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Get device pointers
    d_B             = GPU_storage_obj->d_B;
    d_X             = GPU_storage_obj->d_X;
    matA            = GPU_storage_obj->matA;
    matB            = GPU_storage_obj->matB;
    matC            = GPU_storage_obj->matC;
    spsmDescr_noop  = GPU_storage_obj->spsmDescr_noop;
    spsmDescr_trans = GPU_storage_obj->spsmDescr_trans;


    // Set-up cuSPARSE handle
    cusparseCreate(&handle);

    // Reset memory for X and B on device
    err = cudaMemset(d_B, 0.0, sizeof(double) * m * ncol);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    err = cudaMemset(d_X, 0.0, sizeof(double) * m * ncol);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    // Copy data to device and set values
    debug_info("Copying");
    err = cudaMemcpy(d_B, B, sizeof(double) * m * ncol, cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    sp_status = cusparseDnMatSetValues((cusparseDnMatDescr_t) *matB, d_B);
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    //
    // Solving triangular equation system on device
    //
    auto start = clock();

    sp_status = cusparseSpSM_solve(
        handle, // cuSPARSE handle
        trans ? CUSPARSE_OPERATION_TRANSPOSE
              : CUSPARSE_OPERATION_NON_TRANSPOSE, // Op(A)
        CUSPARSE_OPERATION_NON_TRANSPOSE,         // Op(B)
        &alpha,                    // alpha value for equation system
        *matA,                     // matA descriptor
        *matB,                     // matB descriptor
        *matC,                     // matC descriptor
        CUDA_R_64F,                // value type
        CUSPARSE_SPSM_ALG_DEFAULT, // cuSPARSE algorithm
        trans ? *spsmDescr_trans : *spsmDescr_noop // spsm storage object
    );
    cudaDeviceSynchronize();
    if (checkError(__func__, __LINE__, sp_status) != 0) {
        *status = 1;
        return;
    }

    err = cudaGetLastError();
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }
    debug_info("Time: %.3f", (double)(clock() - start) / CLOCKS_PER_SEC);

    // Copy results back to host
    err = cudaMemcpy(X, d_X, sizeof(double) * m * ncol, cudaMemcpyDeviceToHost);
    if (checkError(__func__, __LINE__, err) != 0) {
        *status = 1;
        return;
    }

    debug_info("Returning");
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

void sparse2gpu(double *V, long *I, long *J, long nnz, long m, long ncol,
                int is_lower, void **GPU_obj, int *status) {
    sparse_solve_init(V, I, J, nnz, m, ncol, is_lower, GPU_obj, status);
};

void dcsrtrsv_solve_gpu(void *GPU_obj, char transA, double *B, long ncol, double *X,
                        int *status) {
    sparse_solve_compute(GPU_obj, transA, B, ncol, X, status);
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
