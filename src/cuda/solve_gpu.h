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

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cusolverSp.h>
#include <cusparse.h>

extern "C" {

/**
 *  \brief C - Wrapper for sparse_solve_init.
 *
 *  Refer to the documentation of sparse_solve_init for details.
 */
void sparse2gpu(double *V, int *I, int *J, long nnz, long m, long max_ncol,
                void **GPU_obj, int *status);

/**
 *  \brief C - Wrapper for sparse_solve_compute.
 *
 *  Refer to the documentation of sparse_solve_compute for details.
 */
void dcsrtrsv_solve_gpu(void *GPU_obj, double *B, int ncol, double *X,
                        int *status);

/**
 *  \brief C - Wrapper for sparse_solve_destroy.
 *
 *  Refer to the documentation of sparse_solve_destroy for details.
 */
void free_sparse_gpu(void **GPU_obj, int *status);
};
// End extern "C"

struct GPU_sparse_storage {
     int *d_cscColPtr;
     int *d_cscRowInd;
     double *d_cscVal;
     int *d_csrRowPtr;
     int *d_csrColInd;
     double *d_csrVal;
     long nnz;
     long m;
     long max_ncol;
     double *d_X;
     double *d_B;
     void   *d_pBuffer_csc;
     void   *d_pBuffer_csr;
     bsrsm2Info_t info_csc;
     bsrsm2Info_t info_csr;
};


// void cholGPU(bool copy, double *matrix, Uint size, double *B, Uint rhs_cols,
//      double *LogDet, double *RESULT);

/**
 *  \brief Initializes storage object with required data on the GPU.
 *  
 *  This function prepares the GPU for iterative solver computation by loading
 *  the required data into the GPU memory.
 *  
 *  \param V A pointer to a vector of matrix values in COO format.
 *  \param I A pointer to a vector of row indices in COO format.
 *  \param J A pointer to a vector of column indices in COO format.
 *  \param nnz The number of non-zero values in the matrix (length of V).
 *  \param m The number of rows and columns in the matrix.
 *  \param maxncol The maximum number of columns in the right-hand side 
 * (RHS) matrix in equation systems.
 *  \param GPUobj A pointer in which the GPU object for iterative solver will be stored.
 *  \param status A pointer to an integer that holds the error code, if any.
 */
void sparse_solve_init(double *V, int *I, int *J, long nnz, long m,
                       long maxncol, void **GPUobj, int *status);

/**
 *  \brief Frees the memory in the GPU object.
 *
 *  This function releases the GPU memory that was allocated for the iterative solver computation.
 *
 *  \param GPU_obj A pointer to the GPU object.
 *  \param status A pointer to an integer that holds the error code, if any.
 */
void sparse_solve_destroy(void **GPU_obj, int *status);

/**
 *  \brief Computes the solution to the equation system defined by the matrix stored in GPU_obj and B.
 *
 *  This function uses the GPU object to solve the equation system. The result is stored in X.
 *
 *  \param GPU_obj A pointer to the GPU object.
 *  \param B A pointer to the RHS matrix of size m x ncol.
 *  \param ncol The number of columns in B and X.
 *  \param X A pointer to the solution matrix of size m x ncol.
 *  \param status A pointer to an integer that holds the error code, if any.
 */
void sparse_solve_compute(void *GPUobj, double *B, int ncol, double *X,
                          int *status);

