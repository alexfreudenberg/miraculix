/*
   Authors
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

/**
 *  \brief C - Wrapper for dense_solve.
 *
 *  Refer to the documentation of dense_solve for details.
 */
void potrs_solve_gpu(double *A, unsigned int input_size, double *B,
                     unsigned int rhs_cols, double *X, double *logdet,
                     int oversubscribe, int *status);
};
// End extern "C"

struct GPU_sparse_storage {
  int *d_I;
  int *d_J;
  double *d_V;
  long nnz;
  long m;
  long ncol;
  int is_lower;
  int *d_X;
  int *d_B;
  cusparseConstSpMatDescr_t *matA;
  cusparseConstDnMatDescr_t *matB;
  cusparseSpSMDescr_t *spsmDescr_noop;
  cusparseSpSMDescr_t *spsmDescr_trans
};

/**
 * Solves an equation system defined by a symmetric, positive-definite matrix A and the right-hand side B.
 * 
 * @param A         Pointer to the input matrix A in row-major order.
 * @param input_size The dimension of the square matrix A.
 * @param B         Pointer to the input matrix B in row-major order.
 * @param rhs_cols  The number of columns in the right-hand side matrix B.
 * @param X         Pointer to the result matrix X in row-major order.
 * @param logdet    Pointer to store the log-determinant of A.
 * @param oversubscribe Controls whether to use managed memory for device oversubscription (1: true, 0: false).
 * 
 * @return          Returns an integer indicating the success (0) or failure (-1) of the solve operation.
 * 
 * @note            The pointer to X can be equal to the one for B. 
 * @note            This function uses the potrs and potrf routines in the cuSOLVER library.
 */
int dense_solve(double* A, unsigned int input_size, double* B, unsigned int rhs_cols,
    double* X, double *logdet, int oversubscribe);

/**
 *  \brief Initializes storage object with required data on the GPU for solving an equation system defined by a sparse symmetric, positive-definite matrix A. 
 *  
 *  This function prepares the GPU for iterative solver computation by loading
 *  the sparse matrix into the GPU memory.
 *  
 *  \param V A pointer to the vector of values of A in COO format.
 *  \param I A pointer to the vector of row indices of A in COO format.
 *  \param J A pointer to the vector of column indices of A in COO format.
 *  \param nnz The number of non-zero values in the matrix (length of V).
 *  \param m The number of rows and columns in the matrix.
 *  \param maxncol The maximum number of columns in the right-hand side 
 * (RHS) matrix in equation systems.
 *  \param GPUobj A pointer in which the GPU object for iterative solver will be stored.
 *  \param status A pointer to an integer that holds the error code, if any.
 * 
 *  \note The compile message can be switched off by setting the environment variable PRINT_LEVEL to -1.
 *  \note This function uses the sparse matrix routines in the cuSPARSE library, in particular the analysis routine of the bsrsm2 collection.
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
 * 
 *  \note This function uses the solve function in the bsrsm2 function collection in the cuSPARSE library. 
 */
void sparse_solve_compute(void *GPUobj, double *B, int ncol, double *X,
                          int *status);

