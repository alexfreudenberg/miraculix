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

#define Uint unsigned int
#define PADDIM 4L
#define THREADS_PER_BLOCK 1024 // 2048 / 32
#define BLOCKS 1024

extern "C" {
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

void sparse2gpu(
     double *V,     // Vector of matrix values (COO format)
     int *I,        // Vector of row indices (COO format)
     int *J,        // Vector of column indices (COO format)
     long nnz,       // Number of nonzero values (length of V)
     long m,         // Number of rows and columns of matrix
     long max_ncol,  // Maximum number of columns of RHS in equation systems
     void **GPU_obj, // Pointer in which GPU object for iterative solver will be
                    // stored
     int *status
);
void dcsrtrsv_solve(void *GPU_obj, // Pointer to GPU object
                   double *B,     // Pointer to RHS matrix of size m x ncol
                   int ncol,      // Number of columns of B and X
                   double *X,     // Solution matrix of size size m x ncol
                   int *status
);
void freegpu_sparse(void *GPU_obj, int *status);

};