
/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Martin Schlather

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


#ifndef RFU_rfutils_h
#define RFU_rfutils_h 1

#include "errors_messages.h"


//#define SCALAR_RU_H 1
#define SCALAR_BASE 0
#define SCALAR_AVX 1
#define SCALAR_NEARFMA  7 // never change number, see haplogeno.R !!
#define SCALAR_KAHAN 9


typedef // benoetigt
struct solve_storage {
  errorstring_type err_msg;
  InversionMethod method, newMethods[SOLVE_METHODS];
  usr_bool sparse;
  int size, actual_size, actual_pivot;
  int  nsuper;
  Long n_main, n_rhs, n_w2, n_U, n_D, n_w3, n_lnz, n_result;
  //   SICH, n_MM, n_VT, n_ work, n_ nnzlindx,
    
  int 
    *pivot_idx, n_pivot_idx, 
    *iwork, n_iwork, //eigen, svd, LU, spam
    *pivotsparse, n_pivotsparse, *xlnz, n_xlnz, //spam
    *snode, n_snode, *xsuper, n_xsuper, *invp, n_invp,   // spam
    *cols, n_cols, *rows, n_rows, *lindx, n_lindx, // spam
    *xja,  n_xja; // chol, eigen, spam
  double 
  *main, *rhs,// diagonal, general -- FORBIDDEN for further use
    *w2, // eigen, svd, LU, QR, pivot
    *U, // eigen, svd, pivot
    *D, // eigen, svd, cholesky, spam, pivot
    *w3, // spam, QR, svd, eigen
    *lnz, // spam, svd
    *result,  // sqrtPosDefFree
    *to_be_deleted; 
} solve_storage;



#define LINEAR_BASE 0
#define LINEAR_AVX 1

void linearX(double *x, double y, Long len, double *out, Long n);

void solve_DELETE0(solve_storage *x);

utilsoption_type *WhichOptionList(bool local);

#endif
