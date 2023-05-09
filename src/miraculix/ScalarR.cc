
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


#include "Basic_miraculix.h"

#if defined compatibility_to_R_h
#include "xport_import.h"
#include "Scalar.h"
#include "MXinfo.h"
#include "options.h"
#include "miraculix.h"


SEXP crossprodInt(SEXP X, SEXP Y, SEXP mode) {
  KEY_type *KT = KEYT_M();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  int n;
  Uint nrow,
    len,
    lenY,
    ncol;
  if (isMatrix(X)) {
    nrow = ncols(X);
    len = nrows(X);
  } else {
    nrow = 1;
    len = LENGTH(X);
  }
  if (isMatrix(Y)) {
    ncol = ncols(Y);
    lenY = nrows(Y);
  } else { 
    ncol = 1;
    lenY = LENGTH(Y);
  }
  if (lenY != len) ERR0("sizes of 'x' and 'y' do not match");
  if (LENGTH(mode) == 0) n = SCALAR_INT_DEFAULT;
  else {
    n = INTEGER(mode)[0];
    if (n < 0) n =  SCALAR_INT_DEFAULT;
  }

  
  SEXP Ans; 
  PROTECT(Ans = allocMatrix(INTSXP, nrow, ncol));
  Uint *x, *y,
    *ans = (Uint*) INTEGER(Ans);
  if (TYPEOF(X) == INTSXP) x = (Uint*) INTEGER(X); else x=(Uint*) LOGICAL(X);
  if (TYPEOF(Y) == INTSXP) y = (Uint *)INTEGER(Y); else y=(Uint*) LOGICAL(Y);

  crossprod_Int(x, y, nrow, ncol, len, cores,  ans);
  
  UNPROTECT(1);
  return Ans;
}



#endif
