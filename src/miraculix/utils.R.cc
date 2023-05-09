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


#include "Basic_RFUlocal.h"

#if defined compatibility_to_R_h
#include "compatibility.lapack.h"
#include "RandomFieldsUtils.h"
#include "kleinkram.h"
#include "zzz_RFU.h"
#include "utils.h"
#include "xport_import_RFU.h"
#include "extern_RFU.h"


#ifdef __cplusplus
extern "C" int *ToIntI(SEXP X, bool *create, bool round);
#endif

int *ToIntI(SEXP X, bool *create, bool round) {
  KEY_type *KT = KEYT();
  if (TYPEOF(X) == INTSXP) {
    *create = false;
    return INTEGER(X);
  }
  if (TYPEOF(X) == LGLSXP) {
    *create = false;
    return LOGICAL(X);
  }
  int len = length(X);
  // TO DO !!  
  //  if (len > 100 || PL > 1)
  //    HELPINFO("Better use 'integer' as storage mode (for one of the arguments).");
  int *y;
  if (*create || KT->ToIntN < len) {
    y = (int *) MALLOC(sizeof(*y) * len);    
    if (y == NULL) ERR1("not enough memory for an %d vector of integers", len);
    if (!*create) {
      FREE(KT->ToIntDummy);
      KT->ToIntDummy = y;
      KT->ToIntN = len;
    }
  } else y = KT->ToIntDummy;
  double *x = (double *) REAL(X);
  if (round) for (int i=0; i<len; i++) y[i] = (int) ROUND(x[i]);
  else for (int i=0; i<len; i++) y[i] = (int) x[i];
  return y;
}
  

int *ToInt(SEXP X) {
  bool ToFalse[1] = { false };
  return ToIntI(X, ToFalse, false);
}


double *ToRealI(SEXP X, bool *create) {
  KEY_type *KT = KEYT();
  if (TYPEOF(X) == REALSXP) { 
    *create = false;
    return REAL(X);
  }
  // TO DO !!
  //HELPINFO("Better use 'double' as storage mode (for one of the arguments).");
  int len = length(X); 
  double *y;
  if (create || KT->ToRealN < len) {
    y = (double *) MALLOC(sizeof(*y) * len);
    if (y == NULL) ERR1("not enough memory for an %d vector of doubles", len);
    if (!create) {
      FREE(KT->ToRealDummy);
      KT->ToRealDummy = y;
      KT->ToRealN = len;
    }
  } else y = KT->ToRealDummy;
  int *x;
  if (TYPEOF(X)==INTSXP) x=INTEGER(X); else x=LOGICAL(X);
  for (int i=0; i<len; i++) y[i] = (double) x[i];
  return y;
}

double *ToReal(SEXP X) {
  bool ToFalse[1] = { false };
 if (TYPEOF(X) == REALSXP) return REAL(X);
  return ToRealI(X, ToFalse);
}


SEXP testStrassen(SEXP A, SEXP B, SEXP truenrowA, SEXP truenrowB) {
  Long nrowA = nrows(A);
  Long nrowB = nrows(B);
  Long ncolA = ncols(A);
  Long ncolB = ncols(B);
  Long nrowR = nrowA < nrowB ? nrowB : nrowA;
  Long ncolR = ncolA > ncolB ? ncolA : ncolB;
  Long size = nrowR * ncolR;
  int *R = (int*) MALLOC(sizeof(int) * size);
  Minus_int(INTEGER(A),	INTEGER(truenrowA)[0], ncolA, nrowA,
	    INTEGER(B), INTEGER(truenrowB)[0], ncolB, nrowB,
	    nrowR, R);
  SEXP Ans;
  PROTECT(Ans=allocMatrix(REALSXP, nrowR, ncolR));
  double *ans = REAL(Ans);
  for (Long i=0; i<size; i++) ans[i] = (double) R[i];
  return Ans;
}

SEXP DivByRow(SEXP M, SEXP V) {
  Long
    l = length(V),
    r = nrows(M),
    c = ncols(M);

  double *m = REAL(M),
    *v = REAL(V);
  
  if (l != c) ERR0("vector does not match matrix");
  for (Long j=0; j<c; j++) {
    double vj = v[j];
    for (Long i=0; i<r; i++) {
      *(m++) /= vj;
    }
  }
  return M;
}


SEXP colMaxs(SEXP M) {
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  Long
    r = nrows(M),
    c = ncols(M);
  if (r == 0) return R_NilValue;
  SEXP Ans;
  if (TYPEOF(M) == REALSXP) {
    PROTECT(Ans = allocVector(REALSXP, c));
    colMaxsI(REAL(M), r, c, REAL(Ans), cores);
  } else {
    bool i = TYPEOF(M) == INTSXP;
    PROTECT(Ans = allocVector(i ? INTSXP : LGLSXP, c));
    int *m, *a;
    if (i) {
      m = INTEGER(M);
      a = INTEGER(Ans);
    } else {
      m = LOGICAL(M);
      a = LOGICAL(Ans);
    }
    colMaxsIint(m, r, c, a, cores);
  }
  UNPROTECT(1);
  return Ans;
}


SEXP rowProd(SEXP M) {
  Long
    r = nrows(M),
     c = ncols(M);
  if (r == 0) return R_NilValue;
  SEXP Ans;
  if (TYPEOF(M) == REALSXP) {
    PROTECT(Ans = allocVector(REALSXP, r));
    rowProdI(REAL(M), r, c, REAL(Ans));
   } else {
    // printf("type = %d", TYPEOF(M));
    RFERROR("transform to double first") ;
  }
  UNPROTECT(1);
  return Ans;
}


SEXP rowMeansX(SEXP M, SEXP Weight) {
  // todo : SSE2 / AVX
  Long
    L = length(Weight),
    r = nrows(M),
    c = ncols(M);
  double *weight = NULL;
  if (L != 0) weight = ToReal(Weight);
  if (r == 0 || c == 0) return R_NilValue;
  if (L != c && L != 0)
    ERR0("Length of 'weight' must equal number of columns of 'x'.");
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, r));
  DEF_VOID(M);
  rowMeansI(VOID(M), r, c, weight, REAL(Ans));
  if (L != 0 && TYPEOF(Weight) != REALSXP) { FREE(weight); }
  UNPROTECT(1);
  return Ans;
}


SEXP dbinorm(SEXP X, SEXP Sigma) { // 12'41
  Long nrow,
    ncol = 2;
  double *x, *y;
  if (TYPEOF(X) == VECSXP) {
    if (length(X) != ncol) BUG;
    SEXP xx = VECTOR_ELT(X, 0);
    nrow = length(xx);
    x = REAL(xx);
    y = REAL(VECTOR_ELT(X, 1));
  } else {
    if (isMatrix(X)) {
      if (ncols(X) != ncol) BUG;
      nrow = nrows(X);
    } else if (isVector(X)) {
      if (length(X) != ncol) BUG;
      nrow = 1;
    } else BUG;
    x = REAL(X);
    y = x + nrow;
  }
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, nrow));
  dbinormI(x, y, nrow, length(Sigma) == 0 ? NULL : REAL(Sigma), REAL(Ans));
  UNPROTECT(1);
  return Ans;
}




SEXP test(SEXP AA, SEXP CC, SEXP X) {
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  Long nrow = nrows(AA),
    ncol = ncols(AA),
    dim = length(X),
    k = MIN(ncol / 2, nrow),
    m = MAX(ncol, nrow);
  
  double
    eps = 1e-14,
    *A = REAL(AA),
    *C = REAL(CC),
    *x = REAL(X),    
    z[2],
    *a[2] = {(double*) MALLOC(m * m * sizeof(double)),
	       (double*) MALLOC(m * m * sizeof(double))};

  if (ncols(CC) != nrows(CC) ||  ncols(CC) != ncol) BUG;
  if (length(X) != nrow) BUG;

  for (int i=0; i<=17; i++) {
    for (int j=0; j<=1; j++) {
      SetLaMode(j == 0 ? LA_INTERN : LA_R, cores);
      switch(i) {
      case 1: z[j] = XkCXtl(A, C, nrow, ncol, nrow / 3, nrow / 4, cores); break;
      case 2: XCXt(A, C, a[j], nrow, ncol, cores); break;
      case 3: AtA(A, nrow, ncol, a[j], cores); break;
      case 4: xA(x, A, nrow, ncol, a[j], cores); break;
      case 5: xA_noomp(x, A, nrow, ncol, a[j]); break;
	//    case : xA(x1, x2,  A, nrow, ncol, a[j]1,  a[j]2); break;
      case 6: z[j] = xAx(x, C, nrow, cores); break;
      case 7: Ax(A, C, nrow, ncol, a[j], cores); break;// C genuegend lang. Reicht.
      //    case 8: Ax(A, x, x2, nrow, ncol, a[j]1,  a[j]2); break;
      case 8: z[j] =xUy(x, C, A, dim, cores); break; // A genuegend lang. Reicht.
      case 9: z[j] =xUxz(x, C, dim, a[j], cores); break;
      case 10: z[j] =x_UxPz(x, C, A, dim,cores); break; // A genuegend lang. Reicht.
      case 11: z[j] =xUx(x, C, dim, cores); break;
      case 12: matmult(A, C, a[j], nrow, ncol, k, cores); break;
      case 13: matmulttransposed(A, C, a[j], ncol, nrow, k, cores); break;
	//case : matmulttransposedInt(int *A, int *B, int *c, ncol, ncol, k); break; 
      case 14: matmult_2ndtransp(A, C, a[j], nrow, ncol, k, cores); break;
      case 15: matmult_2ndtransp(A, C, a[j], nrow, ncol, cores); break;
      case 16: matmult_tt(A, C, a[j], ncol, nrow, k,cores); break;
	//     case 17: z[j]=  scalar(A, C, ncol); break;	
      default: BUG;
      }

      int size = 0;
      switch(i) {
      case 1: case 6: case 8: case 9:case 10: case 11: case 17:
	if (FABS(z[0] - z[1])> eps) { PRINTF("i=%d", i); BUG; }
	break;
      case 2:  size = ncol * ncol;
	break;
      case 3: case 15: size = nrow * nrow;
	break;
      case 4 : case 5: size = ncol;
	break;
      case 7 : size = nrow;
	break;
      case 12: case 13: case 14: case 16: size = nrow * k;
	break;
      default: BUG;
      }
      for (int p=0; p<size; p++)
	if (FABS(a[0][p] - a[1][p]) > eps)  { PRINTF("i=%d, %d", i, p); BUG; }
    }
  }

 FREE(a[0]);
 FREE(a[1]);
  
  return R_NilValue;
}



SEXP quadratic(SEXP A, SEXP x) {
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  SEXP ans;
  int len = length(x);
  if (len != nrows(A) || len != ncols(A)) ERR0("'x' and 'A' do not match.");
  PROTECT(ans = allocVector(REALSXP, 1));
  REAL(ans)[0] = xAx(REAL(x), REAL(A), len, cores);
  UNPROTECT(1);
  return ans;
}

SEXP dotXV(SEXP M, SEXP V) {
  Long
    r = nrows(M),
    c = ncols(M),
    l = length(V)
    ;
  if (l != r) ERR0("X and v do not match");
  if (r == 0) return R_NilValue;
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, r, c));
  dotXVI(REAL(M), r, c,  REAL(V), REAL(Ans));
  UNPROTECT(1);
  return Ans;
}


SEXP debuggingLevel() {
  SEXP ans;
  PROTECT(ans = allocVector(INTSXP, 1));
#ifdef SCHLATHERS_MACHINE
  INTEGER(ans)[0] = 1;
#else
  INTEGER(ans)[0] = 0;
#endif  
  UNPROTECT(1);
  return ans;
}
// for debugging only
SEXP DebugCall() {
  //  return R_NilValue;
  //   KEY_type *KT = KEYT();						
  //  assert((KT->n_data_names == 0) xor (KT->data_names != NULL)); 
  //  assert((KT->n_coord_names == 0) xor (KT->coord_names != NULL));
  //  assert((KT->n_data_idx == 0) xor (KT->data_idx != NULL));	
  //  assert((KT->n_coord_idx == 0) xor (KT->coord_idx != NULL));
  return R_NilValue;
}




#define Nmodi 9
name_type modi = { "1x1", "2x2", "4x4", "8x8", "near", "simple", "precise", "kahan", "1x1p"};


SEXP scalarR(SEXP x, SEXP y, SEXP Mode) { // unused
  Long len = length(x);
  if (length(y) != len) ERR0("x and y differ in length");
  int mode;
  if (length(Mode) == 0) mode = -1;
  else if (INTSXP==TYPEOF(Mode)) mode = INTEGER(Mode)[0];
  else mode = Match((char*) CHAR(STRING_ELT(Mode, 0)), modi, Nmodi);
  SEXP Ans;

  if (isMatrix(x)) {
    Long nc = ncols(x);
    PROTECT(Ans = allocVector(REALSXP, nc * (nc - 1) / 2));
    double *ans = REAL(Ans);
    *ans = scalarX(REAL(x), REAL(y), len, 11); // no PROTECT(needed
    UNPROTECT(1);
  } else {
    PROTECT(Ans = allocVector(REALSXP, 1));
    double *ans = REAL(Ans);
    *ans = scalarX(REAL(x), REAL(y), len, mode); // no PROTECT(needed
    UNPROTECT(1);
  }
  return Ans;
}




SEXP crossprodX(SEXP X, SEXP Y, SEXP mode) {
  if (TYPEOF(X) != TYPEOF(Y)) ERR0("matrices do not have the same type");
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);

  if (cores == 0) BUG;
  
  Long n, nrow,
    len,
    lenY,
    ncol;
  if (isMatrix(X)) {
    nrow = ncols(X);
    len = nrows(X);
  } else {
    nrow = 1;
    len = length(X);
  }
  if (isMatrix(Y)) {
    ncol = ncols(Y);
    lenY = nrows(Y);
  } else {
    ncol = 1;
    lenY = length(Y);
  }
  if (lenY != len) ERR0("sizes of 'x' and 'y' do not match");
  if (length(mode) == 0) n = SCALAR_DEFAULT;
  else {
    n = INTEGER(mode)[0];
    if (n < 0) n =  SCALAR_DEFAULT;
  }
  SEXP Ans; 
  PROTECT(Ans = allocMatrix(REALSXP, nrow, ncol));
  double *ans = REAL(Ans);
 
  if (TYPEOF(X) == REALSXP) {
     double 
       *x = REAL(X),
       *y = REAL(Y);
     if (x == y) AtA(x, len, ncol, ans, cores,
		     KT->global_utils.solve.AtAmode);
     else  matmulttransposed(x, y, ans, len, nrow, ncol, cores);
  } else {
    int *x, *y;
    if (TYPEOF(X) == INTSXP) {
      x = INTEGER(X);
      y = INTEGER(Y);
    } else {
      x = LOGICAL(X);
      y = LOGICAL(Y);
    }
    if (nrow == ncol) {
      assert(sizeof(double) == sizeof(Long));
      AtAInt(x, y, len, ncol, len, ncol, (Long*) ans,
	     cores, &(KT->global_utils.solve));
      Long size = nrow  * ncol,
	sizeM4 = size - 4,
	i = 0;
      for ( ; i<sizeM4; i+=4) {
	ans[i] = (double) ((Long*) ans)[i];
 	ans[i+1] = (double) ((Long*) ans)[i+1];
 	ans[i+2] = (double) ((Long*) ans)[i+2];
 	ans[i+3] = (double) ((Long*) ans)[i+3];
     }
     for ( ; i<size; i++) ans[i] = (double) ((Long*) ans)[i];
    }
    else  ERR0("matrices do not have real-valued type");
  }

  UNPROTECT(1);
  return Ans;
}


#endif
