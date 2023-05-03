/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2015 - 2017 Martin Schlather, Reinhard Furrer, Martin Kroll
 Copyright (C) 2017 - 2020 Martin Schlather
 Copyright (C) 2021 - 2022 Martin Schlather, Alexander Freudenberg
  
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

#include "RandomFieldsUtils.h"
#include "zzz_RFU.h"
#include "kleinkram.h" 
#include "xport_import_RFU.h"
#include "extern_RFU.h"
#include "solve_intern.h"


SEXP doPosDef(SEXP M, SEXP rhs, SEXP logdet, int calculate,
	      solve_storage *Pt, solve_options *Sp, int VARIABLE_IS_NOT_USED cores){
  // rhs_cols == 0 iff RHS = NULL
  int rhs_rows, rhs_cols,
    size = ncols(M);
  if (nrows(M) != size) ERR0("not a square matrix");
  int
    err = NOERROR;
   bool deleteMM = false,
    deleteRHS = false;
  SEXP res;
  solve_storage Pt0,
    *pt = Pt;
  if (pt == NULL) {
    solve_NULL(&Pt0);
    pt = &Pt0;
  }

  if (rhs == R_NilValue) {
    rhs_rows = rhs_cols = 0;
  } else if (isMatrix(rhs)) {
    rhs_rows = nrows(rhs);
    rhs_cols = ncols(rhs);
  } else if ((rhs_rows = length(rhs)) == 0) {
    rhs_cols = 0;
  } else {
    rhs_cols = 1;
  }
  if (rhs_rows > 0 && rhs_rows != size)
    ERR0("vector size does not match the matrix size");
  
  int new_cols = rhs_cols == 0 ? size : rhs_cols;
  Long total = (Long) size * new_cols;

  //  res =  PROTECT(isReal(M) ? duplicate(M): coerceVector(M, REALSXP)); UNPROTECT(1); return res;

  if (rhs_cols==0 || isMatrix(rhs)) {
    res = PROTECT(allocMatrix(REALSXP, size, new_cols));
  } else {
    res =  PROTECT(allocVector(REALSXP, total));
  }


  double *MM=NULL, 
    *RHS = NULL;
  if (TYPEOF(M) != REALSXP) {
    if (TYPEOF(M) != INTSXP && TYPEOF(M) != LGLSXP) 
      GERR0("numerical matrix expected");
    if ((deleteMM = rhs_cols != 0))
      MM = (double*) MALLOC(total * sizeof(*MM));
    else MM = REAL(res);
    if (TYPEOF(M) == INTSXP) {
      for (Long i=0; i<total; i++) 
	MM[i] = INTEGER(M)[i] == NA_INTEGER ? RF_NA : (double) INTEGER(M)[i];
    } else {
      for (Long i=0; i<total; i++) 
	MM[i] = LOGICAL(M)[i] == NA_LOGICAL ? RF_NA : (double) LOGICAL(M)[i];
    } 
  } else MM = REAL(M); 

  if (rhs_cols > 0) {
    if ((deleteRHS = TYPEOF(rhs) != REALSXP)) {
      if (TYPEOF(rhs) != INTSXP && TYPEOF(rhs) != LGLSXP) 
	GERR0("numerical matrix expected");
      Long totalRHS = (Long) rhs_cols * rhs_rows; 
      RHS = (double*) MALLOC(totalRHS * sizeof(*RHS));
      if (TYPEOF(rhs) == INTSXP) {
	for (Long i=0; i<totalRHS; i++) 
	  RHS[i] = INTEGER(rhs)[i] == NA_INTEGER 
	    ? RF_NA : (double) INTEGER(rhs)[i];
      } else if (TYPEOF(rhs) == LGLSXP) {
	for (Long i=0; i<totalRHS; i++) 
	  RHS[i] = LOGICAL(rhs)[i] == NA_LOGICAL
	    ? RF_NA : (double) LOGICAL(rhs)[i];
      } 
    } else RHS = REAL(rhs);
  }


  //printf("length = %d\n", length(logdet) == 0 );
  err = doPosDefIntern(MM, size, true, // no PROTECT( needed
		       rhs_cols == 0 ? NULL : RHS,//rhs_cols == 0 iff RHS = NULL
		       rhs_cols, 
		       (rhs_cols == 0 && TYPEOF(M) == REALSXP) ||
		       (rhs_cols > 0 && TYPEOF(rhs) == REALSXP) ? REAL(res)
		       : NULL, 
		       length(logdet) == 0 ? NULL : REAL(logdet),
		       calculate, pt, Sp, cores);

 ErrorHandling:
  if (deleteMM) { FREE(MM); }
  if (deleteRHS) { FREE(RHS); }
  if (pt != Pt) solve_DELETE0(pt);
  
  UNPROTECT(1);
  if (err != NOERROR) {
    const char *methname[] = {"solvePosDef", "cholesky", "determinant"};
    errorstring_type msg;
    switch (err) {
    case ERRORMEMORYALLOCATION : STRCPY(msg, "memory allocation error"); break;
    case ERRORNOTPROGRAMMEDYET : STRCPY(msg, "not programmed yet"); break;
    case ERRORFAILED : STRCPY(msg, "algorithm has failed"); break;
    case ERRORM : STRCPY(msg, pt->err_msg);
      break;
    default:  STRCPY(msg, "<unknown error>");
    }
    RFERROR2("'%.200s': %.200s.\n", methname[calculate], msg);    
  }

  return res;
}


SEXP SolvePosDefR(SEXP M, SEXP rhs, SEXP logdet){
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  // rhs_cols == 0 iff RHS = NULL
  return doPosDef(M, rhs, logdet, SOLVE, NULL,
		  &(KT->global_utils.solve), cores);
}



SEXP Chol(SEXP M) {
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  solve_options sp;
  MEMCOPY(&sp, &(OPTIONS.solve), sizeof(solve_options));
  sp.Methods[0] = sp.Methods[1] = Cholesky;
  sp.sparse = False; // currently does not work, waiting for Reinhard
  solve_storage Pt;
  solve_NULL(&Pt);
  SEXP Ans;
  PROTECT(Ans = doPosDef(M, R_NilValue, R_NilValue, MATRIXSQRT, &Pt, &sp,
			 cores));

  if (Pt.actual_pivot == PIVOT_DO || Pt.actual_pivot ==  PIVOT_IDX) {    
    // NEVER: FREE(OPTIONS.solve.pivot_idx); See Pivot_Cholesky:
    SEXP Idx, Info1, Info3;
    PROTECT(Idx = allocVector(INTSXP, Pt.n_pivot_idx));
    MEMCOPY(INTEGER(Idx), Pt.pivot_idx,
	    sizeof(*(Pt.pivot_idx)) * Pt.n_pivot_idx);
    setAttrib(Ans, install("pivot_idx"), Idx);
    
    PROTECT(Info1 = allocVector(INTSXP, 1));
    INTEGER(Info1)[0] = Pt.actual_size;
    setAttrib(Ans, install("pivot_actual_size"), Info1);
  
    PROTECT(Info3 = allocVector(INTSXP, 1));
    INTEGER(Info3)[0] = PIVOT_DO;
    setAttrib(Ans, install("actual_pivot"), Info3);
   
    UNPROTECT(3);
    assert(Pt.n_pivot_idx == ncols(M));
  }
  
  solve_DELETE0(&Pt);
  UNPROTECT(1);
  return Ans;
}



SEXP tcholRHS(SEXP C, SEXP RHS) {  
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  int n_protect = 2;
  SEXP Ans, Idx;
  PROTECT(Idx = getAttrib(C, install("pivot_idx")));
  bool pivot = length(Idx) > 0;
  int
    n = isMatrix(RHS) ? ncols(RHS) : 1,
    rows = isMatrix(RHS) ? nrows(RHS) : length(RHS),
    size = ncols(C),
    act_size =size;
  if (pivot) {
    SEXP dummy;
    PROTECT(dummy = getAttrib(C, install("pivot_actual_size")));
    act_size=INTEGER(dummy)[0];
    n_protect++;
  }
  int *pi = pivot ?  (int *) INTEGER(Idx) : NULL;
    if (isMatrix(RHS)) PROTECT(Ans = allocMatrix(REALSXP, size, n));
    else PROTECT(Ans = allocVector(REALSXP, size));
  if (rows < act_size) ERR0("too few rows of RHS");
  sqrtRHS_Chol(REAL(C), size, REAL(RHS), rows, n, REAL(Ans),
	       pivot, act_size, cores, pi);
  UNPROTECT(n_protect);
  return Ans;
}



SEXP chol2mv(SEXP C, SEXP N) {  
  KEY_type *KT = KEYT();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  int n_protect = 2;
  SEXP Ans, Idx;
  PROTECT(Idx= getAttrib(C, install("pivot_idx")));
  bool pivot = length(Idx) > 0;
  int
    n = INTEGER(N)[0],
    size = ncols(C),
    act_size = size;
  if (pivot) {
    SEXP dummy;
    PROTECT(dummy = getAttrib(C, install("pivot_actual_size")));
    act_size = INTEGER(dummy)[0];
    n_protect++;
  }
  Long n_act_size = (Long) n * act_size;
  int *pi = pivot ? INTEGER(Idx) : NULL;
  if (n == 1) PROTECT(Ans = allocVector(REALSXP, size));
  else PROTECT(Ans = allocMatrix(REALSXP, size, n));
  double *gauss = (double *) MALLOC(sizeof(double) * n_act_size);
  if (gauss == NULL) ERR0("memory allocation error");
  GetRNGstate();
  for (Long i=0; i<n_act_size; gauss[i++] = GAUSS_RANDOM(1.0));
  PutRNGstate();
  sqrtRHS_Chol(REAL(C), size, gauss, act_size, n, REAL(Ans),
	       pivot, act_size, cores, pi);
  FREE(gauss);
  UNPROTECT(n_protect);
  return Ans;
}



#endif
