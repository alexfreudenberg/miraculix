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

#include "miraculix.h"
#include "xport_import.h"
#include "MXinfo.h"
#include "options.h"
#include "kleinkram.h"
//#include "haplogeno.h"
#include "Files.h"
#include "Template.h"
//#include "2bit.h"
#include "Vector.matrix.h"

SEXP matrixvector012(SEXP matrix, SEXP vector) {
  int
    n = LENGTH(vector),
    r = nrows(matrix),
    c = ncols(matrix)
    ;
  if (c != n) ERR0("vector and matrix do not match");
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, r));
  matrixvector012I( REAL(matrix), r, c, vector, REAL(ans));
   
  UNPROTECT(1);
  return ans;
}


SEXP vector012matrix(SEXP vector, SEXP matrix) {
  int 
    n = LENGTH(vector),
    r = nrows(matrix),
    c = ncols(matrix)
    ;
  if (r != n) ERR0("vector and matrix do not match");
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, c));
  vector012matrixI(vector, REAL(matrix), r, c,  REAL(ans));
  UNPROTECT(1);
  return ans;
} 



SEXP vectorGeno(SEXP V, SEXP SxI) {
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  tuning_options *tuning = &(global->tuning);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
 

  if (TYPEOF(SxI) == STRSXP)
    return file_intern(SxI, DotFile, 0, global, utils, false, V);

  Long *info = GetInfo(SxI, true);
  if (info == NULL) return R_NilValue;
  Long
    repetV = isMatrix(V) ? ncols(V) : 1,
    snps = info[SNPS],
    individuals = info[INDIVIDUALS];
  if ( LENGTH(V) % snps !=0) ERR0("vector 'V' not of correct length");  

  SEXP Ans;
  if (repetV == 1) Ans = PROTECT(allocVector(REALSXP, individuals));
  else Ans = PROTECT(allocMatrix(REALSXP, individuals, repetV));

  vectorGeno_means_double(SxI, REAL(V), repetV, global, utils, REAL(Ans));
     
  UNPROTECT(1);
  return Ans;
}
 



// utils_options global_utils;
// Ext_get_utilsoption(&global_utils, false);

SEXP genoVector(SEXP SxI, SEXP V) {
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  tuning_options *tuning = &(global->tuning);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);

  if (TYPEOF(SxI) == STRSXP) 
    return file_intern(SxI, FileDot, 0, global, utils, false, V);
  
  Long *info = GetInfo(SxI, true);
  if (info == NULL) return R_NilValue;
  
  Long
    repetV = isMatrix(V) ? ncols(V) : 1,
    snps = info[SNPS],
    individuals = info[INDIVIDUALS];
  
  if ( LENGTH(V) % individuals != 0) ERR0("Vector 'V' not of correct length");
  
  SEXP Ans;
  if (repetV == 1) Ans = PROTECT(allocVector(REALSXP, snps));
  else Ans = PROTECT(allocMatrix(REALSXP, snps, repetV));

  
  genoVector_means_double(SxI, REAL(V), repetV, global, utils, REAL(Ans));

  UNPROTECT(1);
  return Ans;
}



SEXP IsolveRelMat(Long individuals, double *Aorig, double *tau,
		  int ntau,
		  double *Vec,
		  double *beta, int nbeta,
		  int returns, bool destroy,
		  int VARIABLE_IS_NOT_USED cores) {
  // returns: number of return values; could be 3,2,1
  // if more than 1 element is returned, a list is returned
   const char *names[3] = {"rest", "yhat", "rel.matrix"};
  double
     *X = NULL,
    *pA = NULL,
    *pAtau = NULL;

  Long
    protects = 0,
    iP1 = individuals + 1,
    i2 = individuals * individuals;
  SEXP yhat=R_NilValue,
    Ans=R_NilValue,
    rest, RA = R_NilValue;
  PROTECT(rest = allocVector(REALSXP, individuals));
  protects++;
  double *r = REAL(rest);
  MEMCOPY(r, Vec, sizeof(*r) * individuals);

  if (returns == 1) { // only "rest" is returned
    if (destroy) pAtau = Aorig;
    else {
      X =(double *) MALLOC(sizeof(*X) * i2);
      MEMCOPY(X, Aorig, sizeof(*X) * i2);
      pAtau = X;     
    }
    
  } else { // at least "rest" and "yhat" are returned
    SEXP namevec;
    PROTECT(Ans = allocVector(VECSXP, returns));
    protects++;
    PROTECT(namevec = allocVector(STRSXP, returns));
    protects++;
    for (Long k=0; k<returns; k++) SET_STRING_ELT(namevec, k, mkChar(names[k]));
    setAttrib(Ans, R_NamesSymbol, namevec);
    SET_VECTOR_ELT(Ans, 0, rest);
    PROTECT(yhat = allocVector(REALSXP, individuals));
    protects++;
    SET_VECTOR_ELT(Ans, 1, yhat);
  
    if (returns == 2) {
      X = (double *) MALLOC(sizeof(*X) * i2);
      MEMCOPY(X, Aorig, sizeof(*X) * i2);
      pAtau = X;
      pA = Aorig;
    }

    else if (returns == 3) {
      PROTECT(RA = allocMatrix(REALSXP, individuals, individuals));
      protects++;
      SET_VECTOR_ELT(Ans, 2, RA);
      MEMCOPY(REAL(RA), Aorig, sizeof(*Aorig) * i2);
      if (destroy) {
	pAtau = Aorig;
	pA = REAL(RA);
      } else {
	pAtau = REAL(RA);
	pA = Aorig;
      }
    }

    else BUG;
  }

  if (ntau > 0) {
    if (tau[0] <= 0) ERR0("'tau' must be positive");
    if (ntau != 1 && ntau != individuals)
      WARN2("length of tau (%d) differs from the number of individuals (%ld).",
	    ntau, individuals);
    for (Long k=0; k<i2; k+=iP1) pAtau[k] += tau[k % ntau];
  }
  
  solve_storage *Ssolve = (solve_storage*) MALLOC(sizeof(solve_storage));
  Ext_solve_NULL(Ssolve); 
  solve_options Soption;						       
  MEMCOPY(&Soption, &(KEYT_M()->global_utils.solve), sizeof(solve_options));
  Soption.Methods[0] = Cholesky;					      
  Soption.Methods[1] = NoFurtherInversionMethod;
  Rint err = Ext_SolvePosDefSp(pAtau, individuals, true,  r, 1, NULL,
			       Ssolve, &Soption, cores);
  FREE(X);
  if (err != NOERROR) {
    errorstring_type errorstring;
    STRNCPY(errorstring, Ssolve->err_msg, MAXERRORSTRING);
    Ext_solve_DELETE(&Ssolve);
    ERR1("error occurred when solving the system: %.50s)", errorstring);
  }
  Ext_solve_DELETE(&Ssolve);
  
  if (returns == 1) {
    UNPROTECT(protects);
    return rest;
  }
  
  double *y = REAL(yhat);
  xA(r, pA, individuals, individuals, y, cores);
  
 if (nbeta > 0) {
    if (nbeta != 1 && nbeta != individuals)
      WARN2("length of beta (%d) differs from the number of individuals (%ld).",
	    nbeta, individuals);
    for (Long i=0; i<individuals; i++) y[i] += beta[i % nbeta];  // ok
 }

  if (returns == 3 && !destroy) MEMCOPY(REAL(RA), Aorig, sizeof(*Aorig) * i2);
  
  UNPROTECT(protects);
  return Ans;
}


SEXP solveRelMat(SEXP A, SEXP Tau, SEXP Vec, SEXP Beta, SEXP Destroy) {
  KEY_type *KT = KEYT_M();
  int cores = GreaterZero(KT->global_utils.basic.cores);
  Long individuals = nrows(A);
  bool lenBeta = LENGTH(Beta);
 return IsolveRelMat(individuals, REAL(A),//no PROTECT( needd
		      REAL(Tau), LENGTH(Tau), REAL(Vec),
		      REAL(Beta), lenBeta,
		     1 + (int) (lenBeta > 0), LOGICAL(Destroy)[0], cores);
}

#endif
