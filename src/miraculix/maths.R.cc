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

#include <R_ext/Lapack.h>
#include "RandomFieldsUtils.h"
#include "zzz_RFU.h"
#include "utils.h"
#include "xport_import_RFU.h"



double struve_intern(double x, double nu, double factor_Sign, bool expscaled);
SEXP struve(SEXP X, SEXP Nu, SEXP Factor_Sign, SEXP Expscaled) {
  int i,    
    lenx = length(X),
    lennu = length(Nu),
    len = lenx;  
  if (len < lennu) len = lennu;
  SEXP Result;
  PROTECT(Result = allocVector(REALSXP, len));
  double *x = REAL(X),
    *nu  = REAL(Nu),
    factor_sign = REAL(Factor_Sign)[0],
    *result = REAL(Result);
  bool expscaled = LOGICAL(Expscaled)[0];
  for (i=0; i<len; i++)
    result[i]=struve_intern(x[i % lenx], nu[i % lennu], factor_sign, expscaled);
 
  UNPROTECT(1);
  return Result;
}



SEXP I0ML0(SEXP X) {
  SEXP Result;
  PROTECT(Result = allocVector(REALSXP, length(X)));
  double *x = REAL(X),
    *result = REAL(Result);
  int i,    
    lenx = length(X);  
  for (i=0; i<lenx; i++) result[i] = I0mL0(x[i]);

  UNPROTECT(1);
  return Result;
}



typedef double (*primfct1)(double);
typedef double (*primfct3)(double, double, double, whittle_work_type*);
#define CALCULATE(PRIMFCTN)			\
  double *x = REAL(X);				\
  int n = length(X),				\
    deriv = INTEGER(Derivative)[0];					\
  if (deriv < 0 || deriv > 4) ERR0("value of 'derivative' out of range"); \
  PRIMFCTN F = fctns[deriv];						\
									\
  SEXP Ans;								\
  PROTECT(Ans=allocVector(REALSXP, n));					\
  double *ans = REAL(Ans);						\
  for (int i=0; i<n; i++) ans[i] = F

#define RETURN					\
  UNPROTECT(1);					\
  return(Ans);


SEXP gaussr(SEXP X, SEXP Derivative) {  
  static primfct1 fctns[] = {Gauss, DGauss, DDGauss, D3Gauss, D4Gauss};
  CALCULATE(primfct1)(FABS(x[i]));
  RETURN;
}

SEXP WMr(SEXP X, SEXP Nu, SEXP Derivative, SEXP Factor) {  
  static primfct3 fctns[] = {WM, DWM, DDWM, D3WM, D4WM };
  double 
    *nu = REAL(Nu),
    *factor = REAL(Factor);
  int 
    nnu = length(Nu),
    nfactor = length(Factor);  
  CALCULATE(primfct3)(FABS(x[i]), nu[i % nnu], factor[i % nfactor], NULL);
  RETURN;
}
 

SEXP logWM2r(SEXP X, SEXP Nu1, SEXP Nu2, SEXP Factor) {  
  double 
    nu1 = REAL(Nu1)[0],
    nu2 = REAL(Nu2)[0],
    factor = REAL(Factor)[0];
  double *x = REAL(X);				
  //  int n = length(X);	
  if (nu1 <= 0.0 || nu2 <= 0.0) ERR0("'nu' must be positive");
  if (factor < 0.0) ERR0("'factor' must be positive");
 									
  SEXP Ans;								
  PROTECT(Ans=allocVector(REALSXP, 1));					
  double *ans = REAL(Ans);						
  ans[0] = logWM2(FABS(x[0]), nu1, nu2, factor, NULL);
  UNPROTECT(1);					
  return(Ans);
}



SEXP besselk_simd(SEXP X, SEXP Nu) {
  if (length(Nu) != 1) ERR0("Length of nu must be 1.");
  SEXP Ans = R_NilValue;
#if defined AVX2
  int L = length(X);
  PROTECT(Ans = allocVector(REALSXP, L));
  bes_k_simd(REAL(X), REAL(Nu)[0], L, REAL(Ans));
  UNPROTECT(1);
#else
  ERR0("'AVX2' not available.");
#endif  
  return Ans;
}


#endif // compatibility_to_R

