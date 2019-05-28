
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2014 -- 2018  Martin Schlather
Copyright (C) 2014 -- 2015 Florian Skene: SSE2+SSE3

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#include "miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "AutoMiraculix.h"
#include <General_utils.h>
#include "intrinsics.h"


SEXP CreateEmptyCodeVector(Uint snps, Uint individuals, Uint mem,
			   SEXPTYPE type, SEXP Information, bool createP) {
  SEXP Code, info, Infos, P, PMT;
  Uint
    sizeP = snps + INFO_P_P;
  
  PROTECT(Code = allocVector(type, mem));
  if (type == REALSXP) {
    double * A = REAL(Code);
    for (Uint k=0; k<mem; A[k++] = 0.0);
  } else if (type == INTSXP) {
    int* A = INTEGER(Code);
    for (Uint k=0; k<mem; A[k++] = 0.0);
  } else BUG;
  
  PROTECT(Infos = allocVector(VECSXP, INFO_LAST + 1));
  for(Uint i=0; i<=INFO_LAST; SET_VECTOR_ELT(Infos, i++, R_NilValue));

  if (createP) {
    PROTECT(PMT = allocVector(REALSXP, individuals));
    double *pd = REAL(PMT);
    for (Uint i=0; i < individuals; pd[i++] = RF_NA);
    
    PROTECT(P = allocVector(REALSXP, sizeP));
    pd = REAL(P);
    for (Uint i=0; i < sizeP; pd[i++] = 0.0);

    SET_VECTOR_ELT(Infos, INFO_P, P);
  } 
 
  PROTECT(info = allocVector(INTSXP, INFO_INFO_LAST + 1));
  int *pI = INTEGER(info);
  for (Uint i=0; i<=INFO_INFO_LAST; pI[i++] = NA_INTEGER);
  pI[WHAT] = GENOMATRIX;
  pI[SNPS] = snps;
  pI[INDIVIDUALS] = individuals;  
  pI[SNPxIND] = true;
  pI[MEM] = mem;
  
  SET_VECTOR_ELT(Infos, INFO_INFO, info);
  setAttrib(Code, Information, Infos); 
  UNPROTECT(3 + createP * 2);
  return Code;
}

double DoCentering(double *ans, Uint individuals, bool centred, bool normalized,
		 double sumGeno) {  
  double
    nd = (double) individuals,
    nSq = nd * nd,
    factor = RF_NAN,
    sum_pi_qi = RF_NAN;
  if (centred || normalized) {
    double
      sum_P = RF_NAN,
      *sums = (double *) MALLOC(sizeof(double) * individuals),
      total = 0.0;
    
    Uint k = 0,
      nP1 = individuals + 1;
    for (Uint i=0; i<individuals; i++) {
      double dummy = 0.0;
      for (Uint j=0; j<individuals; j++) dummy += ans[k++];
      sums[i] = nd * dummy;
      total += dummy;
    }

    if (centred) {
      for (Uint i=0; i<individuals; i++) {
	Uint idx = i * nP1;
	for (Uint j=i; j<individuals; j++, idx++) {
	  ans[idx] = nSq * ans[idx] - sums[i] - sums[j] + total;
	}
      }
    }
    
    if (normalized) {
      sum_P = nd * sumGeno;
      factor = sum_pi_qi = sum_P - 0.5 * total;
      if (!centred) factor /= nSq;
    } else factor = nSq;
    
    //    print("sumpq=%10g  %10g; sumGeno=%10g\n", sum_pi_qi, sum_pi_qi * nSq, sumGeno);
    
    //    print("sigma2 = %10g\n", sum_pi_qi);
    
    for (Uint i=0; i<individuals; i++) {
      Uint idx = i * nP1;
      for (Uint j=i; j<individuals; j++) {
	ans[i + j * individuals] = (ans[idx++] /= factor);
      }
    }
    
    FREE(sums);
  }
  return sum_pi_qi /nSq ;
}

  
SEXP Sigma2 = R_NilValue;
void  DoCentering(SEXP Ans, Uint individuals, double sumGeno) {
  double
    sum_pi_qi = DoCentering(REAL(Ans), individuals,
			    GLOBAL.relationship.centered == True,
			    GLOBAL.relationship.normalized == True, sumGeno);

  if (GLOBAL.relationship.returnsigma) {
    if (Sigma2 == R_NilValue) Sigma2 = install("sigma2");
    setAttrib(Ans, Sigma2,  ScalarReal(sum_pi_qi));
  }
  
}
  
  
