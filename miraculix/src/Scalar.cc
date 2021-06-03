
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2019 -- 2019 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*

Makefile must be:

PKG_LIBS =  $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)  -march=native  -mssse3 
PKG_CXXFLAGS =  $(SHLIB_OPENMP_CXXFLAGS)  -march=native -mssse3 

 */


#include "Scalar.h"
#include <General_utils.h>
#include "intrinsics.h"
#include "miraculix.h"


#define U(N) v1[N] * v2[N]

#define scalarNbyN(AtOnce, AndSoOn)				      \
  Uint scalarUint##AtOnce##by##AtOnce( Uint * v1, Uint * v2, Uint N){ \
    Uint						\
    *endv1 = v1 + (N / AtOnce) * AtOnce,		\
      *end = v1 + N,					\
      sum = 0.0;							\
    for(; v1 < endv1; v1 += AtOnce, v2 += AtOnce)			\
      sum += U(0) + U(1) + U(2) + U(3) AndSoOn;				\
    for(; v1 < end; v1++, v2++) sum +=  v2[0] * v1[0];			\
    return sum;								\
}

scalarNbyN(4,)
scalarNbyN(8, + U(4) + U(5) + U(6) + U(7))
scalarNbyN(16, + U(4) + U(5) + U(6) + U(7) + U(8) + U(9) + U(10) + U(11) + U(12) + U(13) + U(14) + U(15))

#if defined AVX2
Uint scalarUintAVX2(Uint * V1, Uint * V2, Uint N){
  Uint steps = N / UnitsPerBlock;
  BlockType0    
    *v1 = (BlockType0 *) V1,
    *endv1 = v1 + steps,
    *v2 = (BlockType0 *) V2;
  BlockUnitType sum;
  ZERO(sum.vi);

  Uint zaehler = 0;

  for(; v1 < endv1; v1 ++, v2 ++) {
    zaehler++;
    BlockType dummy, s1, s2;
    LOADU(s1, v1);
    LOADU(s2, v2);    
    MULT32(dummy, s1, s2);
    ADD32(sum.vi, sum.vi, dummy);
  }
  
  Uint
    totalsum = (sum.u32[0] + sum.u32[1] + sum.u32[2] + sum.u32[3] +
		sum.u32[4] + sum.u32[5]	+ sum.u32[6] + sum.u32[7]
		),
    total = steps * UnitsPerBlock,
    *end = V1 + N;  				
  V1 += total;
  V2 += total;
  for(; V1 < end; V1++, V2++) totalsum += V2[0] * V1[0];			
  return totalsum;								
}
#else
Uint scalarUintAVX2(Uint * V1, Uint * V2, Uint N) { BUG; }
#endif


//bool pr = true;
Uint scalarUint(Uint *x, Uint *y, Uint len, Uint n) {
 //  if (pr) { prUintf("mode = %d\n", n); pr = false; }
 // 0 : 7.9
// 1:  7.55
// 2: 7.8
// 3:7.58
//4: 7.5 
// 5: 7.4!
//6:7.4
//7: 7.9
// 8: "ewige" schleife
  switch(n) {
  case SCALAR_INT_8 : return scalarUint8by8(x, y, len); 
  case SCALAR_INT_16 : return scalarUint16by16(x, y, len);
  case SCALAR_INT_AVX2 : return scalarUintAVX2(x, y, len);
  default : ERR("method not available"); 
  }
  return RF_NAN;
}
  

 

void matmulttransposedUint(Uint *A, Uint *B, double *c, Uint m, Uint l, Uint n) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),  
// saving result in C
  if (A == B && n == l) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 20)
#endif
    for (Uint i=0; i<l; i++) {    
      double *C = c + i;
      Uint *Aim = A + i * m;
      for (Uint j=i; j<n; j++)
	c[j + i * l] = C[j * l] = (double) SCALARUINT(Aim, B + j * m, m);
    }
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) 
#endif
    for (Uint i=0; i<l; i++) {    
      double *C = c + i;
      Uint *Aim = A + i * m;
      for (Uint j=0; j<n; j++) C[j * l] = (double) SCALARUINT(Aim, B + j * m,m);
    }
  }
}



SEXP crossprodInt(SEXP X, SEXP Y, SEXP mode) {
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
    len = length(X);
  }
  if (isMatrix(Y)) {
    ncol = ncols(Y);
    lenY = nrows(Y);
  } else {
    ncol = 1;
    lenY = length(Y);
  }
  if (lenY != len) ERR("sizes of 'x' and 'y' do not match");
  if (length(mode) == 0) n = SCALAR_INT_DEFAULT;
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

  
  if (x == y) {
    assert(nrow == ncol);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 20)
#endif
    for (Uint i=0; i<nrow; i++) {    
      Uint *C = ans + i,
	*Aim = x + i * len;
      for (Uint j=i; j<nrow; j++)
	ans[j + i * nrow] = C[j * nrow] = SCALARUINT(Aim, x + j * len, len);
    }
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) 
#endif
    for (Uint i=0; i<nrow; i++) {    
      Uint *C = ans + i,
	*Aim = x + i * len;
      for (Uint j=0; j<ncol; j++)
	C[j * nrow] = SCALARUINT(Aim, y + j * len, len);
    }
  }

  
  UNPROTECT(1);
  return Ans;
}



