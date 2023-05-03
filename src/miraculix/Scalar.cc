
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


/*

Makefile must be:

PKG_LIBS =  $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)  -march=native  -mssse3 
PKG_CXXFLAGS =  $(SHLIB_OPENMP_CXXFLAGS)  -march=native -mssse3 

 */


#include "Basic_miraculix.h"
#include "xport_import.h"
#include "Scalar.h"
#include "MXinfo.h"
#include "options.h"
 

#define U(N) v1[N] * v2[N]

#define scalarNbyN(AtOnce, AndSoOn)				      \
  Long scalarInt##AtOnce##by##AtOnce( int * v1, int * v2, int N){ \
    int						\
    *endv1 = v1 + (N / AtOnce) * AtOnce,		\
      *end = v1 + N;							\
    Long sum = 0;							\
    for (; v1 < endv1; v1 += AtOnce, v2 += AtOnce)			\
      sum += U(0) + U(1) + U(2) + U(3) AndSoOn;				\
    for (; v1 < end; v1++, v2++) sum +=  v2[0] * v1[0];			\
    return sum;								\
}

scalarNbyN(4,)
scalarNbyN(8, + U(4) + U(5) + U(6) + U(7))
scalarNbyN(16, + U(4) + U(5) + U(6) + U(7) + U(8) + U(9) + U(10) + U(11) +
	   U(12) + U(13) + U(14) + U(15))



//bool pr = true;
Long scalarInt(int *x, int *y, int len, int n) {
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
  case SCALAR_INT_8 : return scalarInt8by8(x, y, len); 
  case SCALAR_INT_16 : return scalarInt16by16(x, y, len);
  case SCALAR_INT_AVX2 : return scalarIntAVX2(x, y, len);
  default : ERR0("coding not available"); 
  }
  return NA_LONG;
}

void matmulttransposedInt(int *A, int *B, double *c, int m, int l, int n,
			   int VARIABLE_IS_NOT_USED cores) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),  
// saving result in C
  if (A == B && n == l) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20)
#endif
    for (Long k=0; k<l; k++) {    
      double *C = c + k;
      int *Aim = A + k * m;
      for (Long j=k; j<n; j++)
	c[j + k * l] = C[j * l] = (double) SCALARINT(Aim, B + j * m, m);
    }
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
    for (Long k=0; k<l; k++) {    
      double *C = c + k;
      int *Aim = A + k * m;
      for (Long j=0; j<n; j++)
	C[j * l] = (double) SCALARINT(Aim, B + j * m,m);
    }
  }
}

void crossprod_Int(int *x, int *y, int nrow, int ncol, int len,
		   int VARIABLE_IS_NOT_USED cores, int *ans) {
  if (x == y) {
    assert(nrow == ncol);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20)
#endif
    for (Long k=0; k<nrow; k++) {    
      int *C = ans + k,
	*Aim = x + k * len;
      for (Long j=k; j<nrow; j++) {
	Long s = SCALARINT(Aim, x + j * len, len);
	if (s > MAXINT || s < MAXINT) BUG;
	ans[j + k * nrow] =  C[j * nrow] = (int) s;
      }
	
    }
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
    for (Long k=0; k<nrow; k++) {    
      int *C = ans + k,
	*Aim = x + k * len;
      for (Long j=0; j<ncol; j++) {
	Long s = SCALARINT(Aim, y + j * len, len);
	if (s > MAXINT || s < MAXINT) BUG;
	C[j * nrow] = (int) s;
      }
    }
  }
}

