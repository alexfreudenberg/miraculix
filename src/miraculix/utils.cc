
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
#include "compatibility.SEXP.h"
#include "compatibility.lapack.h"
#include "kleinkram.h"
#include "zzz_RFU.h"
#include "utils.h"
//#include "xport_import_RFU.h"
#include "extern_RFU.h"
  

AVAILABLE_SIMD

void AtAIntBlock2x2(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		    int VARIABLE_IS_NOT_USED cores, int scalarVersion,
		    Long m, Long n);

void AtAIntBlock3x3(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		    int VARIABLE_IS_NOT_USED cores, int scalarVersion,
		    Long m, Long n);

void Strassen(int *A, int  *B,
	      Long nrow, Long ncol, // a and b have same size 
	      Long ld,
	      Long ldC,
	      Long *C,
	      int cores, solve_options *options);

void AtAInt(int *a, int *b,
	    Long nrow, Long ncol, // a and b have same size
	    Long ld, Long ldC,
	    Long *C, // result
	    int VARIABLE_IS_NOT_USED cores,
	    solve_options *options) {
  //  printf("AtaInt:Mode=%d %ld %ld; %ld %ld %ld %ld; C=%ld\n", options->AtAmode, (Long) a, (Long) b, nrow, ncol, ld, ldC, (Long) C);
  ASSERT_SOLVE(options);
  Long 
    m = options->AtAnrow,
    n = options->AtAncol;
  switch(options->AtAmode) {
  case 10:// tiles, aber ungekreuzt
    AtAInt(a, b, nrow, ncol, ld, ldC, C, cores, m, n); 
    break;
  case 12: // tiles und gekreuzt 2x2
    if (!avx2Avail) ERR("method for matrix multiplication needs AVX2");
    AtAIntBlock2x2(a, b, nrow, ncol, ld, ldC, C, cores, SCALAR_VERSION, m, n);
    break;
  case 13: // tiles und gekreuzt 3x3
    if (!avx2Avail) ERR("method for matrix multiplication needs AVX2");
    AtAIntBlock3x3(a, b, nrow, ncol, ld, ldC, C, cores, SCALAR_VERSION, m, n);
    break;
  case 20:// strassen
    Strassen(a, b, nrow, ncol, ld, ldC, C, cores, options);
    return;
  case 33: // trial -- hat gar nichts gebracht! :
    // tiles und gekreuzt 3x3. inline + sonst optimierung
    if (!avx2Avail) ERR("method for matrix multiplication needs AVX2");
    AtAIntSpeedBlock3x3(a, b, nrow, ncol, ld, ldC, C, cores, m, n);
    break;

  default :
    // mode 0...4:, nur bestenfalls AVX, aber ohne tiling
    // 1 : AVX
    // 100 : primitiv
    // 101 : primitiv pipe/SIMD triggering
    if (cores == 0)  BUG;
    AtAInt(a, b, nrow, ncol, ld, ldC, C, cores, options->AtAmode);
  }
}



void AtAInt(int *a, 
	    Long nrow, Long ncol, // a and b have same size
	    Long *C, // result
	    int VARIABLE_IS_NOT_USED cores,
	    solve_options *options) {
  AtAInt(a, a, nrow, ncol, nrow, ncol, C, cores, options);
}


#define repet 4
#define MUPDATE(NAME,SIGN)						\
  void NAME##_##Long(Long *A, Long nrowA, Long ncolA, Long ldA, \
		     Long ldB, Long *B) {			\
    Long *a = A, *b=B;							\
    Long endA = nrowA - repet + 1;					\
    for (Long j=0; j<ncolA; j++, b+=ldB, a+=ldA) {	\
      Long i = 0;							\
      for ( ; i<endA; i+=repet) {					\
	a[i] SIGN##= b[i];						\
	a[i+1] SIGN##= b[i+1];						\
	a[i+2] SIGN##= b[i+2];						\
	a[i+3] SIGN##= b[i+3];						\
	/*a[i+4]SIGN=b[i+4]; a[i+5] SIGN= b[i+5];a[i+6] SIGN= b[i+6]; a[i+7] SIGN= b[i+7];*/ \
      }									\
      for (; i<nrowA; i++) a[i] SIGN##= b[i];				\
    }									\
  }
MUPDATE(Minus,-)
MUPDATE(Plus,+)
MUPDATE(Equal,)


#define MOP(NAME, SIGN, INT)						\
  void NAME##_##INT(INT *A, Long nrowA, Long ncolA, Long ldA,	\
		    INT *B, Long nrowB, Long ncolB, Long ldB,\
		  Long ldC, INT *C) {		\
/* C must be max(nrowA, nrowB) * max(ncolA, ncolB) */\
  Long ncolR = ncolA < ncolB ? ncolB : ncolA;\
  Long nrowR = nrowA < nrowB ? nrowB : nrowA;\
  Long bytes = sizeof(INT) * nrowA;\
  MEMSET(C, 0, ncolR * nrowR * sizeof(INT));\
  INT  *r = C;\
  INT *a = A,\
    *b = B;\
  for (Long j=0; j<ncolA; j++, a+=ldA, r+=ldC)\
    MEMCOPY(r, a, bytes);\
  if (b == NULL) return;			\
						\
  r = C;					\
  Long endB = nrowB - repet +1;\
  for (Long j=0; j<ncolB; j++, b+=ldB, r+=ldC) {\
    Long i = 0;\
    for ( ; i<endB; i+=repet) {\
      r[i] SIGN##= b[i];\
      r[i+1] SIGN##= b[i+1];\
      r[i+2] SIGN##= b[i+2];\
      r[i+3] SIGN##= b[i+3];\
      /*r[i+4]SIGN##=b[i+4]; r[i+5] SIGN##= b[i+5];r[i+6] SIGN##= b[i+6]; r[i+7] SIGN##= b[i+7];*/\
    }\
    for (; i<nrowB; i++) r[i] SIGN##= b[i];\
  }\
}

MOP(Minus,-,int)
MOP(Plus,+,int)



#define a11 A
#define a12 A + nrow1
#define a21 A + nrow1 *  ld
#define a22 A + nrow1 + ld * ncol1

#define A11 a11, ncol1, ncol1, ld
#define A12 a12, nrow2, ncol1, ld
#define A21 a21, nrow1, ncol2, ld
#define A22 a22, nrow2, ncol2, ld

#define b11 B
#define b12 B + ld * nrow1
#define b21 B + nrow1
#define b22 B + nrow1 + ld * ncol1
#define B11 b11, nrow1, ncol1, ld
#define B12 b12, nrow1, ncol2, ld
#define B21 b21, nrow2, ncol1, ld
#define B22 b22, nrow2, ncol2, ld

#define c11 C
#define c12 C + ldC * ncol1
#define c21 C + ncol1
#define c22 C + ncol1 + ldC * ncol1


#define C11 c11, ncol1, ncol1, ldC
#define C12 c12, ncol1, ncol2, ldC
#define C21 c21, ncol2, ncol1, ldC
#define C22 c22, ncol2, ncol2, ldC

#define CC11 dummyNeu, LenNeu, ldC, c11
#define CC12 dummyNeu, LenNeu, ldC, c12
#define CC21 dummyNeu, LenNeu, ldC, c21
#define CC22 dummyNeu, LenNeu, ldC, c22

#define NONE NULL, nrow1, ncol1, ld

#define DUMMY_DO nrow1, dummyNeu
#define DUMMY dummyNeu + sizeR, LenNeu - sizeR, DUMMY_DO /* memory for saving */

#define CHECK_DUMMY				\
  if (sizeR > LenNeu) ERR3(msg, __LINE__, sizeR, LenNeu )	\

#define CALC(NAME, X, OP, Y)					\
  /* NROW == ldC (NAME) */				\
  int *NAME = (int*) dummyNeu;					\
  CHECK_DUMMY;							\
  LenNeu -= sizeR;					\
  dummyNeu += sizeR;					\
  OP##_int(X, Y, nrow1, NAME);


#define DO(C, OP, A) OP##_Long(C, A)

#define STRASSEN(M1, M2, TRUEROW, TRUECOL,  RES) \
  /*  //printf("%d%s ", level,#M1); */					\
  strassen(M1, M2, TRUEROW, TRUECOL, ncol1, cores, options,  RES, level + 1)


#define showC					\
	/* // */  printf("\n");			\
  for (Long ii=0; ii<trueNcol; ii++) {		\
    for (Long jj=0; jj<trueNcol; jj++) {	\
    	/* // */  printf("%ld ", C[ii + jj * ldC]);		\
    }						\
  /* // */  printf("\n");					\
  }

#define showM(NR)							\
  /* // */ printf("*** M%d ***\n", NR);					\
  for (Long ii=0; ii<nrow1; ii++) {					\
    Long *where = NR == 7 ? c11 : NR == 6 ? c22 : dummyNeu;		\
    Long LD = NR == 7 ? ldC : NR == 6 ? ldC : nrow1;			\
  for (Long jj=0; jj<ncol1; jj++) {					\
    /* // */printf("%ld ", where[ii + jj * LD]);			\
  }									\
 /* // */ printf("\n");							\
  }  \
  /* // */  printf("\n");

#undef showC
#define showC
#undef showR
#define showR(NR)


void strassen(int *A, int *B, Long trueNrow, Long trueNcol,
	      Long ld, int cores, solve_options *options,
	      Long *dummy, Long Lendummy, 
	      Long ldC, Long *C, int level) {
  ASSERT_SOLVE(options);
  if (trueNrow <= options->StrassenMin || trueNcol <= options->StrassenMin) {
    solve_options opt;
    MEMCOPY(&opt, options, sizeof(solve_options));
    opt.AtAmode = 1;
    // opt.AtAmode = 12;
   //    printf(". . ");
    //printf("strassen->AtA=%ld %ld %ld %ld dummyLen=%ld\n", trueNrow, trueNcol, ld, ldC, Lendummy);
    AtAInt(A, B, trueNrow, trueNcol, ld, ldC, C, cores, &opt);
   return;
  }
  
 Long
    *dummyNeu, LenNeu;
  const Long
    nrow2 = trueNrow / 2,
    nrow1 = trueNrow - nrow2,
    ncol2 = trueNcol / 2,
    ncol1 = trueNcol - ncol2,
    sizeR = nrow1 * ncol1;
  const char* msg = "dummy too small in line %d (sizeR=%ld, Avail=%ld)";
  
  //  printf("\nstrassen= %ld x %ld; ld=%ld ldC=%ld dummyLen=%ld %ld level=%d ", trueNrow, trueNcol, ld, ldC, Lendummy, sizeR, level);

  //// M7
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M71, A12, Minus, A22);
  CALC(M72, B21, Plus, B22);
  STRASSEN(M71, M72,
	   nrow2, ncol1,  // true, Calc.nrow, Calc.ncol
	   CC11);
  //  showM(7);
   //  printf("Attt\n"); showC;

  //// M6
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M61, A21, Minus, A11);
  CALC(M62, B12, Plus, B11);
  STRASSEN(M61, M62,
	   nrow1, ncol2 /* !! */ ,  // true
	   CC22);
  //  showM(6)
 
 
  //// M1
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M11, A22, Plus, A11);
  CALC(M12, B22, Plus, B11);
  CHECK_DUMMY; // memory for saving available ?
  STRASSEN(M11, M12,
	   nrow1, ncol1,  // true
	   DUMMY); // das hier geht schief!!!
  DO(C11, Plus, DUMMY_DO);
  DO(C22, Plus, DUMMY_DO);
  // showM(1) //(ok)
 
 
  //// M2
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M21, A21, Plus, A22);
  CALC(M22, B11, Minus, NONE);
  CHECK_DUMMY; // memory for saving available ?
  STRASSEN(M21, M22,
	   nrow1, ncol1 ,  // true
	   DUMMY);
  DO(C21, Equal, DUMMY_DO);
  DO(C22, Minus, DUMMY_DO);

  
  //// M3
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M31, A11, Minus, NONE);
  CALC(M32, B12, Minus, B22);
  CHECK_DUMMY; // memory for saving available ?
  STRASSEN(M31, M32,
	   nrow1, ncol1,  // true
	   DUMMY);
  DO(C12, Equal, DUMMY_DO);
  DO(C22, Plus, DUMMY_DO);
  //
  //  showM(3)
  //
  //  showC;

  //// M4
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M41, A22, Plus, NONE);
  CALC(M42, B21, Minus, B11);
  CHECK_DUMMY; // memory for saving available ?
 STRASSEN(M41, M42,
	   nrow1, ncol1,  // true
	   DUMMY);
  DO(C11, Plus, DUMMY_DO);
  DO(C21, Plus, DUMMY_DO);
  //  showM(4); 
 
 
  //// M5
  dummyNeu = dummy;
  LenNeu = Lendummy;
  CALC(M51, A11, Plus, A12);
  CALC(M52, B22, Minus, NONE);
  CHECK_DUMMY; // memory for saving available ?
  STRASSEN(M51, M52,
	   nrow1, ncol1,  // true
	   DUMMY);
  DO(C11, Minus, DUMMY_DO);
  DO(C12, Plus, DUMMY_DO);
  //  showM(5); 
   
  //showC;
  return;
  
}


void Strassen(int *A, int  *B,
	      Long nrow, Long ncol, // a and b have same size 
	      Long ld,
	      Long ldC,
	      Long *C,
	      int VARIABLE_IS_NOT_USED cores, solve_options *options) {
  // calculates A^T * B with same sizes of A and B
  ASSERT_SOLVE(options);
  const double factor = options->StrassenFactor;
  Long 
    trueNrow = nrow,
    trueNcol = ncol;
  
  if (factor * (double) nrow <= (double) ncol) trueNcol = nrow;
  else if (factor * (double) ncol <=  (double) nrow) trueNrow = ncol;

  Long
    sizeP2 = 1L << ((int) (LOG((double) (trueNrow > trueNcol ? trueNrow
					 : trueNcol)) / LOG2 - 1e-10) + 1L),
    dummy_size = sizeP2  * sizeP2 + 00,
    *dummy = (Long*) CALLOC(sizeof(Long), dummy_size);

  //  if (sizeP2 != (trueNrow > trueNcol ? trueNrow : trueNcol))
  //  ERR("currently only power of 2 available");

  strassen(A, B, trueNrow, trueNcol, ld, cores, options,
	   dummy, dummy_size, ldC, C, 0);
  if (trueNrow != nrow && trueNcol == ncol) {
    Long *intermediate = (Long*) MALLOC(sizeof(Long) * ncol * ncol);
    int *a = A+trueNrow,
      *b = B+trueNrow;
    
    for (Long i=trueNrow; i<nrow; i+=trueNrow, a+=trueNrow, b+=trueNrow) {      
      strassen(a, b, 
	       i <= nrow - trueNrow ? trueNrow : nrow % trueNrow, 
	       trueNcol,
	       ld, cores, options,
	       dummy, dummy_size, ncol, intermediate, 0);
      Plus_Long(C, ncol, ncol, ldC, ncol, intermediate);
      }
    
    FREE(intermediate);
    
  } else if (trueNrow == nrow && trueNcol != ncol) {
    ERR("Strassen only programmed for inner lengths greater than outer");
  } 

  FREE(dummy);
   //showC; BUG;
}
  


#define algn_general(X)  ((1U + (uintptr_t) (((uintptr_t) X - 1U) / BytesPerBlock)) * BytesPerBlock)

double static inline *algn(double *X) {
  assert(algn_general(X)>=(uintptr_t)X); return (double *) algn_general(X);
}

#if defined SSE41 || defined AVX2
int static inline *algnInt(int *X) {
  assert(algn_general(X)>=(uintptr_t)X); return (int *) algn_general(X);
}
#endif



void colMaxsI256(double *M, Long r, Long c, double *ans, int cores);
void colMaxsIint256(int *M, Long r, Long c, int *ans, int cores);


void colMaxsIint(int *M, Long r, Long c, int *ans, int cores) {
  if (r < 32
#if defined AVX2
      || !avx2Avail
#elif defined  SSE41
      || !sse41Avail
#endif      
       ) {
    for (Long i=0; i<c; i++) {
      int *m = M + r * i,
	dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
      ans[i] = dummy;
    }
    return;
  }

  if (avx2Avail) {
    colMaxsIint256(M, r, c, ans, cores);
    return;
  }
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif    
  for (Long i=0; i<c; i++) {
     int dummy,
      *m = M + r * i;
#if defined SSE41 || defined AVX2    
     int *start = algnInt(m),
       *end = m + r;
    uintptr_t End = (uintptr_t) (end - integers);
    if ((uintptr_t) start < End) {
      BlockType *m0 = (BlockType0*) start,
	Dummy = LOAD((BlockType0*) m0);
      for (m0++ ; (uintptr_t) m0 < End; m0++) {
	Dummy = MAXINTEGER(Dummy, LOAD(m0));
      }
      int *d = (int *) &Dummy;
      dummy = d[0];
      dummy = MAX(dummy, d[1]);
      dummy = MAX(dummy, d[2]);
      dummy = MAX(dummy, d[3]);
#if defined AVX2
      dummy = MAX(dummy, d[4]);
      dummy = MAX(dummy, d[5]);
      dummy = MAX(dummy, d[6]);
      dummy = MAX(dummy, d[7]);
#endif // AVX
      for ( ; m<start; m++) dummy = MAX(dummy, *m);
      m = (int *) m0;
      for ( ; m<end; m++) dummy = MAX(dummy, *m);
    } else {
      dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
    }
#else // not SSE4
    dummy = m[0];    
    for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
#endif    
    ans[i] = dummy;
  }
}

 
      
void colMaxsI(double *M, Long r, Long c, double *ans, int cores) {
  if (r < 16
#if defined AVX2
      || !avx2Avail
#elif defined  SSE2
      || !sse2Avail
#endif      
      ) {
    for (Long i=0; i<c; i++) {
      double *m = M + r * i,
	dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
      ans[i] = dummy;
    }
    return;
  }

  if (avx2Avail) {
    colMaxsI256(M, r, c, ans, cores);
    return;
  }
 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif  
  for (Long i=0; i<c; i++) {
    double dummy,
      *m = M + r * i;
#if defined SSE2
    double *start = algn(m),
      *end = m + r;
    uintptr_t End = (uintptr_t) (end - doubles);
    if ((uintptr_t) start < End) {
      Doubles * m0 = (Doubles*) start,
	Dummy = (Doubles) LOAD((BlockType0*) m0);
      for (m0++ ; (uintptr_t) m0 < End; m0++) {
	Dummy = MAXDOUBLE(Dummy, (Doubles) LOAD((BlockType0*) m0));
      }
      double *d = (double *) &Dummy;
      dummy = d[0];
      dummy = MAX(dummy, d[1]);
#if defined AVX
      dummy = MAX(dummy, d[2]);
      dummy = MAX(dummy, d[3]);
#endif
      for ( ; m<start; m++) dummy = MAX(dummy, *m);
      m = (double *) m0;
      for ( ; m<end; m++) dummy = MAX(dummy, *m);
    } else {
      dummy = m[0];    
      for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
    }
#else
    dummy = m[0];    
    for (Long j=1; j<r; j++) dummy = MAX(dummy, m[j]);
#endif    
    ans[i] = dummy;
  }
}

void rowProdI(double *m, Long r, Long c, double *ans) {
  Long r4 = r / 4;
  MEMCOPY(ans, m, sizeof(*ans) * r);
  m += r;
  for (Long ic=1; ic<c; ic++) {
    double *a = ans;
    for (Long ir=0; ir<r4; ir++) {
      *(a++) *= *(m++);
      *(a++) *= *(m++);
      *(a++) *= *(m++);
      *(a++) *= *(m++);
    }
    for (Long ir=r4 * 4; ir<r; ir++) *(a++) *= *(m++);
  }
}



void rowMeansI(void *M, int sexp, Long r, Long c, double *weight, double *ans) {
  for (Long j=0; j<r; j++) ans[j] = 0.0;
  if (weight == NULL) {    
#define for1					\
    for (Long i=0; i<c; i++, m+=r) {			\
      for (Long j=0; j<r; j++) ans[j] += (double) m[j];	\
    }
  
    if (sexp== REALSXP) { double *m = (double*) M; for1; }
    else {
      int *m = (int*) M;
      for1;
    }    
  } else {    
 #define for2							\
    for (Long i=0; i<c; i++, m+=r) {				\
      double dummy = weight[i]; /* load1(weight); MULTDOUBLE */ \
      for (Long j=0; j<r; j++) ans[j] += (double) m[j] * dummy;	\
    }

    if (sexp == REALSXP) { double *m = (double*) M; for2; }
    else {
      int *m = (int*) M;
      for2;
    }
    
   }
  double invc = 1.0 / (double) c;
  for (Long j=0; j<r; j++) ans[j] *= invc;
}

void dbinormI(double *x, double *y, Long nrow, double *sigma, double *ans) {  
  //  Long nrow4 = nrow - 4;
  if (sigma == NULL) {
    double invtwopi = 1.0 / TWOPI;
    /*
      minushalfX[4] ={-0.5, -0.5, -0.5, -0.5},
      invtwopiX [4] = {invtwopi, invtwopi, invtwopi, invtwopi};
      Long i=0;
      
      #define atonce 4
    __m256d minushalf4 = LOADuDOUBLE(minushalfX),
       invtwopi4 = LOADuDOUBLE(invtwopiX);
      
    for (; i<nrow4; i+=atonce) {
      __m256d x4 = LOADuDOUBLE(x + i);
      double *xx4 = (double *) &x4;
      x4 = MULTDOUBLE(x4, x4);
      {
	__m256d y4 = LOADuDOUBLE(y + i);
	y4 = MULTDOUBLE(y4, y4);
	x4 = ADDDOUBLE(x4, y4);
      }
      x4 = MULTDOUBLE(minushalf4, x4);
      xx4[0] = EXP(xx4[0]);
      xx4[1] = EXP(xx4[1]);
      xx4[2] = EXP(xx4[2]);
      xx4[3] = EXP(xx4[3]);
      x4 = MULTDOUBLE(x4, invtwopi4);
      STOREuDOUBLE(ans + i, x4);
    }
    */
    for (Long i=0; i<nrow; i++) 
      ans[i] = EXP(-0.5 * (x[i] * x[i] + y[i] * y[i])) * invtwopi;
  } else {
    double 
      sigma1 = sigma[0],
      sigma4 = sigma[3],
      inv2piSrtS = 1.0 / (TWOPI * SQRT(sigma1 * sigma4)),
      invS1half = 0.5 / sigma1,
      invS4half = 0.5 / sigma4;
    if (sigma[1] == 0.0 && sigma[2] == 0.0) {
      for (Long i=0 ; i<nrow; i++)
	ans[i] = EXP(- (x[i] * x[i] * invS1half + y[i] * y[i] * invS4half) )
	  * inv2piSrtS;
    } else BUG;
  }
 }




void dotXVI(double *M, Long r, Long c, double *V, double *Ans) {
// bringt nix
  //#ifdef DO_PARALLEL
  //#p ragma omp parallel for num_threads(cores) 
  //#endif  
  for (Long i=0; i<c; i++) {
    double 
      *ans = Ans + r * i,
      *v = V,
      *m = M + r * i;
#if defined SSE2_DONOTUSE_AS_SLOWER
    double *end = m + r - doubles;
    for ( ; m < end; m += doubles, ans += doubles, v += doubles)
      STOREuDOUBLE(ans, MULTDOUBLE(LOADuDOUBLE(m), LOADuDOUBLE(v)));
    end += doubles;
    for (; m < end; m++) *ans = *m * *v;
#else
    for (Long j=0; j<r; j++)  ans[j] = m[j] * v[j];
 #endif    
  }
}


Long scalarprodInt(int * v1, int * v2, Long N){
  int *endv1 = v1 + N;
  Long sum = 0;
  for(; v1!= endv1; v1++, v2++) sum +=  v2[0] * v1[0];
  return sum;
}
 
 
Long scalarprodInt2by2(int * v1, int * v2, Long N){
  int *endv1 = v1 + (N / 2) * 2,
    *end = v1 + N;
  Long sum = 0;
  for(; v1 < endv1; v1 += 2, v2 += 2) sum += v2[0] * v1[0] + v2[1] * v1[1];
  if (v1 < end) sum += v2[0] * v1[0]; 
  return sum;
}
 
 
Long scalarprodInt4by4(int * v1, int * v2, Long N){
  // printf("4by4 %d %d %d\n", sse, sse2, avx);
  int*endv1 = v1 + (N / 4) * 4,
    *end = v1 + N;
  Long sum = 0;
   for(; v1 < endv1; v1 += 4, v2 += 4)
    sum += v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2]+ v2[3] * v1[3];
  for(; v1 < end; v1++, v2++) sum += v2[0] * v1[0];        
  return sum;
}

 
Long scalarprodInt8by8(int * v1, int * v2, Long N){
  int
    *endv1 = v1 + (N / 8) * 8,
    *end = v1 + N;
  Long sum = 0;
  for(; v1 < endv1; v1 += 8, v2 += 8)
    sum += v2[0] * v1[0] + v2[1] * v1[1]+ v2[2] * v1[2] + v2[3] * v1[3] +
      v2[4] * v1[4] + v2[5] * v1[5]+ v2[6] * v1[6]+ v2[7] * v1[7];
  for(; v1 < end; v1++, v2++) sum +=  v2[0] * v1[0];        
  return sum;
}


Long scalarprodInt16by16(int * v1, int * v2, Long N){
  int
    *endv1 = v1 + (N / 16) * 16,
    *end = v1 + N;
  Long sum = 0;
  for(; v1 < endv1; v1 += 16, v2 += 16)
    sum += v2[0] * v1[0] + v2[1] * v1[1]+ v2[2] * v1[2] + v2[3] * v1[3] +
      v2[4] * v1[4] + v2[5] * v1[5]+ v2[6] * v1[6]+ v2[7] * v1[7] +
      v2[8] * v1[8] + v2[9] * v1[9]+ v2[10] * v1[10] + v2[11] * v1[11] +
      v2[12] * v1[12] + v2[13] * v1[13]+ v2[14] * v1[14]+ v2[15] * v1[15];
  for(; v1 < end; v1++, v2++) sum +=  v2[0] * v1[0];        
  return sum;
}


Long avx_scalarprodInt8(int * x, int * y, Long L);
Long avx_scalarprodInt2(int * x, int * y, Long L);
Long scalarXint(int *x, int *y, Long len, Long n) {
  // parallel lohnt i.A. nicht: 28.2.20121 alles was parallel ist, rausgeworfen
  assert(n >= 0);
  //  __m128 a, b;  a = _mm_add_ps ((__m128) a, (__m128) b);
  //   __m128i c, d; c = _mm_add_epi16 ((__m128i) c, (__m128i) d);
  //  __m256d e, f;   e = _mm256_add_pd (e,  f);
   //  printf("n=%d %d ",n, 0);
  switch(n) {
    //  printf("%d\n", n);
  case SCALAR_AVX : // == 1
    //   printf(" # %d ", avx);
    if (avx2Avail) return avx_scalarprodInt8(x, y, len); // best one kernel
    break;
  case 2 : return scalarprodInt(x, y, len);
  case 3 : return scalarprodInt2by2(x, y, len); 
  case 4  : return scalarprodInt4by4(x, y, len); 
  case 5 : return scalarprodInt8by8(x, y, len); 
  case 6 : return scalarprodInt16by16(x, y, len); 
  case 8 :
    if (avx2Avail) return avx_scalarprodInt2(x, y, len);  //best
    break;
 
 
  case SCALAR_BASE :
  default :
    {}
  }
  //  return scalarprodInt(x, y, len);
  return scalarprodInt4by4(x, y, len);
}
  



double scalarprod(double * v1, double * v2, Long N){
  double *endv1 = v1 + N,
    sum = 0;
  for(; v1!= endv1; v1++, v2++) sum +=  v2[0] * v1[0];
  return sum;
}
 
 
double scalarprod2by2(double * v1, double * v2, Long N){
  double *endv1 = v1 + (N / 2) * 2,
    *end = v1 + N,
    sum = 0;
  for(; v1 < endv1; v1 += 2, v2 += 2) sum += v2[0] * v1[0] + v2[1] * v1[1];
  if (v1 < end) sum += v2[0] * v1[0]; 
  return sum;
}
 
 
double scalarprod4by4(double * v1, double * v2, Long N){
  // printf("4by4 %d %d %d\n", sse, sse2, avx);
  double*endv1 = v1 + (N / 4) * 4,
    *end = v1 + N,
    sum = 0;
  for(; v1 < endv1; v1 += 4, v2 += 4)
    sum += v2[0] * v1[0] + v2[1] * v1[1] + v2[2] * v1[2]+ v2[3] * v1[3];
  for(; v1 < end; v1++, v2++) sum += v2[0] * v1[0];        
  return sum;
}

 
double scalarprod8by8(double * v1, double * v2, Long N){
  double
    *endv1 = v1 + (N / 8) * 8,
    *end = v1 + N,
    sum = 0.0;
  for(; v1 < endv1; v1 += 8, v2 += 8)
    sum += v2[0] * v1[0] + v2[1] * v1[1]+ v2[2] * v1[2] + v2[3] * v1[3] +
      v2[4] * v1[4] + v2[5] * v1[5]+ v2[6] * v1[6]+ v2[7] * v1[7];
  for(; v1 < end; v1++, v2++) sum +=  v2[0] * v1[0];        
  return sum;
}


void avx_scalarprodM(double * x, double * y, Long len, double *res);
double avx_scalarprodDnearfma(double * x, double * y, Long len);
double avx_scalarprodD(double * x, double * y, Long L);
double avx_scalarprodDopt(double * x, double * y, Long L);
double avx_scalarprodDP(double * x, double * y, Long L) ;
double avx_scalarprodDK(double * x, double * y, Long L);
double scalarX(double *x, double *y, Long len, Long n) {
  // parallel lohnt i.A. nicht: 28.2.20121 alles was parallel ist, rausgeworfen
  assert(n >= 0);
  //  __m128 a, b;  a = _mm_add_ps ((__m128) a, (__m128) b);
  //   __m128i c, d; c = _mm_add_epi16 ((__m128i) c, (__m128i) d);
  //  __m256d e, f;   e = _mm256_add_pd (e,  f);
   //  printf("n=%d %d ",n, avx);
  switch(n) {
    //  printf("%d\n", n);
  case SCALAR_AVX :
    //  printf(" # %d ", avx);
    if (avxAvail) return avx_scalarprodD(x, y, len); // best one kernel
    break;
  case 2 : return scalarprod(x, y, len);
  case 3 : return scalarprod2by2(x, y, len); 
  case 4 : return scalarprod4by4(x, y, len); 
  case 5 : return scalarprod8by8(x, y, len); 
    //  case 5 :
    //#ifdef FMA_AVAILABLE
    //   return avx_scalarprodDfma(x, y, len);
    //#endif    
  case SCALAR_NEARFMA :
    if (avxAvail) return avx_scalarprodDnearfma(x, y, len);
    break;
  case 8 :
    if (avxAvail) return avx_scalarprodDP(x, y, len);  //best
    break;
  case SCALAR_KAHAN :
    if (avxAvail) return avx_scalarprodDK(x, y, len); // kahan   
    break;

  case SCALAR_BASE :
  default :
    {}
  }
  return scalarprod4by4(x, y, len);
}
  

void avx_linearprodD(double * v1,  double v2, Long N, double *inout);
void linearprod2by2(double * v1,  double v2, Long N, double *inout){
  double *endv1 = v1 + (N / 2) * 2,
    *end = v1 + N;
  for(; v1 < endv1; v1+=2, inout+=2) {
      inout[0] += v2 * v1[0];
      inout[1] += v2 * v1[1];
  }
  if (v1 < end) inout[0] += v2 * v1[0];
}
 

void linearX(double *x, double y, Long len, double *inout, Long n) {
  switch(n) {
  case LINEAR_AVX :
    if (avxAvail) { avx_linearprodD(x, y, len, inout); return; }
    break; 
  case LINEAR_BASE:
  default :
    {}
  }
  linearprod2by2(x, y, len, inout); 
}
  
