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
#include "compatibility.lapack.h"
#include "kleinkram.h"
#include "options_RFU.h"
#include "utils.h"
//#include "xport_import_RFU.h"
#include "extern_RFU.h"


#if defined AVX2


ASSERT_SIMD(avx2_fctns, avx2);


#define repet 8
#define atonce (integers * repet)
#define SET_0(NR) sum##NR = _mm256_setzero_si256()
#define P_0(NR) prod##NR = _mm256_setzero_si256()


#define ADDI(NR)							\
  prod0 = MULT32(LOADU( (__m256i*) (x + i + NR * integers)),		\
    LOADU( (__m256i*) (y + i + NR * integers)));			\
  sum0 = ADD32(sum0, prod0)


Long avx_scalarprodInt8(int * x, int * y, Long len){
  BUG; // es muss doch integers / 2 nicht integers sein!!!!
  Long i = 0,
    sum = 0;
  __m256i SET_0(0), P_0(0);

  // BUG;
  if ( len >= atonce) {
    if (!false && ((uintptr_t) x - (uintptr_t) ( y)) % BytesPerBlock == 0) {
      // bringt leider nicht viel auf avx2 bei Skalarprodukten!!
      // rest <- 
      Long odd = (integers - ((uintptr_t) x % BytesPerBlock) / sizeof(int))
	% integers;

       for (; i < odd; i++, x++, y++) sum += *x * *y;
      len -= odd;
      if (odd != 0) { printf("odd=%ld %ld\n", odd, (uintptr_t) x % BytesPerBlock);  }

      const Long lenM = len - (atonce - 1);  
      for (; i < lenM; i += atonce) {
	ADDI(0); ADDI(1); ADDI(2); ADDI(3);
	ADDI(4); ADDI(5); ADDI(6); ADDI(7); 
      }
    } else {
      const Long lenM = len - (atonce - 1);  
      for (; i < lenM; i += atonce) {
	ADDI(0); ADDI(1); ADDI(2); ADDI(3);
	ADDI(4); ADDI(5); ADDI(6); ADDI(7); 
      }
    }
  }  
  const Long end = len - integers + 1;
  for (; i < end; i += integers) { ADDI(0); } 
  int *D  = (int *) &sum0;
  sum += (Long) D[0] + D[1] + D[2] + D[3] + D[4] + D[5] + D[6] + D[7];
  for (; i < len; i++) sum += x[i] * y[i];
  return sum;
}



Long avx_scalarprodInt8X(int * x, int * y, Long len){
  BUG; // es muss doch repet = 4 sein!!!!
  Long i = 0,
    lenM = len - (atonce - 1);  
  __m256i SET_0(0), P_0(0);
 
  if ( len >= atonce) {
    for (; i < lenM; i += atonce) {
      ADDI(0); ADDI(1); ADDI(2); ADDI(3); ADDI(4); ADDI(5); ADDI(6); ADDI(7); 
    }
  }
  lenM = len - integers + 1;
  for (; i < lenM; i += integers) { ADDI(0); } 
  int *D  = (int *) &sum0;
  Long sum = (Long) D[0] + D[1] + D[2] + D[3] + D[4] + D[5] + D[6] + D[7];
   for (; i < len; i++) sum += x[i] * y[i];
  return sum;
}

#undef repet
#define repet 4
Long avx_scalarprodInt2(int * x, int * y, Long len){
  Long i = 0,
    lenM = len - (atonce - 1);  
    
  __m256i SET_0(0), SET_0(1), SET_0(2), SET_0(3),
    *X = (__m256i*) x,
    *Y = (__m256i*) y;

  if ( len >= atonce) {
    for (; i < lenM; i += atonce) {
      const __m256i a0 = LOADU(X++);
      const __m256i b0 = LOADU(Y++);
      const __m256i a1 = LOADU(X++);
      const __m256i b1 = LOADU(Y++);
      const __m256i a2 = LOADU(X++);
      const __m256i b2 = LOADU(Y++);
      const __m256i a3 = LOADU(X++);      
      const __m256i b3 = LOADU(Y++);
      const __m256i c0 = MULT32(a0, b0);
      const __m256i c1 = MULT32(a1, b1);
      const __m256i c2 = MULT32(a2, b2);
      const __m256i c3 = MULT32(a3, b3);
      sum0 = ADD32(c0, sum0);
      sum1 = ADD32(c1, sum1);
      sum2 = ADD32(c2, sum2);
      sum3 = ADD32(c3, sum3);
     }
  }
  lenM = len - integers + 1;
  __m256i P_0(0);
  for (; i < lenM; i += integers) { ADDI(0); }
  int *D;
  Long sum = 0;
  D = (int *) &sum0; sum += (Long) D[0]+D[1]+D[2]+D[3]+D[4]+D[5]+ D[6]+D[7];
  D = (int *) &sum1; sum += (Long) D[0]+D[1]+D[2]+D[3]+D[4]+D[5]+ D[6]+D[7];
  D = (int *) &sum2; sum += (Long) D[0]+D[1]+D[2]+D[3]+D[4]+D[5]+ D[6]+D[7];
  D = (int *) &sum3; sum += (Long) D[0]+D[1]+D[2]+D[3]+D[4]+D[5]+ D[6]+D[7];
    
  for (; i < len; i++) sum += x[i] * y[i];
  return sum;
}



#undef repet
#define repet 1
void scalarXintBlock2x2(int * x, int * y, Long len, Long ld, Long ncol,
		     Long *C){
  BUG; // FEHLER?!: MULT32 liest auch nur 128 BIT !!
  Long i = 0,
    lenM = len - (atonce - 1);  
  __m256i
    SET_0(00), SET_0(01), 
    SET_0(10), SET_0(11), 
   *X0 = (__m256i*) x,
    *Y0 = (__m256i*) y,
    *X1 = (__m256i*) (x + ld),
    *Y1 = (__m256i*) (y + ld)
    ;

  for (; i < lenM; i += atonce) {
    const __m256i a0 = LOADU(X0++);
    const __m256i b0 = LOADU(Y0++);
    const __m256i a1 = LOADU(X1++);
    const __m256i b1 = LOADU(Y1++);
   
    const __m256i c00 = MULT32(a0, b0);
    const __m256i c01 = MULT32(a0, b1);
    const __m256i c10 = MULT32(a1, b0);
    const __m256i c11 = MULT32(a1, b1);
    
    sum00 = ADD32(c00, sum00);
    sum01 = ADD32(c01, sum01);
    sum10 = ADD32(c10, sum10);
    sum11 = ADD32(c11, sum11);
  }
  
  int *D, *X, *Y;
  Long sum = 0;

#define ADDUP(I, J)							\
  D = (int*) &sum##I##J; sum = (Long) D[0]+D[1]+D[2]+D[3]+D[4]+D[5]+D[6]+D[7];\
  /*  printf("ADD=%d %d; %d%d D=%ld \n", i, len, I, J, sum); */		\
  X = x + (I) * ld;					\
  Y = y + (J) * ld;					\
  for (Long i0=i; i0 < len; i0++) sum += X[i0] * Y[i0];			\
  C[(I) + (J) * ncol] += sum;						\

  ADDUP(0,0); 
  ADDUP(0,1); 
  ADDUP(1,0); 
  ADDUP(1,1); 
   
  return;
}




void AtAIntBlock2x2(int *a, int *b,
		    Long nrow, Long ncol, // a and b have same size
		    Long ld, Long ldC, Long *C, // result
		    int VARIABLE_IS_NOT_USED cores,
		    int scalarVersion,		    
		    Long m, Long n) { // tile sizes

  // C =  A^T %*% B
 Long
   bytes = sizeof(*C) * ncol,
    tileRows = m <= 0 || m > nrow ? nrow : m,// for rowPos
   tileCols = n <= 0 || n > ncol ? ncol : n; // of result
  Long miniRows = 2;
  Long miniCols = 2;
  
  if (ldC == ncol) MEMSET(C, 0,  bytes * ncol);
  else for (Long i=0; i<ncol; i++) MEMSET(C + ldC * i, 0, bytes);

 #ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif
  for (Long tR=0; tR<ncol; tR += tileCols) {
    //    printf("tR = %ld %ld %d >= %ld\n", tR, tileCols, miniRows, ncol);
    Long Rend = tR + tileCols < ncol ? tR + tileCols : ncol,
      i_end = Rend - miniRows + 1;      
    for (Long tC=a==b ? tR : 0; tC<ncol; tC += tileCols) {
      Long Cend = tC + tileCols < ncol ? tC + tileCols : ncol,
	j_end = Cend - miniCols + 1;	  
      for (Long rowPos = 0; rowPos < nrow; rowPos += tileRows) {
	Long curLen = tileRows < nrow - rowPos? tileRows : nrow - rowPos;
	Long i = tR;
	for ( ; i<i_end; i+= miniRows) {
	  Long j = a==b && tC == tR ? i : tC;
	  int	    
	    *A = a + i * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<j_end; j+=miniCols, B+=miniCols * ld) {
	    if (j > ncol - miniCols || i > ncol - miniRows ||
		curLen  > nrow - rowPos) {
	      printf("j=%ld >= %ld - %ld; i=%ld >= %ld - %ld; Len=%ld > %ld - %ld\n",
		     j, ncol, miniCols,
		     i, ncol, miniRows,
		     curLen, nrow, rowPos);
	      ERR("seg-fault verhindert");
	    }
	    scalarXintBlock2x2(A, B, curLen, ld, ldC,
			    C + i + ldC * j);
	    // printf("i=%ld,%ld; nrow=%ld,%ld start=%ld, tR=%ld,%ld C=%ld; %ld=%ld+%ld*%ld\n", i, j, nrow, ncol, start, tR, tC, C [i + ldC * j], i + ldC * j, i, ldC, j);
	  }

	  //	  if (false)
	  if (j < Cend) {
	    for (Long ii=0 ; ii<miniRows; ii++) {
	      Long jj = j;
	      A = a + (i+ii) * ld + rowPos;
	      B = b + jj * ld + rowPos;
	      for ( ; jj<Cend; jj++, B+=ld) 
		C[i + ii + ldC * jj] +=
		  scalarXint(A, B, curLen, scalarVersion);
	    }
	  }
	}

	//	if (false)
	for (Long ii=i; ii<Rend; ii++) {
	  Long j = a==b && tC == tR ? ii : tC;
	  int
	    *A = a + ii * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<Cend; j++, B+=ld) 
	    C[ii + ldC * j] += scalarXint(A, B, curLen, scalarVersion);
	}
	
      }
    }
  }
  
  if (a==b) {
    for (Long i=0; i<ncol; i++) 
      for (Long j=i; j<ncol; j++) 
	C[j + ldC * i] = C[j * ldC + i];
  }


  //  for (Long i=0; i<ncol; i++) {
  //    for (Long j=0; j<ncol; j++)
  //      printf("%d:%d:%d ", (int) (j * ldC + i), (int) C[j * ldC + i],(int) C[2]);
  //    printf("\n");
  //  }

}
    
void scalarXintBlock3x3(int * x, int * y, Long len, Long ld, Long ncol,
		     Long *C){
    BUG; // FEHLER?!: MULT32 liest auch nur 128 BIT !!

    
  Long i = 0,
    lenM = len - (atonce - 1);  
   
  __m256i
    SET_0(00), SET_0(01), SET_0(02),
    SET_0(10), SET_0(11), SET_0(12),
    SET_0(20), SET_0(21), SET_0(22),
    *X0 = (__m256i*) x,
    *Y0 = (__m256i*) y,
    *X1 = (__m256i*) (x + ld),
    *Y1 = (__m256i*) (y + ld),
    *X2 = (__m256i*) (x + 2 * ld),
    *Y2 = (__m256i*) (y + 2 * ld)
    ;


  //  if (false)  // crazy
  for (; i < lenM; i += atonce) {
    const __m256i a0 = LOADU(X0++);
    const __m256i b0 = LOADU(Y0++);
    const __m256i a1 = LOADU(X1++);
    const __m256i b1 = LOADU(Y1++);
    const __m256i a2 = LOADU(X2++);
    const __m256i b2 = LOADU(Y2++);
    
    const __m256i c00 = MULT32(a0, b0);
    const __m256i c01 = MULT32(a0, b1);
    const __m256i c02 = MULT32(a0, b2);
    const __m256i c10 = MULT32(a1, b0);
    const __m256i c20 = MULT32(a2, b0);
    const __m256i c11 = MULT32(a1, b1);
    const __m256i c12 = MULT32(a1, b2);
    const __m256i c21 = MULT32(a2, b1);
    const __m256i c22 = MULT32(a2, b2);
    
    sum00 = ADD32(c00, sum00);
    sum01 = ADD32(c01, sum01);
    sum02 = ADD32(c02, sum02);
    sum10 = ADD32(c10, sum10);
    sum11 = ADD32(c11, sum11);
    sum12 = ADD32(c12, sum12);
    sum20 = ADD32(c20, sum20);
    sum21 = ADD32(c21, sum21);
    sum22 = ADD32(c22, sum22);
  }
  
  int *D, *X, *Y;
  Long sum = 0;

#define ADDUP(I, J)							\
  D = (int*) &sum##I##J; sum = (Long) D[0]+D[1]+D[2]+D[3]+D[4]+D[5]+D[6]+D[7];\
  /*  printf("ADD=%d %d; %d%d D=%ld \n", i, len, I, J, sum); */		\
  X = x + (I) * ld;					\
  Y = y + (J) * ld;					\
  for (Long i0=i; i0 < len; i0++) sum += X[i0] * Y[i0];			\
  C[(I) + (J) * ncol] += sum;						\
 
  ADDUP(0,0); 
  ADDUP(0,1); 
  ADDUP(0,2); 
  ADDUP(1,0); 
  ADDUP(1,1); 
  ADDUP(1,2); 
  ADDUP(2,0); 
  ADDUP(2,1); 
  ADDUP(2,2); 
  
  return;
}





void AtAIntBlock3x3(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		 int VARIABLE_IS_NOT_USED cores,
		    int scalarVersion,
		 Long m, Long n) { // tile sizes
 // C =  A^T %*% B
 Long
    bytes = sizeof(*C) * ncol,
   tileRows = m <= 0 || m > nrow ? nrow : m, // for rowPos
   tileCols = n <= 0 || n > ncol ? ncol : n; // of result
 Long miniRows = 3;
  Long miniCols = 3;
  
  if (ldC == ncol) MEMSET(C, 0,  bytes * ncol);
  else for (Long i=0; i<ncol; i++) MEMSET(C + ldC * i, 0, bytes);

 #ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif
  for (Long tR=0; tR<ncol; tR += tileCols) {
    //    printf("tR = %ld %ld %d >= %ld\n", tR, tileCols, miniRows, ncol);
    Long Rend = tR + tileCols < ncol ? tR + tileCols : ncol,
      i_end = Rend - miniRows + 1;      
    for (Long tC=a==b ? tR : 0; tC<ncol; tC += tileCols) {
      Long Cend = tC + tileCols < ncol ? tC + tileCols : ncol,
	j_end = Cend - miniCols + 1;	  
      for (Long rowPos = 0; rowPos < nrow; rowPos += tileRows) {
	Long curLen = tileRows < nrow - rowPos? tileRows : nrow - rowPos;
	Long i = tR;
	for ( ; i<i_end; i+= miniRows) {
	  Long j = a==b && tC == tR ? i : tC;
	  int	    
	    *A = a + i * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<j_end; j+=miniCols, B+=miniCols * ld) {
	    if (j > ncol - miniCols || i > ncol - miniRows ||
		curLen  > nrow - rowPos) {
	      printf("j=%ld >= %ld - %ld; i=%ld >= %ld - %ld; Len=%ld > %ld - %ld\n",
		     j, ncol, miniCols,
		     i, ncol, miniRows,
		     curLen, nrow, rowPos);
	      ERR("seg-fault");
	    }
	    scalarXintBlock3x3(A, B, curLen, ld, ldC,
			    C + i + ldC * j);
	    //printf("X %ld %ld\n", nrow,  C [ i + ldC * j]);
	    // printf("i=%ld,%ld; nrow=%ld,%ld start=%ld, tR=%ld,%ld C=%ld; %ld=%ld+%ld*%ld\n", i, j, nrow, ncol, start, tR, tC, C [i + ldC * j], i + ldC * j, i, ldC, j);
	  }

	  //	  if (false)
	  if (j < Cend) {
	    for (Long ii=0 ; ii<miniRows; ii++) {
	      Long jj = j;
	      A = a + (i+ii) * ld + rowPos;
	      B = b + jj * ld + rowPos;
	      for ( ; jj<Cend; jj++, B+=ld) {
		C[i + ii + ldC * jj] += scalarXint(A, B, curLen, scalarVersion);
		//	printf("Y %ld\n",  C [ i + ii + ldC * j]);
	      }
	    }
	  }
	}

	//	if (false)
	for (Long ii=i; ii<Rend; ii++) {
	  Long j = a==b && tC == tR ? ii : tC;
	  int
	    *A = a + ii * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<Cend; j++, B+=ld) {
	    C[ii + ldC * j] += scalarXint(A, B, curLen, scalarVersion);
	    //	printf("Z %ld\n",  C [+ ii + ldC * j]);
	  }
	}
	
      }
    }
  }
  

  if (a==b) {
    for (Long i=0; i<ncol; i++) 
      for (Long j=i; j<ncol; j++) 
	C[j + ldC * i] = C[j * ldC + i];
  }


  //  for (Long i=0; i<ncol; i++) {
  //    for (Long j=0; j<ncol; j++)
  //      printf("%d:%d:%d ", (int) (j * ldC + i), (int) C[j * ldC + i],(int) C[2]);
  //    printf("\n");
  //  }

  
}
    




#define algn_general(X)  ((1U + (uintptr_t) (((uintptr_t) X - 1U) / BytesPerBlock)) * BytesPerBlock)

#if defined SSE41 || defined AVX2
int static inline *algnInt(int *X) {
  assert(algn_general(X)>=(uintptr_t)X); return (int *) algn_general(X);
}
#endif


void colMaxsIint256(int *M, Long r, Long c, int *ans, int cores) {
  if (r < 32
#if defined AVX2
      || !avx2Avail
#elif defined  SSE41
      || !sse41Avail
#endif      
       ) {
    for (int i=0; i<c; i++) {
      int *m = M + r * i,
	dummy = m[0];    
      for (int j=1; j<r; j++) dummy = MAX(dummy, m[j]);
      ans[i] = dummy;
    }
    return;
  }
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif    
  for (int i=0; i<c; i++) {
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
#endif // AVX2
      for ( ; m<start; m++) dummy = MAX(dummy, *m);
      m = (int *) m0;
      for ( ; m<end; m++) dummy = MAX(dummy, *m);
    } else {
      dummy = m[0];    
      for (int j=1; j<r; j++) dummy = MAX(dummy, m[j]);
    }
#else // not SSE4
    dummy = m[0];    
    for (int j=1; j<r; j++) dummy = MAX(dummy, m[j]);
#endif    
    ans[i] = dummy;
  }
}




 
//
#define scalarXintLocal avx_scalarprodInt8
inline Long scalarXintLocaXl(int * x, int * y, Long len){
  BUG; // es muss doch integers / 2 nicht integers sein!!!!
  Long i = 0,
    lenM = len - (atonce - 1);  
  __m256i SET_0(0), P_0(0);
 
  if ( len >= atonce) {
    for (; i < lenM; i += atonce) {
      ADDI(0); ADDI(1); ADDI(2); ADDI(3); ADDI(4); ADDI(5); ADDI(6); ADDI(7); 
    }
  }
  lenM = len - integers + 1;
  for (; i < lenM; i += integers) { ADDI(0); } 
  int *D  = (int *) &sum0;
  Long sum = (Long) D[0] + D[1] + D[2] + D[3] + D[4] + D[5] + D[6] + D[7];
  for (; i < len; i++) sum += x[i] * y[i];
  return sum;
}


void AtAIntSpeedBlock3x3(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		 int VARIABLE_IS_NOT_USED cores,
		 Long m, Long n) { // tile sizes
  // C =  A^T %*% B
 Long
    bytes = sizeof(*C) * ncol,
   tileRows = m <= 0 || m > nrow ? nrow : m, // for rowPos
   tileCols = n <= 0 || n > ncol ? ncol : n; // of result
  Long miniRows = 3;
  Long miniCols = 3;
  
  if (ldC == ncol) MEMSET(C, 0,  bytes * ncol);
  else for (Long i=0; i<ncol; i++) MEMSET(C + ldC * i, 0, bytes);

 #ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif
  for (Long tR=0; tR<ncol; tR += tileCols) {
    Long Rend = tR + tileCols < ncol ? tR + tileCols : ncol,
      i_end = Rend - miniRows + 1;      
    for (Long tC=a==b ? tR : 0; tC<ncol; tC += tileCols) {
      Long Cend = tC + tileCols < ncol ? tC + tileCols : ncol,
	j_end = Cend - miniCols + 1;	  
      for (Long rowPos = 0; rowPos < nrow; rowPos += tileRows) {
	Long curLen = tileRows < nrow - rowPos? tileRows : nrow - rowPos;
	Long i = tR;
	for ( ; i<i_end; i+= miniRows) {
	  Long j = a==b && tC == tR ? i : tC;
	  int	    
	    *A = a + i * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<j_end; j+=miniCols, B+=miniCols * ld) {
	    scalarXintBlock3x3(A, B, curLen, ld, ldC, C + i + ldC * j);
	  }

	  if (j < Cend) {
	    for (Long ii=0 ; ii<miniRows; ii++) {
	      Long jj = j;
	      A = a + (i+ii) * ld + rowPos;
	      B = b + jj * ld + rowPos;
	      for ( ; jj<Cend; jj++, B+=ld) {
		C[i + ii + ldC * jj] += scalarXintLocal(A, B, curLen);
	      }
	    }
	  }
	}

	for (Long ii=i; ii<Rend; ii++) {
	  Long j = a==b && tC == tR ? ii : tC;
	  int
	    *A = a + ii * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<Cend; j++, B+=ld) {
	    C[ii + ldC * j] += scalarXintLocal(A, B, curLen);
	  }
	}
	
      }
    }
  }
  

  if (a==b) {
    for (Long i=0; i<ncol; i++) 
      for (Long j=i; j<ncol; j++) 
	C[j + ldC * i] = C[j * ldC + i];
  }


  //  for (Long i=0; i<ncol; i++) {
  //    for (Long j=0; j<ncol; j++)
  //      printf("%d:%d:%d ", (int) (j * ldC + i), (int) C[j * ldC + i],(int) C[2]);
  //    printf("\n");
  //  }

  
}
    



#else

void colMaxsIint(int *M, Long r, Long c, int *ans, int cores);
void colMaxsIint256(int *M, Long r, Long c, int *ans, int cores) {colMaxsIint(M, r, c, ans, cores); }


Long avx_scalarprodInt8(int * x, int * y, Long len){
  ERR("avx_scalarprodInt8:avx2 not available");
  return NA_LONG;
}
Long avx_scalarprodInt2(int * x, int * y, Long len){
  ERR("avx_scalarprodInt2:avx2 not available");
  return NA_LONG;
}

void scalarXintBlock2x2(int * x, int * y, Long len, Long ld, Long ncol,
		     Long *C){
   ERR("scalarXintBlock2x2:avx2 not available");
}
void scalarXintBlock3x3(int * x, int * y, Long len, Long ld, Long ncol,
		     Long *C){
  ERR("scalarXintBlock3x3:avx2 not available");
}

void AtAIntSpeedBlock3x3(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		 int VARIABLE_IS_NOT_USED cores,
		 Long m, Long n) { // tile sizes
  ERR("AtAIntSpeedBlock3x3: avx2 not available");
}

void AtAIntBlock2x2(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		 int VARIABLE_IS_NOT_USED cores, int scalarVersion,
		 Long m, Long n)  { // tile sizes
  ERR("AtAIntBlock2x2: avx2 not available");
}

void AtAIntBlock3x3(int *a, int *b,
		 Long nrow, Long ncol, // a and b have same size
		 Long ld, Long ldC, Long *C, // result
		    int VARIABLE_IS_NOT_USED cores, int scalarVersion,
		 Long m, Long n)  { // tile sizes
  ERR(":AtAIntBlock3x3 avx2 not available");
}



SIMD_MISS(avx2_fctns, avx2);

#endif

