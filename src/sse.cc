
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



//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// additional code for SSE2:

#include <smmintrin.h>
#include <inttypes.h>
#include "intrinsics.h"

#include "miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "AutoMiraculix.h"
#include <General_utils.h>
#include "general.relmatrix.h"

#include "error.h"

#define Method MethodSSE
#define Information InformationSSE
#if defined SSSE3 or defined SSE2
typedef union {
  //__m128 vf;
  __m128i vi;
  //  __m128d vd;
  uint64 u64[2];
  //  double d8[2];
  //  float f4[4];
  //uint32_t u4[4];
} __uni16;
#endif


#define bits 4L
#define blocklength 16L
#define snpsPerCompressed 16L
#define snpsPer128 32L
#define repetSSE2 10L
#define atonceSSE2 3L
#define atonce_SSE3 4L

SEXP sse_check() {
  SEXP dummy;
  PROTECT(dummy=allocVector(INTSXP, 1));
  INTEGER(dummy)[0] =
#ifdef SSE2
    1
#else
    0
#endif
    ;
  UNPROTECT(1);
  return dummy;
}


SEXP Information = R_NilValue,
  Method = R_NilValue;
void static Init() {
  Information = install("information");
  Method = install("method");
}


SEXP static create_codevector(Uint what, Uint snps, Uint individuals) {
  Uint method = GLOBAL.relationship.method;
  //  if (method == Hamming2) { if (snps % 96) BUG; }
  //  else if (method == Hamming3) { if (snps % 128) BUG; }
  //  else BUG;
  Uint
    f = 2 * (method == Hamming2 ? repetSSE2 * atonceSSE2 : atonce_SSE3), // wechsel 64 bit auf 128; 3 bzw 4
    //                                       auf ein Mal; sse2: 10x wiederholt
    mem = individuals * (1L + (snps - 1L) / (f * snpsPerCompressed)) *f,
    mem2 = 2 * mem; // M, M_T2

  SEXP Code =  CreateEmptyCodeVector(snps, individuals, mem2, REALSXP,
				     Information, true),
    Infos = getAttrib(Code, Information);
  Uint *pi = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  pi[WHAT] = what; // first SNP in SNP sequence
  pi[MEM] = mem; // not mem2!!
  

//if (PRESERVE) {BUG;  R_PreserveObject(Code);}
  return(Code);
}


static uint64 coding1[] = {0, 5, 10},    //0, 1, 2
  coding2[] = {0, 3, 15},    //0, 1, 2
  coding3[] = {0, 6, 15},    //0, 1, 2
  coding4[] = {0, 3, 15};    //0, 1, 2
int rev_coding1[16] = {0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0},
    rev_coding3[16] = {0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 2};
  
void matrix_codingHI(int *SNPxIndiv, int individuals, int snps,
		      uint64 *M, uint64 *M_T, double *p,
		    uint64 *code1, uint64  *code2, Uint n_compressed) {
  
  double *p2 = p + INFO_P_P;
  // coding of M and M^T
  for (int i = 0; i<individuals; i++) {
    int count=0,
      *pSNPxI = SNPxIndiv + i * snps;
    uint64 *Mptr = ((uint64*) M) + i * n_compressed,
      *M_Tptr = ((uint64*) M_T) + i * n_compressed;
    // *Mptr = *M_Tptr = 0; already done
    for (int s = 0; s<snps; s++) {
      int shift = 60-bits*count;
      *Mptr |= code1[pSNPxI[s]] << shift;
      *M_Tptr |= code2[pSNPxI[s]] << shift;
      p2[s] += (double) pSNPxI[s];
      count = (count + 1) % snpsPerCompressed;
      if (count == 0) {
	Mptr++;
	M_Tptr++;
      }
    }
  }
     
  double inv_individuals = 1.0 / (double) individuals;
  double
    ppt2 = 0.0,
    sump = 0.0,
    sum_pi_qi2 = 0.0;
  for (int s=0; s<snps; s++) {
    sump += p2[s];   
    double ps = (p2[s] *= inv_individuals),
      p22 = ps * ps;
    ppt2 += p22;
    sum_pi_qi2 += ps - 0.5 * p22;
  }
  p[INFO_P_PPT] = ppt2;
  p[INFO_P_SUMPQ] = sum_pi_qi2;
  p[INFO_P_SUMP] = sump;
}


SEXP matrix_codingH(SEXP SNPxIndiv){
  if (Information == R_NilValue) Init();
   bool hamming2 = GLOBAL.relationship.method == Hamming2;
  int
    snps = nrows(SNPxIndiv),
    individuals = ncols(SNPxIndiv);
  SEXP Code = create_codevector(GENOMATRIX, snps, individuals),
     Infos = getAttrib(Code, Information);
  Uint *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
     n_compressed = info[MEM] / individuals;
  uint64 
    *M = (uint64*) REAL(Code),
    *M_T = (uint64*) (REAL(Code) + info[MEM]);
  matrix_codingHI(INTEGER(SNPxIndiv), individuals, snps, M, M_T,
		   REAL(VECTOR_ELT(Infos, INFO_P)),
		   hamming2 ? coding1 : coding3,
		 hamming2 ? coding2 : coding4,
		 n_compressed);
  if (PL > 1) {
    Uint mem = info[MEM];
    PRINTF("Data: %d individuals and %d SNPs\nSize of M and M_T each: %d MB\n",
	   individuals, snps,  (int) (sizeof(uint64) * mem / 1048576));
  }
    
  return Code;
}




SEXP matrixH_get(SEXP SNPxIndiv) {
  if (Information == R_NilValue) Init();
  bool hamming2 = INTEGER(getAttrib(SNPxIndiv, Method))[0] == Hamming2;
  int *rev = hamming2 ? rev_coding1 : rev_coding3;
  SEXP
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  Uint 
    individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    n_compressed = info[MEM] / individuals,
    endfor = info[MEM];
  SEXP Ans;
  PROTECT(Ans=allocMatrix(INTSXP, snps, individuals));
  //if (PL > 1) PRINTF("Data: %d individuals and %d SNPs\n", individuals, snps);
  Uint *ans = (Uint *) INTEGER(Ans);
  uint64 *M = (uint64 *) REAL(SNPxIndiv);
  for (Uint a=0; a<endfor; a+=n_compressed) {
    uint64 *Ma = M + a;
    for (Uint s=0; s<snps; s++) {
      *(ans++) = rev[(Ma[s / snpsPerCompressed]
		      >> (60 - 4  * (s % snpsPerCompressed))) & 15L];
    }
  }
  UNPROTECT(1);
  return Ans;
}



#define LOADU _mm_loadu_si128
#define SET8 _mm_set1_epi8
#define SET16 _mm_set1_epi16
#define ADD64 _mm_add_epi64
#define AND _mm_and_si128
#define OR _mm_or_si128
#define SHR64 _mm_srli_epi64
#define ZERO  _mm_setzero_si128()
void matrix_mult2(double *M, double *M_T, int individuals, int snps,
		  double *ergb) {
  // M*M^T:
  // outer for-loop represents rows of M (only within partition, if > 1 cores
  // inner for-loop represents columns of M^T starting at diagonal element
  int
    f = repetSSE2 * atonceSSE2, // 10 * 3
    compr128 = (1L + (snps - 1L) / (f * snpsPer128)) * f,
    loops = compr128 / f;
  const __m128i
    m1 = SET8(0x55),
    m2 = SET8(0x33),
    m4 = SET8(0x0f),
    m8 = SET16(0x00ff);
  uint64 
    OOO1 = INT64_C(0x0001000100010001);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 20)
#endif 
  for (int i = 0; i<individuals; i++) {
    __m128i
      *mpp = ((__m128i*) M) + i * compr128,
      *m_tpp = ((__m128i*) M_T) + i * compr128;
    for (int j = i; j<individuals; j++) {
      __m128i *mptr = mpp; //mptr set back to first element of row in 	
      uint64
	tot = 0;
      for (int s=0; s<loops; s++) {
	__uni16 acc;
	acc.vi = ZERO;
	// do-while runs 10x before continuing
	for (int t=0; t<repetSSE2; t++) {
	  __m128i
	    half1 = AND(LOADU(mptr++), LOADU(m_tpp++)),
	    half2 = AND(SHR64(half1, 1), m1);
	  half1 = AND(half1, m1);
	  half1 = OR(half1, half2);
	  __m128i count1 = AND(LOADU(mptr++), LOADU(m_tpp++));
	  count1 = ADD64(count1, half1);	  
	   __m128i count2 = AND(LOADU(mptr++), LOADU(m_tpp++));
	  count2 = ADD64(count2, half2);
#define dummy1 half1
#define dummy2 half2
	  dummy1 = AND(count1, m2);
	  dummy2 = AND(SHR64(count1, 2), m2);
	  count1 = ADD64(dummy1, dummy2);
	  dummy1 = AND(count2, m2);
	  dummy2 = AND(SHR64(count2, 2), m2);
	  dummy1 = ADD64(dummy1, dummy2);
	  count1 = ADD64(count1, dummy1);
	  dummy1 = AND(count1, m4);
	  dummy2 = AND(SHR64(count1, 4), m4);
	  dummy1 = ADD64(dummy1, dummy2);
	  acc.vi = ADD64(acc.vi, dummy1);
	}
	acc.vi = ADD64(AND(acc.vi, m8),
		      AND(SHR64(acc.vi, 8), m8));
	tot += ((acc.u64[0] + acc.u64[1]) * OOO1) >> 48;
      }
      ergb[j + i * individuals] = ergb[i + j * individuals] = (double) tot;
    }
  }
}
 
#define SHR16 _mm_srli_epi16
#define SHUFFLE8 _mm_shuffle_epi8
#define ADD8 _mm_add_epi8
void matrix_mult3(double *M, double *M_T, int individuals, int snps,
		   double *ergb) {   
  // M*M^T:    
  int
    compr128 = (1L + (snps - 1L) / ( atonce_SSE3 * snpsPer128)) * atonce_SSE3,
    loops = 1 + (compr128 -1) / atonce_SSE3;
  
  const unsigned _LUT[] = {0x02010100, 0x03020201, 0x03020201, 0x04030302};
  const __m128i LUT = _mm_load_si128((__m128i*) _LUT), // _mm_set_epi32(,,,)
    mask = SET8(0x0F),
    zero = _mm_setzero_si128();
  
  __uni16  intermed;
  intermed.vi =_mm_setzero_si128();
   
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 20)
#endif
  for (int i = 0; i<individuals; i++) {
    __m128i
      *mpp = (__m128i*) M + i * compr128,
      *m_tpp = (__m128i*) M_T + i * compr128;
    for (int j = i; j<individuals; j++) {
      int tot = 0;
      __m128i *mptr = mpp;     
      for (int k = 0; k < loops; k++) {    // '/4' hinzugefuegt
	__m128i v0, lo, count0, count1;
#define COUNT(count0)						\
	v0 = AND(LOADU(mptr++), LOADU(m_tpp++));		\
	lo = AND(mask, v0);					\
	v0 = AND(mask, SHR16(v0,4)); /* hi */			\
	lo = SHUFFLE8(LUT, lo);					\
	v0 = SHUFFLE8(LUT, v0);					\
	count0 = ADD8(lo, v0);
	
	COUNT(count0);
	COUNT(v0);
	count0 = ADD8(count0, v0);
	COUNT(count1);
	COUNT(v0);
	count1 = ADD8(count1, v0);
	count0 = ADD8(count0, count1);
	
	//bitwise difference between 'tota' and all zeros
	intermed.vi = _mm_sad_epu8(count0, zero);
	// adds the two bitcounts which are saved as unsigned 64bit integers
	// to continuous bitcount
	tot += (int) intermed.u64[0] + (int)intermed.u64[1];
      }
      ergb[j + i * individuals] = ergb[i + j * individuals] = (double) tot;
    }
  }
}


SEXP matrix_mult(SEXP M) {
  SEXP Ans,
    Infos = getAttrib(M, Information);
  Uint *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  int
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    total = info[MEM];
  PROTECT(Ans = allocMatrix(REALSXP, individuals, individuals));
  double *m = REAL(M),
    *ans = REAL(Ans);
    
  int meth = INTEGER(getAttrib(M, Method))[0];
  if (meth == Hamming2) matrix_mult2(m, m + total, individuals, snps, ans);
  else matrix_mult3(m, m + total, individuals, snps, ans);
 
  //  double
  //    *PP = REAL(VECTOR_ELT(Infos, INFO_P)),
  //    inv_sum_pi_qi = 1.0 / PP[INFO_P_SUMPQ];

  DoCentering(Ans, individuals, 
	      REAL(VECTOR_ELT(Infos, INFO_P))[INFO_P_SUMP]);
  
  
  UNPROTECT(1);
  return Ans;
}




