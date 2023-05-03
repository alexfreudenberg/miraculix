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


#define BitsPerCode 2
#define MY_CODING Plink
#define MY_LDABITALIGN MY_LDABITALIGN_PLINK

#define MY_VARIANT 32
#define BytesPerPartUnit 2


#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "TemplateUint.h"
#include "bitUint.h"
#include "plink.h"

#include "Haplo.h"
#include "MX.h"


#define RIGHT ((unit_t) 0x5555555555555555)
#define LEFT ((unit_t) 0xAAAAAAAAAAAAAAAA)
#define PLINKTRAFO							\
  , rPlink = twobit & RIGHT, lPlink = twobit & LEFT,			\
    mPlink=lPlink & rPlink << 1;			     \
  twobit = mPlink | ((mPlink xor lPlink) >> 1)

Uint SUM_MISSINGS[256];


Long LdaPlink(Long snps, Long LDAbitalign) {
  return LDAalignBitSnps(snps * BitsPerCode); 
}

Long LdaOrigPlink(Long snps, Long LDAbitalign) {
  if (LDAbitalign != BitsPerByte && LDAbitalign != 0) BUG;
  return ROUND_GEQ(snps * BitsPerCode, BitsPerByte);
}


void InitPlink() {
  assert(sizeof(unit_t) <= 8); // def RIGHT, LEFT
  for (Uint i3=0; i3<4; i3++) {
    Uint V3 = i3;
    Uint missings3 = (Uint) (i3 == 1U);
    for (Uint i2=0; i2<4; i2++) {
      Uint V2 = 4U * V3 + i2;
      Uint missings2 = missings3 + (Uint) (i2 == 1U);
      for (Uint i1=0; i1<4; i1++) {
	Uint V1 = 4U * V2 + i1;
	Uint missings1 = missings2 + (Uint) (i1 == 1U);
	for (Uint i0=0; i0<4; i0++) {
	  Uint V0 = 4U * V1 + i0;
	  Uint missings0 = missings1 + (Uint) (i0 == 1U);
	  SUM_MISSINGS[V0] = missings0;
	}
      }
    }
  }
}


get_matrix_start(_4Byte, Plink)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix_end(_4Byte,PLINKTRAFO)  

get_matrix_start(_1Byte, Plink)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix_end(_1Byte,PLINKTRAFO)   


coding_start(_4Byte,Plink)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
coding_end(_4Byte, HUMAN2PLINK(*mm) )
 
coding_start(_1Byte,Plink)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
coding_end(_1Byte, HUMAN2PLINK(*mm) )
 
 


coding2trans_start(_4Byte, Plink)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif								
coding2trans_do(_4Byte, HUMAN2PLINK)

coding2trans_start(_1Byte, Plink)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding2trans_do(_1Byte, HUMAN2PLINK)



vectorGeno_start(Ulong, Ulong, Ulong, Plink);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
vectorGeno_fctn(Ulong, Ulong, Ulong, PLINKTRAFO)


vectorGeno_start(double, double, double, Plink);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
vectorGeno_fctn(double, double, double, PLINKTRAFO)


vectorGeno_start(LongDouble, LongDouble, LongDouble, Plink);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
vectorGeno_fctn(LongDouble, LongDouble, LongDouble, PLINKTRAFO)


sumGeno_start(Plink) 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
sumGeno_fctn(PLINK2HUMAN)



void getPlinkMissings(Uchar* plink, SEXP SxI, basic_options *opt) {
  int cores = GreaterZero(opt->cores);
  Long null, eins, 
     *info = GetInfo(SxI);
  Long columns, ldaByte,
    snps = info[SNPS],
    indiv = info[INDIVIDUALS];    
  Long *P = NULL;
  Long length; // in Bytes
  if (plink == NULL) {
    ldaByte = info[LDA] * BytesPerUnit;
    length = info[LDA] * info[INDIVIDUALS];
    P = (Long*) Align(SxI, opt);
    columns = indiv;
    null = 0;
    eins = 1;
  } else {
    ldaByte = DIV_GEQ(indiv, CodesPerByte);
    length = ldaByte * snps;
    P = (Long*) plink;
    columns = snps;
    null = 1;
    eins = 0;    
  }
  Long missings = 0;
  Long lenLong = length / sizeof(Long);
  Long i = 0;
  
  for(; i < lenLong; i+=sizeof(Long), P++) {
    Long anymiss = ( (P[i] xor 0x5555555555555555) & 0x5555555555555555 &
		      (P[i] >> 1) );
    if (anymiss) {
      Uchar *p = (Uchar*) (P + i);
      for (Ulong j=0; j<sizeof(Long); j++) missings += SUM_MISSINGS[p[j]];
    }
  }
  Uchar*p = (Uchar*) P;
  for(; i < length; i++, p++) missings += SUM_MISSINGS[p[i]];
  
  SEXP Miss = R_NilValue;
  info[MISSINGS] = missings;
  if (missings) {
    Miss = PROTECT(allocMatrixExtra(LONGSXP, 2, missings));
    setAttrib(SxI, Missings, Miss);
    Long *pM = LONG(Miss);

    
    Long mc = 0;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
    for (Long col=0; col < columns; col++) {
      Uchar *code = plink + ldaByte * col;
      for (Long row=0; row<ldaByte; row++, code++) {
	Long m = SUM_MISSINGS[code[row]]; // assuming only few missings ...
	if (m) {
	  for(Uint k=0; k<CodesPerByte; k++) {
	    if ((m & 0x03) == 0x01) {
	      pM[mc + null] = col;
	      pM[mc + eins] = row * CodesPerByte + k;
	      mc += 2;
	    }
	    m >>= BitsPerCode;
	  }
	} // m
      } // row
    } // col
    
    UNPROTECT(1);
  } // missings
}

// up to now miniIndivChunk 4 staticSize 32, mini 4


#define staticSize 32
#define miniSnpChunkFactor 4

/*
// miniSnpChunkFactor 4, row = -4
1 : 2.830 ( 26.899)
2 : 2.173 ( 20.304) 2.245 ( 20.924) 
4: 1.573 ( 14.349) 
8:1.115 (  9.757) 1.095 (  9.624) 1.138 (  9.957)
10:1.125 (  9.903)1.121 (  9.823) 1.182 ( 10.497)
16:0.953 (  8.185)0.992 (  8.488) 0.951 (  8.197) 
   0.977 (  8.337)0.955 (  8.188) 0.951 (  8.181)
   0.964 (  8.190)0.964 (  8.178) 0.964 (  8.214) 
32:0.888 (  7.484) 0.888 (  7.484) 
64: 0.922 (  7.862) 0.915 (  7.807)0.920 (  7.810) 
128: 0.917 (  7.774)
256: 0.907 (  7.65
512:1.782 ( 16.428) 
1024:2.193 ( 20.435)
 */

/*
 static 10: row =-4 
 1:1.200 ( 10.626) 1.158 ( 10.136) 1.196 ( 10.606)
 2: 1.063 (  9.152) 1.126 (  9.933)1.080 (  9.440)
 3: 1.065 (  9.321) 1.072 (  9.362)1.088 (  9.426
 4:  1.073 (  9.319)  1.058 (  9.250)1.053 (  9.065) 
 7:1.344 ( 11.291)  1.128 (  9.925)1.211 ( 10.645) 
 10:1.125 (  9.903)1.121 (  9.823) 1.182 ( 10.497)
 40:1.184 ( 10.375)1.297 ( 10.770) 1.195 ( 10.541)
 400: 1.548 ( 14.084)

 static 32: wait
 2: 0.902 (  7.619) 0.894 (  7.614)
 4:0.888 (  7.484) 0.888 (  7.484) 
 8: 0.958 (  8.172)
40:1.114 (  9.803)

*/



// row = -120
//#define miniIndivChunk 4
//#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3)
// 13.607 (134.684) 14.653 (145.047) 14.099 (139.495) 

//#define miniIndivChunk 8
//#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7)
//   11.578 (114.293)  12.328 (121.752)  13.067 (129.134)  12.592 (124.381)

//#define miniIndivChunk 12
//#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7) SZ X(8) SZ X(9) SZ X(10) SZ X(11)
//   9.854 ( 97.135) 11.074 (109.180)  11.075 (109.186) 11.128 (109.955) 


// row = -240
// #define miniIndivChunk 8
// 22.930 (227.752) 23.786 (236.256) 22.819 (226.631) 22.341 (221.874) 
// #define miniIndivChunk 12, 
// 19.811 (196.704) 24.043 (238.450) 22.760 (225.573)  22.548 (223.428) 20.150 (199.923)
// #define miniIndivChunk 16, 
// 21.602 (214.604)  20.087 (199.287) 21.612 (214.435) 20.586 (204.267)


//#define miniIndivChunk 24
//#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7) SZ X(8) SZ X(9) SZ X(10) SZ X(11) SZ X(12) SZ X(13) SZ X(14) SZ X(15) SZ X(16) SZ X(17) SZ X(18) SZ X(19) SZ X(20) SZ X(21) SZ X(22) SZ X(23)
// 9.524 ( 93.788)  11.823 (116.668)  11.172 (110.222) staticSize 32  mini 4
// 10.682 (105.412)  11.146 (109.841)staticSize 32  mini 8
//  10.891 (107.580)11.532 (113.978) 11.932 (117.855) staticSize 32  mini 2


#define POINTER_TO_CODE(N) \
  Ulong *snp##N = (Ulong*) (c + colIdxB[iChunk + N] * ldaInByte)
#define SNP_POINTER_TO_CODE MULTI(POINTER_TO_CODE,;)

#define V_GET_FROM_B(N) const double v##N = valueB[iChunk + N]
#define GET_V_FROM_B MULTI(V_GET_FROM_B,;)

// the next could be improved by separate code or at least, apply to the last only (assuming that colIdxB is rowwise ordered
#define TMP_LEN 8
#define SNPS_GET_FROM_POINTER(N) Ulong snps##N;				\
  if (snp##N + f + 1 <= snpsEnd) snps##N = snp##N[f];		\
  else {					\
    Uchar *S = (Uchar*) (snp##N + f);					\
    Uchar TMP[TMP_LEN] = {						\
      S + 0 < (Uchar*) snpsEnd ? S[0] : (Uchar) 0,			\
      S + 1 < (Uchar*) snpsEnd ? S[1] : (Uchar) 0,			\
      S + 2 < (Uchar*) snpsEnd ? S[2] : (Uchar) 0,			\
      S + 3 < (Uchar*) snpsEnd ? S[3] : (Uchar) 0,			\
      S + 4 < (Uchar*) snpsEnd ? S[4] : (Uchar) 0,			\
      S + 5 < (Uchar*) snpsEnd ? S[5] : (Uchar) 0,			\
      S + 6 < (Uchar*) snpsEnd ? S[6] : (Uchar) 0,			\
      S + 7 < (Uchar*) snpsEnd ? S[7] : (Uchar) 0			\
    };									\
    snps##N = ((Ulong*) TMP)[0];					\
  }								         

#define GET_SNPS_FROM_POINTER MULTI(SNPS_GET_FROM_POINTER,;)

#define PROD(N) ((snps##N & CodeMask) <= 1 ? 0.0		\
		 : (snps##N & CodeMask) == CODE_ONE ? v##N : (2.0 * v##N))
#define SUM_UP MULTI(PROD,+)

#define SNPS_SHR2(N) snps##N >>= BitsPerCode
#define SHR2_SNPS MULTI(SNPS_SHR2,;)

#define CALCULATE					\
  SNP_POINTER_TO_CODE;					\
  GET_V_FROM_B;						\
  for (int f=0; f<genuineMiniFactor; f++) {		\
    double *tmp0 = tmp + f * CodesPerLong;		\
    GET_SNPS_FROM_POINTER;					\
    for (int iCompr=0; iCompr < CodesPerLong; iCompr++) {	\
      tmp0[iCompr] += SUM_UP;					\
      SHR2_SNPS;						\
    }								\
  }								\
  iChunk+=miniIndivChunk;					\
  nonzeros-=miniIndivChunk

sparseTGeno_start(OrigPlink)
  // NOTE: for better readability only, redefine nrow & ncol by their 
  //       typical meaning in (sparse) %*% t(code), namely snps and individuals
  const Long snps = nrow_code;       
  const int max_switch = 8; // 16 has 5%-10% gain wrt 8 for "many" nonzeros
  
  const int chunkNsnp = CodesPerLong * miniSnpChunkFactor;
  const Ulong *snpsEnd = (Ulong*) (code + ncol_code * ldaInByte);
  for (Long iRchunk=0; iRchunk < snps; iRchunk += chunkNsnp) {
    //    printf("%d ", iRchunk);
    const Long byteIdx_snp = iRchunk / CodesPerByte;
    const int lenRchunk = chunkNsnp <= snps - iRchunk ? chunkNsnp
      : (int) (snps - iRchunk);
    const int genuineMiniFactor = DIV_GEQ(lenRchunk, CodesPerLong);
    const Uchar *c = code + byteIdx_snp;
    double *ans = Ans + iRchunk * ldAns;
    assert(sizeof(Ulong) == TMP_LEN);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static,staticSize) ordered 
#endif 
    for (int j=0; j<nIdx; j++) {
      int nonzeros = rowIdxB[j+1] - rowIdxB[j]; // CSR/CSR/Yale
      if (nonzeros == 0) continue;
      double *a = ans + j;     
      double tmp[chunkNsnp] = {0};
      int iChunk=rowIdxB[j];      
      while (nonzeros > 0) {
	switch(nonzeros > max_switch ? max_switch : nonzeros) {
	case 1 : {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI
#endif
#define miniIndivChunk 1
#define MULTI(X,SZ) X(0)
	  CALCULATE;
	}
	  break;
	case 2 : {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI
#endif
#define miniIndivChunk 2
#define MULTI(X,SZ) X(0) SZ X(1)
	  CALCULATE;
	}
	  break;
	case 3 : {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI
#endif
#define miniIndivChunk 3
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2)
	  CALCULATE;
	}
	  break;
	case 4 : case 7 : {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI	
#endif
#define miniIndivChunk 4
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3)
	  CALCULATE;
	}
	  break;
	case 5 :  { // 5% gain
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI	
#endif
#define miniIndivChunk 5
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4)
	  CALCULATE;
	}
	  break;
	case 6 :  { // 5% gain
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI	
#endif
#define miniIndivChunk 6
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5)
	  CALCULATE;
	}
	  break;
	case 8 : case 9 : case 10 : case 11 :    {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI	
#endif
#define miniIndivChunk 8
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7)
	  CALCULATE;
	}
	  break;
	case 12 : case 13 : case 14 : case 15 :  {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI	
#endif
#define miniIndivChunk 12
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7) SZ X(8) SZ X(9) SZ X(10) SZ X(11)
	  CALCULATE;
	}
	  break;
	
	case 16 : {
#if defined MULTI	  
#undef miniIndivChunk
#undef MULTI	
#endif
#define miniIndivChunk 16
#define MULTI(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7) SZ X(8) SZ X(9) SZ X(10) SZ X(11) SZ X(12) SZ X(13) SZ X(14) SZ X(15)
	  CALCULATE;
	}
	break;
	default : BUG;
	} // switch	
      } // while
      for (int k=0; k < lenRchunk; k++) a[k * ldAns] = tmp[k];
    }// j
  }// iRchunk
}

