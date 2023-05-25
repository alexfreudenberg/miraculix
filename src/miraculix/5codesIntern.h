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


#ifndef miraculix_5codesIntern_H
#define miraculix_5codesIntern_H 1


#ifdef DO_PARALLEL
#include "omp.h"
#endif


#define AtOnce 4
#define SumAtOnce 4
#define SLICELEN 100
#define RoughRowChunk 35000

/* 

   RoughRowChunk:

Wageningen/run_gcc snps=150000 indiv=150000 repetV=32 cores=10 floatLoop=0 meanSubst=0 missingsFully0=1 centered=0 trials=10 variant=256 coding=5 SxInotVmean=0 mode=3 SNPmode=0 cmp=0 sparse=0 rowMode=-4 colMode=0 fillMode=0

150000   1 : 7.5
 75000   2 : 5.4
 50000   3 : 5.2
 37500   4 : 5.1
 30000   5 : 5.2
 20000   7 : 5.2
 15000  10 : 5.4
 10000  15 : 5.6

*/

//#define AlteVersion 1
#define AlteVersion 0

#define coreFactor 5 // 0.25 // oder 0.2
  
#define gV5_start0(NR, TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)	\
  int cores = GreaterZero(opt->cores); /* NOT Long !! */      		\
  if ((opt->Cprintlevel && sizeof(TRDTYPE) == 8) || opt->Cprintlevel > 3) { \
    char msg[200];							\
    SPRINTF(msg, "%s Bit %s-%s required for %ld x %ld x %ld", 		\
	    #NR, #TRDTYPE, #ENDTYPE, rows, cols, repetV);		\
    PRINTF("%s", msg); for (int i=(int) STRLEN(msg); i<57; i++) PRINTF(" "); \
    PRINTF(" -> %3ld Bit, %2d cores\n",					\
	   sizeof(AVXTYPE) * BitsPerByte, cores);			\
  }									\
  MEMSET(Ans, 0, ldAns * repetV * sizeof(ENDTYPE));			\
  TRDTYPE *externalV = NULL;						\
  const Long colsCpB = DIV_GEQ(cols, CpB);				\
  const Long colsCpBm1 = colsCpB - 1;					\
  const Long colBlocks = DIV_GEQ(colsCpB,  AtOnce);			\
  Long blockSliceLen = DIV_GEQ(colBlocks, cores * coreFactor);		\
  if (false) printf("blockSliceLen=%ld\n", blockSliceLen);		\
  blockSliceLen = MAX(1, MIN(SLICELEN, blockSliceLen));			\
  Long blocks = DIV_GEQ(colBlocks, blockSliceLen); /* (ca. cores * 5) repetV */; \
  /* for large repetV, use repetV only for splitting */			\
  /* blocks = blocks > cores ? cores : blocks > 1 ? blocks : 1;	*/	\
  const Long blocksM1 = blocks - 1L;					\
  const Long blockRest = colBlocks - blockSliceLen * blocksM1;		\
  const Long sliceLen = blockSliceLen * AtOnce;				\
  const Long rest = blockRest * AtOnce;					\
  if (false) PRINTF("0 < rest =%ld blcokrest=%ld<=%ld ------%ld   cores=%d atOnce=%d blocks=%ld\n", rest, blockRest, blockSliceLen, colBlocks, cores, AtOnce, blocks); \
  assert(rest > 0 && blockRest <= blockSliceLen);			\
  const Long blocksXrows = blocks * rows;				\
  const Long HashSize = 243L;						\
  const Long Hash_BytesPcol = HashSize * sizeof(AVXTYPE);		\
  const Long HashSizePerV = (colBlocks * AtOnce) * HashSize;		\
  const Long NVblocks = DIV_GEQ(repetV, VatOnce);			\
  AVXTYPE *F = (AVXTYPE*) CALLOC(HashSizePerV * NVblocks, sizeof(AVXTYPE)); \
  const Long blocksP1Xrows = blocksXrows + rows;			\
  const Long TmpSize = blocksP1Xrows * NVblocks				\
    +  AlteVersion * (blocksP1Xrows * NVblocks); \
  AVXTYPE *Tmp = (AVXTYPE*) CALLOC(TmpSize, sizeof(AVXTYPE));		\
  AVXTYPE *sum = Tmp + AlteVersion * (blocksP1Xrows * NVblocks); \
  const Long VARIABLE_IS_NOT_USED nr_trds_in_avx =			\
    BytesPerBlock / sizeof(TRDTYPE);					\
  if (false) {PRINTF("rows=%ld, cols=%ld\n", rows, cols);}		\
  if (!true) if (repetV % VatOnce) ERR0("'VatOnce not a multiple of repetV"); \
  if (!true) assert(sizeof(AVXTYPE) % sizeof(TRDTYPE) == 0);		\
    
  //  printf("blocksM1=%ld, rr=%ld repet = %ld cores=%d; colsCpBm1=%ld %s\n", blocksM1,rr, repetV, cores, colsCpBm1, #NR);



#if defined AVX2 || defined AVX512F
#define gV5_start(NR, TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)	\
  gV5_header(NR, TYPE, SNDTYPE, TRDTYPE, ENDTYPE) {		\
  const Long VatOnce = (sizeof(AVXTYPE) / sizeof(TRDTYPE));	\
  const Long restV = repetV % VatOnce;				\
  if (restV > 0 &&						\
      ((VatOnce == 8 && restV <= 4) || (VatOnce == 4 && restV <= 2))) {	\
    vector_##TYPE##_t vD = restV == 1 ? genoVector5v32_##TYPE		\
      : restV == 2 ? genoVector5v128_##TYPE				\
      : genoVector5v256_##TYPE;						\
    vD(code, rows, cols, lda, coding, variant, pseudoFreq,		\
       V, restV, ldV, opt, tuning, Ans, ldAns);				\
    repetV -= restV;							\
    V += restV * ldV;							\
    Ans += restV * ldAns;						\
  }									\
  gV5_start0(NR, TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)
#else
#define gV5_start(NR, TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)	\
  gV5_header(NR, TYPE, SNDTYPE, TRDTYPE, ENDTYPE) {		\
  const Long VatOnce = (sizeof(AVXTYPE) / sizeof(TRDTYPE));	\
  gV5_start0(NR, TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)
#endif


#define gV5_CreateHash(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)	\
  for (Long rr=0; rr < NVblocks; rr++) {				\
    AVXTYPE *f = F + rr * HashSizePerV;					\
    Long sEnd = MIN(VatOnce, repetV - rr * VatOnce);			\
    for (Long i=0; i<=colsCpBm1; i++) {					\
      TYPE *v_i = V + rr * VatOnce * ldV + i * CpB;			\
      AVXTYPE *hash0 = f + i * HashSize;				\
      int hastRest = CpB;						\
      SNDTYPE x[CpB] = {0, 0, 0, 0, 0};/* hier compiler fehler bei collapse!?! */ \
      if (i < colsCpBm1) {} else {					\
	hastRest = (int) (cols - colsCpBm1 * hastRest);			\
	assert(hastRest > 0);						\
	for (int k=hastRest; k<CpB; k++) x[k] = 0;			\
      }									\
      SNDTYPE mean = 0;							\
      for (Long s=0; s<sEnd; s++) {					\
	TYPE *v = v_i + ldV * s;					\
	TRDTYPE *hash = ((TRDTYPE*) hash0) + s;				\
	for (int k=0; k<hastRest; k++) {x[k] = (SNDTYPE) v[k]; }	\
	if (pseudoFreq != NULL) {					\
	  SNDTYPE meanLD = 0.0;						\
	  for (int k=0; k<hastRest; k++)				\
	    meanLD += (SNDTYPE) pseudoFreq[i * CpB + k] * (2 * x[k]); \
	  mean = meanLD;						\
	}								\
	if (externalV != NULL) {					\
	} else {							\
	  for (Uint i4=0; i4<3; i4++) {					\
	    Uint V4 = i4;						\
	    SNDTYPE fct4 = (SNDTYPE) i4 * x[4] - mean;			\
	    for (Uint i3=0; i3<3; i3++) {				\
	      Uint V3 = 3 * V4 + i3;					\
	      SNDTYPE fct3 = fct4 + (SNDTYPE) i3 * x[3];	\
	      for (Uint i2=0; i2<3; i2++) {				\
                Uint V2 = 3 * V3 + i2;					\
		SNDTYPE fct2 = fct3 + (SNDTYPE) i2 * x[2];	\
	   	for (Uint i1=0; i1<3; i1++) {				\
		  Uint V1 = 3 * V2 + i1;      				\
		  SNDTYPE fct1 = fct2 + (SNDTYPE) i1 * x[1];	\
		  							\
		  Uint V0 = 3 * V1 * VatOnce;				\
		  SNDTYPE fct0 = fct1;			\
		  hash[V0] = (TRDTYPE) fct0;			\
		  V0+=VatOnce;						\
		  fct0 += x[0];						\
		  hash[V0] = (TRDTYPE) fct0;			\
		  hash[V0 + VatOnce] = (TRDTYPE) (fct0 + x[0]);	\
		}							\
	      }								\
	    }								\
	  }								\
	}								\
      }									\
    }									\
  }									\
  if (false) PRINTF("Hash end\n");					\
  const Long rowBlocks = MAX(1, DIV_LEQ(rows, RoughRowChunk));		\
  const Long RowChunk =  DIV_GEQ(rows, rowBlocks);			\
  if (false) PRINTF("RowChunk = %ld, rows=%ld, %f\n", RowChunk, rows, rows / (double) RowChunk); \
  for (Long bStart=0; bStart<rows; bStart+=RowChunk) {			\
    Long bEnd = MIN(bStart + RowChunk, rows);				\

/* 
Aufwand fuer Erstellung der Hashtabelle kann deutlich reduziert wreden,
wenn SNDTYPE==TRDTYPE
(i)  i = 0
     AVX-Vektor X0 mit Werten x[0] x VatOnce; 
(ii) 3 AVX-Vektoren: (a) 0  (b) Xi (c) 2 * Xi
(iii) i++
(iv) AVX-Vektor Xi mit Werten x[i] x VatOnce erstellen
(v)
     (a) alle Vektoren kopieren auf die direkt nachfolgenden Speicherplaetze
     (b) das Kopierte plus Xi
     (c) Ergebnis kopieren auf die direkt nachfolgenden Speicherplaetze
     (c) das Kopierte plus Xi
(vi) weiter mit (iii) bis i=4 erreicht

Von V lokal eine Kopie halten? 
V umsortieren, so dass in Reihenfolge ausgelesen wird?
*/
  

#define ccc 0

#define gV5_LoopA(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
  for (Long C=0; C <= blocksM1; C++) {					\
    for (Long rr=0; rr < NVblocks; rr++) { /* all repetV, as size is small*/ \
      AVXTYPE *ff = F + rr * HashSizePerV + sliceLen * C * HashSize;	\
      AVXTYPE *t = Tmp + rr * blocksP1Xrows + rows * ( C);		\
      unit_t *c = code + lda * sliceLen * C;				\
      Long nrCols = C == blocksM1 ? rest : sliceLen;			\
      /* assert(nrCols > 0); */						\
      for (Long i = 0 ; i<nrCols; i+=AtOnce) {/* cols of compressed */	\
        Uchar *pC0 = (Uchar*) (c + i * lda);				\
	AVXTYPE ff0[HashSize];						\
	if (true) MEMCOPY(ff0, ff + (i + 0L) * HashSize, Hash_BytesPcol); \
	else MEMSET(ff0, ccc, Hash_BytesPcol);				\
	Uchar *pC1 = (Uchar*) (c + (i+1L) * lda);			\
	AVXTYPE ff1[HashSize];						\
	if (true) MEMCOPY(ff1, ff + (i+1L) * HashSize, Hash_BytesPcol); \
	else MEMSET(ff1, ccc, Hash_BytesPcol);				\
	Uchar *pC2 = (Uchar*) (c + (i+2L) * lda);			\
	AVXTYPE ff2[HashSize];						\
	if (true) MEMCOPY(ff2, ff + (i+2L) * HashSize, Hash_BytesPcol); \
	  else MEMSET(ff2, ccc, Hash_BytesPcol);			\
	Uchar *pC3 = (Uchar*) (c + (i+3L) * lda);			\
	AVXTYPE ff3[HashSize];						\
	if (true) MEMCOPY(ff3, ff + (i+3L) * HashSize, Hash_BytesPcol); \
	else MEMSET(ff3, ccc, Hash_BytesPcol);				\
	for (Long b=bStart; b<bEnd; b++) { /* rows of compressed */	\
  
/*
  double bbb = i;							
  if (false) for (Long jjj=0; jjj<1000; jjj++) bbb += EXP((double) i * jjj); 
  char ccc = bbb > 0;							
*/


#if AtOnce == 4
#if defined AVX2
#define gV5_LoopB(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
    const AVXTYPE fC0 = AVXLOAD((TRDTYPE*)(ff0 + (int) pC0[b]));	\
    const AVXTYPE fC1 = AVXLOAD((TRDTYPE*)(ff1 + (int) pC1[b]));	\
    const AVXTYPE fC2 = AVXLOAD((TRDTYPE*)(ff2 + (int) pC2[b]));	\
    const AVXTYPE fC3 = AVXLOAD((TRDTYPE*)(ff3 + (int) pC3[b]));	\
    AVXTYPE tb = AVXLOAD((TRDTYPE*)(t + b));				\
    const AVXTYPE g0 = AVXADD(fC0, fC1);				\
    const AVXTYPE g1 = AVXADD(fC2, fC3);				\
    const AVXTYPE h  = AVXADD( g0,  g1);				\
    AVXSTORE((TRDTYPE*)(t + b),  AVXADD( tb,  h));	\
  }
#else
#define gV5_LoopB(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)	\
    t[b] +=  (ff0[(int) pC0[b]] + ff1[(int) pC1[b]])		\
      + (ff2[(int) pC2[b]] + ff3[(int) pC3[b]]);		\
  }
#endif

#elif AtOnce == 3
#if defined AVX2
#define gV5_LoopB(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
   const AVXTYPE fC0 = AVXLOAD((TRDTYPE*)(ff0 + (int) pC0[b]));	\
    const AVXTYPE fC1 = AVXLOAD((TRDTYPE*)(ff1 + (int) pC1[b]));	\
    const AVXTYPE fC2 = AVXLOAD((TRDTYPE*)(ff2 + (int) pC2[b]));	\
    AVXTYPE tb = AVXLOAD((TRDTYPE*)(t + b));				\
    const AVXTYPE g0 = AVXADD(fC0, fC1);				\
    const AVXTYPE h  = AVXADD(g0,  fC2);				\
    tb  = AVXADD( tb,  h);						\
    AVXSTORE((TRDTYPE*)(t + b), tb);					\
  }
#else
#define gV5_LoopB(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
    t[b] += ff0[(int) pC0[b]] + ff1[(int) pC1[b]] + ff2[(int) pC2[b]];	\
  }
#endif

#else
#error AtOnce unsound
#endif // AtOnce


#if SumAtOnce == 4
#if defined AVX2
#define Sum4ColsUp(AVXTYPE, TRDTYPE);					\
  const AVXTYPE fC0 = AVXLOAD((TRDTYPE*)(t0 + j));			\
  const AVXTYPE fC1 = AVXLOAD((TRDTYPE*)(t1 + j));			\
  const AVXTYPE fC2 = AVXLOAD((TRDTYPE*)(t2 + j));			\
  const AVXTYPE fC3 = AVXLOAD((TRDTYPE*)(t3 + j));			\
  const AVXTYPE g0 = AVXADD(fC0, fC1);					\
  const AVXTYPE g1 = AVXADD(fC2, fC3);					\
  const AVXTYPE h  = AVXADD( g0,  g1);					\
  AVXSTORE((TRDTYPE*)(t0 + j), h)
#else
#define Sum4ColsUp(AVXTYPE, TRDTYPE);	\
  t0[j] = (t0[j] + t1[j]) + (t2[j] + t3[j])
#endif
#else
#error SumAtOnce unsound
#endif // AtOnce
    
    
#define gV5_MainLoop(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
  gV5_LoopA(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)			\
    gV5_LoopB(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
    }}									\
  }}									\
   if (false) PRINTF("calc end blocks=%ld\n", blocks);			\
    

#if AlteVersion == 0

#define gV5_SumUp(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
  for (Long rr=0; rr < NVblocks; rr++) {				\
    AVXTYPE *tmp = Tmp + rr * blocksP1Xrows;				\
    Long level = rows;							\
    int tmpCols = blocks;						\
    while (tmpCols > 1) {						\
      Long tmpC4 = DIV_GEQ(tmpCols - 1, SumAtOnce);			\
      for (Long k=0; k < tmpC4; k++) {					\
	const Long kS = k * SumAtOnce;					\
	AVXTYPE *t0 = tmp + MIN(blocksXrows, (kS + 0) * level);		\
	AVXTYPE *t1 = tmp + MIN(blocksXrows, (kS + 1) * level);		\
	AVXTYPE *t2 = tmp + MIN(blocksXrows, (kS + 2) * level);		\
	AVXTYPE *t3 = tmp + MIN(blocksXrows, (kS + 3) * level);		\
	for (Long j=0; j<rows; j++)  {					\
	  Sum4ColsUp(AVXTYPE, TRDTYPE);					\
	}								\
      } /* for k */							\
      level *= SumAtOnce;						\
      tmpCols = DIV_GEQ(tmpCols, SumAtOnce);				\
    } /* while */							\
  } /* for */								\
if (false) PRINTF("end sum up\n");				       

#else

#define gV5_SumUp(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
  for (Long rr=0; rr < NVblocks; rr++) {				\
    Long rowsVatonce = rows * VatOnce;					\
    TRDTYPE *s = (TRDTYPE*) (sum + rr * blocksP1Xrows);			\
    for (Long k =0; k < blocks; k+=SumAtOnce) {				\
      AVXTYPE *T = (Tmp + rr * blocksP1Xrows);				\
      TRDTYPE *t0 = (TRDTYPE*) (T + MIN(blocks, k + 0) * rows);		\
      TRDTYPE *t1 = (TRDTYPE*) (T + MIN(blocks, k + 1) * rows);		\
      TRDTYPE *t2 = (TRDTYPE*) (T + MIN(blocks, k + 2) * rows);		\
      TRDTYPE *t3 = (TRDTYPE*) (T + MIN(blocks, k + 3) * rows);		\
      assert((Uchar*) (t3 + rowsVatonce) <= (Uchar*) (T + blocksP1Xrows)); \
      for (Long j=0; j<rowsVatonce; j++) {				\
	s[j] += (TRDTYPE) ((ENDTYPE) t0[j] + (ENDTYPE) t1[j] +		\
			   (ENDTYPE) t2[j] + (ENDTYPE) t3[j]);		\
      }									\
    }									\
  }									\
  if (false) PRINTF("!!! Old version of summing up (%ld, %ld)!!!\n", sizeof(TRDTYPE), sizeof(TRDTYPE));	        
   BUG; // no bug at all, just do not use this SumUp

    
#endif

    

#define gV5_Sort(TYPE, SNDTYPE, TRDTYPE, AVXTYPE, ENDTYPE)		\
  for (Long rr=0; rr < NVblocks; rr++) {				\
    AVXTYPE *tmp = sum + rr * blocksP1Xrows;				\
    Long sEnd = MIN(VatOnce, repetV-rr * VatOnce);			\
    for (Long b=0; b<rows; b++) {					\
      const TRDTYPE *ST= ((TRDTYPE *) tmp) + b * VatOnce;		\
      ENDTYPE *a = Ans + b;						\
      for (Long s=0; s<sEnd; s++) {					\
	a[ldAns * (rr * VatOnce + s)] = (ENDTYPE) ST[s];		\
      }									\
    }									\
  } /* V */								\
  if (false) PRINTF("final end\n");					\
  FREE(F);								\
  FREE(Tmp);								\
  } // end function
 



#if defined AVXLOAD
#undef AVXLOAD
#undef AVXSTORE
#undef AVXADD
#endif
#define AVXLOAD LOADuFLOAT
#define AVXSTORE STOREuFLOAT
#define AVXADD ADDFLOAT
gV5_start(MY_VARIANT, floatD, LongDouble, float, Floats, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_CreateHash(floatD, LongDouble, float, Floats, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
gV5_MainLoop(floatD, LongDouble, float, Floats, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_SumUp(floatD, LongDouble, float, Floats, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_Sort(floatD, LongDouble, float, Floats, double)



#if defined AVXLOAD
#undef AVXLOAD
#undef AVXSTORE
#undef AVXADD
#endif
#define AVXLOAD LOADuDOUBLE
#define AVXSTORE STOREuDOUBLE
#define AVXADD ADDDOUBLE
gV5_start(MY_VARIANT, double, LongDouble, double, Doubles, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
gV5_CreateHash(double, LongDouble, double, Doubles, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_MainLoop(double, LongDouble, double, Doubles, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_SumUp(double, LongDouble, double, Doubles, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_Sort(double, LongDouble, double, Doubles, double)
  


#endif
