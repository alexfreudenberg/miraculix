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
#define MY_VARIANT 32
#define BytesPerPartUnit 2


#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "TemplateUint.h"
#include "Haplo.h"

#include "MX.h"
#include "1bit.h"
#include "bitBase.h"
#include "bitUint.h"
#include "2bit.h"
#include "plink.h"


#define nr_genotypes 3
static Uint geno_code[nr_genotypes] = {0, 1, 2};
  

#define nr_resultsAND 3
#define nr_resultsOR 4
static Uint
  result_codeAND[nr_resultsAND] = {0, 1, 2},
  result_codeOR[nr_resultsOR] = {0, 1, 2, 3},
  result_valueAND[nr_resultsAND] = {0, 1, 4},
  result_valueOR[nr_resultsOR] = {0, 0, 0, 2};


table_type *TABLE2AND = NULL,
  *TABLE2OR = NULL;

Long Lda2Bit(Long snps, Long LDAbitalign){
  // to get compatibility with 1Bit concerning memory allocation
 // 1Bit has two columns with each empty space above
 // 2bit has these two columns joint with a single empty space
 // at the end
  return LDAalignBitSnps(snps * BitsPerCode); 
}


void Init2() {
  assert(PartUnitsPerBlock == 2); // OK
  initiate_tableI(&TABLE2AND, TABLE2_SIZE, CodesPerPartUnit, BitsPerCode,
		  result_codeAND, result_valueAND, nr_resultsAND);
  initiate_tableI(&TABLE2OR, TABLE2_SIZE, CodesPerPartUnit, BitsPerCode,
		  result_codeOR, result_valueOR, nr_resultsOR);
}

void allele_sum2(unit_t *Code, Long snps, Long individuals, Long lda,
		  basic_options VARIABLE_IS_NOT_USED *opt, Ulong *ans) {
  MEMSET(ans, 0, sizeof(*ans) * snps);					
  for (Long i=0; i<individuals; i++) {					
    unit_t *code = (unit_t*) (Code + i * lda);				
    for (Long s=0; s<snps; s++) ans[s] += ExtractCode(code, s);		
  }									
}

sumGeno_start(2) 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
sumGeno_fctn()
 


#define TWOBITTRAFO 

get_matrix_start(_4Byte, 2)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix_end(_4Byte,TWOBITTRAFO)   

get_matrix_start(_1Byte, 2)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix_end(_1Byte,TWOBITTRAFO)   
	
coding_start(_4Byte,2)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
coding_end(_4Byte,*mm)
 
coding_start(_1Byte, 2)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
coding_end(_1Byte,*mm)



coding2trans_start(_4Byte,2)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif								
coding2trans_do(_4Byte,HUMAN2TWOBIT)

coding2trans_start(_1Byte,2)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding2trans_do(_1Byte,HUMAN2TWOBIT)



void zeroGeno2(SEXP SxI, Uint *Snps, Uint lenSnps, Uint *Indiv, Uint lenIndiv,
	       basic_options *opt) {			  
  int cores = GreaterZero(opt->cores);
  Long
    *info = GetInfo(SxI),
    SxIsnps = info[SNPS],
    SxIindiv = info[INDIVIDUALS];
  Long lda = info[LDA];
  unit_t *C0 = Align(SxI, opt);
    
  ASSERT_LITTLE_ENDIAN;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
 for (Long i=0; i<lenIndiv; i++) {
    if (Indiv[i] >= SxIindiv) continue;
    unit_t *C1 = C0 + Indiv[i] * lda;
    for (Long s=0; s<lenSnps; s++) {
      if (Snps[s] >= SxIsnps) continue;
       Long j;
       unit_t blend;
       PositionCode(Snps[s], &j, &blend);
      C1[j] &= ~blend;
    }
  }
}




void TwoBithaplo2geno2(unit_t *M, Long snps, Long individuals, Long ldM, 
		       int VARIABLE_IS_NOT_USED cores,  unit_t *Ans, Long ldAns) {
  // note that M == Ans is allowed !!
  Long start_individual=0,
    end_individual=individuals;    
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=start_individual; i<end_individual; i++) {
    UNIT_CODING(unit_t, M, ldM,
		0, snps, start_individual, end_individuals,
		Ans, ldAns, FROMHAPLO);
  }
}



vectorGeno_start(Ulong, Ulong, Ulong, 2);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
vectorGeno_fctn(Ulong, Ulong, Ulong, TWOBITTRAFO)


vectorGeno_start(double, double, double, 2);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
vectorGeno_fctn(double, double, double, TWOBITTRAFO)


vectorGeno_start(LongDouble, LongDouble, LongDouble, 2);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
vectorGeno_fctn(LongDouble, LongDouble, LongDouble, TWOBITTRAFO)




#define genoVector2_start(TYPE,TRDTYPE, ENDTYPE)			\
  genoVector_header(TYPE,TRDTYPE, ENDTYPE, 2) {				\
  int cores = GreaterZero(opt->cores);					\
  Long units = Units(rows),						\
    unitsM1 = units - 1L,						\
    size = rows * repetV,						\
    rest = rows - unitsM1 * CodesPerUnit,				\
    slicesM1 = cores - 1L,						\
    slicesCPU = (slicesM1 + 1) * CodesPerUnit;				\
  for (Long s=0; s<size; Ans[s++] = 0.0);				\
  if (coding != TwoBitGeno && coding != TwoBitGenoTransposed && \
      coding != Plink && coding != PlinkTransposed) BUG;		\
  ENDTYPE f1 = coding == TwoBitGeno ? 1 : 0;				\
  ENDTYPE f2 = coding == TwoBitGeno ? 2 : 1;				\
  ENDTYPE f3 = coding == TwoBitGeno ? (ENDTYPE) RF_NA : 2 /* =0 if END = Long*/


#define genoVector2_fctn(TYPE,TRDTYPE, ENDTYPE)			\
  for (Long r=0; r < repetV; r++) {					\
    for (Long metaB=0; metaB <= slicesM1; metaB++) {			\
      Long /* NOT Long !! */						\
	bstart = ((metaB * units) / slicesCPU) * CodesPerUnit,		\
	bendM1 = metaB == slicesM1					\
	? unitsM1							\
	: (((metaB+1L)*units) / slicesCPU) * CodesPerUnit -1L;/*-1 possible!*/ \
      for (Long i=0; i<cols; i++) {					\
	unit_t *pC = code + i * lda;					\
	ENDTYPE f0 = (ENDTYPE) V[i + r * ldV],			\
	  *a = Ans + r * ldAns,						\
	  factor[4] = {0, f1 * f0, f2 * f0, f3 * f0};			\
	ENDTYPE aPF = 0;						\
	if (pseudoFreq != NULL) aPF = 2 * f0 * (ENDTYPE) pseudoFreq[i];	\
	for (Long b=bstart; b<=bendM1; b++) {				\
	  unit_t c = pC[b];						\
	  Long endfor = b < unitsM1 ? CodesPerUnit : rest;		\
	  for (Long j=0; j<endfor; j++) {				\
	    a[j + b * CodesPerUnit] += (ENDTYPE) (factor[c & CodeMask] - aPF); \
	    c >>= BitsPerCode;						\
	  }								\
	}								\
      }									\
    }									\
  }								     

genoVector2_start(Ulong, Ulong, Ulong);
#ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
genoVector2_fctn(Ulong, Ulong, Ulong)
}

genoVector2_start(double, double, double);
//printf("geno\n");
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
genoVector2_fctn(double, double, double)
}

genoVector2_start(LongDouble, LongDouble, LongDouble);
//printf("geno\n");
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif	
genoVector2_fctn(LongDouble, LongDouble, LongDouble)
}



#define gV2_Select(TYPE)					 \
  gV_vG_start(TYPE, 2);						 \
  if (gV) vD = genoVector2_##TYPE;				 \
  else vD = coding == Plink ? vectorGenoPlink_##TYPE : vectorGeno2_##TYPE; \
  gV_vG_end;							\
}

gV2_Select(LongDouble)
gV2_Select(Ulong)


gV_vG_start(double, 2);
//printf("%s v=%d %d %d rows=%ld cols=%ld\n", CODING_NAMES[coding], variant, getAttribPointer(SxI, Next)==R_NilValue, gV, rows, cols);

assert(variant != 257 || !gV);
if (gV) vD = variant < 256 ? genoVector2_double : genoVector2v256_double; 
  else {
    if (variant < 256 || SubstractMeanSxI)
      vD = coding == Plink ? vectorGenoPlink_double : vectorGeno2_double;	
    else vD = variant == 257 && (coding == Plink ||
			       coding == PlinkTransposed ||
			       coding == OrigPlink ||
			       coding == OrigPlinkTransposed
			       )
	   ? vectorGenoPlinkMatrix256_double
	   : vectorGeno2v256_double;
  }
  gV_vG_end;
}


sparseTGeno_header(2Bit) {
  BUG;
  // unsinns aufruf!!:
   sparseTGenoPlink(code, nrow_code, ncol_code,
		    ldaInByte, coding,
		    valueB, nIdx, rowIdxB, colIdxB, tSparse,
		    global, utils, 
		    Ans, ldAns);
  
}
 
