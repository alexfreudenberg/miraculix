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

// _MM_ALIGN16 Uint ccc;
#define BitsPerCode 3
#define MY_VARIANT 32
#define PlainInteger32 1
#define MY_LDABITALIGN 64
#define NO_SSE2 1 // needed as MY_LDA must be given

#define BytesPerPartUnit 2


//#include <stdio.h>
#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "MX.h"
#include "haplogeno.h"
#include "TemplateUint.h"
#include "Haplo.h"
#include "bitUint.h"

#define nr_genotypes 3
#define nr_results 4
static Long rev_geno_code[7] = {0, NA_LONG, NA_LONG, 1, NA_LONG, NA_LONG, 2}; 
static Uint geno_code[nr_genotypes] = {0, 3, 6}, // 0, 1, 2
  result_code[nr_results] = {0, 3, 2, 6};
static Uint result_value[nr_results] = {0, 1, 2, 4};


#define TABLE_SIZE  two_codingBitsPerPartUnit
table_type *TABLE3 = NULL;

void initiate_tableI(table_type **table, int tableSize,
		     int codesPerPartUnit, int bitsPerCode,
		     Uint *resultCode,
		     Uint *resultValue, int NrResults) { // hash table

  int d;
  int *nx = (int *) CALLOC(codesPerPartUnit, sizeof(int));
  
  *table = (table_type*) CALLOC(tableSize, sizeof(table_type));
  table_type* pt = *table;

  while (true) {
    int shift;
    Uint value = 0,
      sum = 0;
    
    for (d=shift=0; d<codesPerPartUnit; d++, shift+=bitsPerCode) {
      value |= resultCode[nx[d]] << shift;
      //    printf("sh=%d ",shift);
      sum += resultValue[nx[d]];
    }

       
    assert(value < (Uint) tableSize && pt[value] == 0);
    pt[value] = (char) sum;
    d = 0;
    nx[d]++;
    while (nx[d] >= NrResults) {
      nx[d] = 0;
      if (++d >= codesPerPartUnit) break;
      nx[d]++;
    }
    if (d >= codesPerPartUnit) break;
  }
  FREE(nx);
}



void Init3() {
  assert(PartUnitsPerBlock == 2); // OK
   initiate_tableI(&TABLE3, TABLE_SIZE, CodesPerPartUnit,
		  BitsPerCode, result_code, result_value, nr_results);
}


#define coding3_start(UINT)						\
  coding_header(UINT,3) {						\
  if (start_snp % CodesPerUnit != 0) BUG;				\
  unit_t *a = Ans + start_snp / CodesPerUnit; 

coding3_start(_4Byte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif								
  for (Long i=start_individual; i<end_individual; i++) {
    UNIT_CODING(_4Byte, M, ldM,
	       start_snp, end_snp, start_individual, end_individual,
	       a, ldAns, FROMINPUT)    
      } // see bit23intern.h
}

coding3_start(_1Byte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif								
for (Long i=start_individual; i<end_individual; i++) {
  UNIT_CODING(_1Byte, M, ldM,
	      start_snp, end_snp, start_individual, end_individual, a,
	      ldAns, FROMINPUT)
    }
}


Long getValue3(unit_t *M, Long S) {
  return rev_geno_code[ExtractCode(M, S)];
}



void allele_sum3(unit_t *Code, Long snps, Long individuals, Long lda,
		  basic_options VARIABLE_IS_NOT_USED *opt, Ulong *ans) {
  assert(TABLE3 != NULL);  
  MEMSET(ans, 0, sizeof(*ans) * snps);					
  for (Long i=0; i<individuals; i++) {					
    unit_t *code = (unit_t*) (Code + i * lda);				
    for (Long s=0; s<snps; s++) ans[s] += getValue3(code, s);		
  }									
}


#define get_matrix3_end(UINT)					\
  for (Long i=0; i < individuals; i++) {			\
    UINT *A = (UINT*) (Ans + i * ldAns);		\
    unit_t *c = code + i * lda;		\
    for (Long s=0; s<snps; s++) A[s] = (UINT) getValue3(c, s);	\
  }								\
}

get_matrix_start(_4Byte,3)
  assert(TABLE3 != NULL);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix3_end(_4Byte)

  
get_matrix_start(_1Byte,3)
  assert(TABLE3 != NULL);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif	  
get_matrix3_end(_1Byte)




SEXP get_matrix23_start(Long snps, Long individuals, 
		      SEXP VARIABLE_IS_NOT_USED G) {
  return allocMatrix(INTSXP, snps, individuals);
}


void zeroGeno3(SEXP SxI, Uint *Snps, Uint lenSnps, Uint *Indiv, Uint lenIndiv,
	       basic_options *opt) {
  int cores =  GreaterZero(opt->cores);
  unit_t  *M = Align(SxI, opt);
  Long
    *info = GetInfo(SxI),
    SxIsnps = info[SNPS],
    SxIindiv = info[INDIVIDUALS],
    lda = info[LDA];
  ASSERT_LITTLE_ENDIAN;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
 for (Long i=0; i<lenIndiv; i++) {					
    if (Indiv[i] >= SxIindiv) continue;
    unit_t *Ma = M + Indiv[i] * lda ;             	
    for (Long s=0; s<lenSnps; s++) {					
      if (Snps[s] >= SxIsnps) continue;     
      Long j;
      unit_t mask;
      PositionCode(Snps[s], &j, &mask);
      Ma[j] &= ~mask;
    }							
  }							       
}

			

sumGeno_start(3)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) { /* ok */
    unit_t *pcode = code + lda * i;
    Ulong sum = 0L,						
      counter = 0;						
   for (Long j=0; j<units; j++) {				
      unit_t s = pcode[j];						
      for (Uint u=0; u<CodesPerUnit; u++) {			
	sum += rev_geno_code[s & CodeMask];
	s >>= BitsPerCode;					
	if (++counter >= CodesPerPartUnit) {			
	  s >>= deltaBitsPartUnit;				
	  counter = 0;						
	}							
      }								
    }
    if (sums != NULL) sums[i] = sum;
    total += sum;
  }						  
  return total;						  
}

  

void TwoBithaplo2geno3(unit_t *M, Long snps, Long individuals, Long oldlda,
		       int VARIABLE_IS_NOT_USED cores,
		       unit_t *Ans,
		       Long ldAns) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) { 
    UNIT_CODING(unit_t, M, oldlda,
		0, snps, 0, individuals, Ans,
		ldAns, FROMHAPLO);
  }
}
