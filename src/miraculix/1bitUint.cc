/*
 Authors 
 Guido Moerkotto
 maintained and modified by Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Guido Moerkotte, Martin Schlather

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


#define BitsPerCode 1
#define MY_VARIANT 32

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"

#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "TemplateUint.h"
#include "MX.h"
#include "bitUint.h"
#include "1bit.h"



Long Lda1Bit(Long snps, Long LDAbitalign){
  return (LDAbitalign ? 2 : 1) * LDAalignBitSnps(snps * BitsPerCode);
}

#define coding1_end(UINT) \
  for (Long i=0; i<deltaIndiv; i++) {\
    unit_t  *code = ans + i * ldAns,\
      *dode = Dode + i * ldAns;\
    UINT *mm = (UINT*) (M + i * ldM);	\
    bool w = false;\
    for (Long j=0; j<allUnits; j++) {\
      unit_t C = 0,\
	D=0;\
      int end = j == allUnitsM1 ? rest : BitsPerUnit;		\
      for (int shft = 0; shft < end; mm++, shft+=BitsPerCode) {	\
	w |= *mm > 2;\
	C |= (unit_t) ((*mm == 1) << shft);	\
	D |= (unit_t) ((*mm == 2) << shft);	\
      }\
      code[j] = C;\
      dode[j] = D;\
    }\
    total += w;\
  }\
  if (total > 0) ERR0("1bit#64: values outside 0,1,2");\
}

coding_start(_4Byte,1)
unit_t  *Dode = ans + ldAns * individuals;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
coding1_end(_4Byte)


coding_start(_1Byte,1)
unit_t  *Dode = Ans + ldAns * individuals;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
coding1_end(_1Byte)



sumGeno_start(1)
  unit_t *T = code + lda * individuals;
  for (Long i=0; i<individuals; i++) {
    unit_t *c = code + i * lda;
    unit_t *TT = T + i * lda;
    Ulong sum = 0UL,
      twice = 0UL;
    for (Long j=0; j<units; j++) { 				
      unit_t s = c[j],
	t = TT[j];						
      for (Long u=0; u<CodesPerUnit; u++) {			
	sum += s & CodeMask;
	twice += t & CodeMask;
	s >>= BitsPerCode;					
	t >>= BitsPerCode;					
      }					
    }
    if (sums != NULL) sums[i] = sum + (twice << 1);
    total += sum + (twice << 1);
  }						  
  return total;			  
}


#define get_matrix1_end(UINT)			\
  for (Long i=0; i<individuals; i++) {\
    UINT *a = (UINT*) (Ans + i * ldAns);	\
    unit_t *C = code + i * lda,	       \
      *D = D0 + i * lda;				\
    for (Long j=0; j<=allUnitsM1; j++, C++, D++) {	\
      unit_t c = *C,					\
	d = *D;								\
      int end_k = j==allUnitsM1 ? rest_k : CodesPerUnit;		\
      for (int k=0; k < end_k; k++) {					\
	*(a++) = (c & CodeMask) + ((d & CodeMask) << 1);		\
	c >>= BitsPerCode;			    \
	d >>= BitsPerCode;			    \
      }						    \
    }						    \
  }						    \
  }


get_matrix_start(_4Byte,1)			   
  unit_t *D0= code + individuals * lda;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix1_end(_4Byte)


get_matrix_start(_1Byte,1)			   
  unit_t *D0= code + individuals * lda;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
get_matrix1_end(_1Byte)



void zeroGeno1(SEXP SxI, Uint *Snps, Uint lenSnps, Uint *Indiv, Uint lenIndiv,
	       basic_options *opt) {
  int cores = GreaterZero(opt->cores);
  unit_t *code = Align(SxI, opt);
  Long *info = GetInfo(SxI),
    SxIsnps = info[SNPS],
    SxIindiv = info[INDIVIDUALS],
    lda = info[LDA];
  ASSERT_LITTLE_ENDIAN;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<lenIndiv; i++) {
    if (Indiv[i] >= SxIindiv) continue;
    unit_t *C1 = code + Indiv[i] * lda,
      *C2 = C1 + lda * SxIindiv;
    for (Long s=0; s<lenSnps; s++) {
      if (Snps[s] >= SxIsnps) continue;
      Long j;
      unit_t blend;
      PositionCode(Snps[s], &j, &blend);
      C1[j] &= ~blend;
      C2[j] &= ~blend;
    }
  }
}


// LDA in Units!

