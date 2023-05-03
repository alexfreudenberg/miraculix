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



#define BitsPerCode 32
#define MY_VARIANT 32

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "TemplateUint.h"
//#include "TwoBitHaplo.h"
#include "MX.h"
#include "Scalar.h"
#include "Haplo.h"


Long LdaPlain(Long snps, Long LDAbitalign) {
  //   return snps * (LDAbitalign ? LDAbitalign : MY_LDABITALIGN) / MY_LDABITALIGN;
  return LDAalignBitSnps(snps * BitsPerCode); 
}


Long LdaPlainDoubled(Long snps, Long LDAbitalign) {
  //   assert(LDAbitalign || snps == 1);
  // return LDAbitalign ? 2 * snps * LDAbitalign / MY_LDABITALIGN : 1;
  return (LDAbitalign == 0 ? 1 : 2) * LDAalignBitSnps(snps * BitsPerCode); 
}


#define coding_plain_end(UINT)					\
  for (Long i=start_individual; i<end_individual; i++) {	\
    UINT *pM = (UINT*) (M + (i - start_individual) * ldM);	\
    _4Byte *pAns = (_4Byte*) (Ans + i * ldAns);			\
    for (Long s=start_snp; s<end_snp; pM++) {			\
      pAns[s++] = (_4Byte) *pM;					\
    }								\
  }							       


coding_header(_4Byte,Plain) {				     
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  coding_plain_end(_4Byte)
    }

coding_header(_1Byte,Plain) {			     
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  coding_plain_end(_1Byte)
    }



get_matrix_header(_4Byte,Plain) {
  Long bytes = MIN(lda, ldAns) * (BitsPerCode / BitsPerByte);
  if (lda == ldAns) {
    if (code != Ans) MEMCOPY(Ans, code, (Long) individuals * bytes);
  } else {
    for (Long i=0; i<individuals; i++) {
      unit_t *m = (code + i * lda),
	*a =  Ans + i * ldAns;
      MEMCOPY(a, m, bytes);
    }
  }
}

get_matrix_header(_1Byte,Plain) {
  for (Long i=0; i<individuals; i++) {
    _4Byte *m = (_4Byte*) (code + i * lda);
    _1Byte *a = (_1Byte*) (Ans + i * ldAns);
    for (Long j=0;  j<snps; j++) a[j] = (_1Byte) m[j];
  }
}




sumGeno_start(Plain)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) {
    unit_t *pcode = code + i * lda;
    Long sum = 0L;					  
    for (Long j=0; j<snps; j++, pcode++) sum += *pcode;
    if (sums != NULL) sums[i] = sum;
    total += sum;
  }
  return total;						  
}

void TwoBithaplo2genoPlain(unit_t *code, Long snps, Long individuals, Long lda,
			   int VARIABLE_IS_NOT_USED cores,
			   unit_t *Ans, Long ldAns) {
  MEMSET(Ans, 0, totalMem(snps, individuals, lda, FourByteGeno,
			  false) * BytesPerUnit);
  assert(BitsPerCode == 32);
#define ByteCodeMask 3

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) {
    unit_t
      *pC = code + lda * i,
      *a = Ans + i * ldAns;
    for (Long s = 0; s < snps; s++)
      a[s] = GetTwoBitHaplo(pC, s);
  }
}

void crossprod_Plain(unit_t * Code, Long snps, Long individuals,
		     Long VARIABLE_IS_NOT_USED lda,
		     int cores, double *A){
  stopIfNotInt(snps | individuals);
  stopIfNotSame(_4Byte, int);
  matmulttransposedInt((int*) Code, (int*) Code, A,
		       (int) snps, (int) individuals,
		       (int) individuals, cores); 
}
