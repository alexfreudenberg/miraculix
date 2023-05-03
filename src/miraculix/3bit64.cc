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



// gcc -mavx -o hello_avx hello_avx.c

// _MM_ALIGN16 Uint ccc;
#define BitsPerCode 3
#define MY_VARIANT 64
#define PlainInteger64 1
#define MY_LDABITALIGN 64
#define NO_SSE2 1 // needed as MY_LDA must be given


//#include <stdio.h>
#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "MX.h"
//#include "Haplo.h"
#include "haplogeno.h"
#include "bitUint.h"


extern table_type *TABLE3;

Long CodesPerBlock3() { return CodesPerBlock; }

Long Lda3(Long snps, Long LDAbitalign) {
  return LDAalignBlocks(Blocks(snps));
}

Long BitsPerCode3() { return BitsPerCode; }




void crossprod3(unit_t * M, Long snps, Long individuals, Long lda,
		int VARIABLE_IS_NOT_USED cores, double *A) {
  Long ldAns = individuals;
  Long blocks = Blocks(snps);
  table_type *table = TABLE3;
  assert(table != NULL);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20) 
#endif	  
  for (Long i=0; i<individuals; i++) {   
    unit_t * Mi = M + i * lda,
      *Mj = Mi;    
    double *ans = A + i,
      *res_ptr_r = ans + i * ldAns;
    for ( Long j  =  i; j<individuals; j++, res_ptr_r++, Mj += lda) {
      Long sum = 0;
      BlockType *MJ = (BlockType0*) Mj,
	*MI = (BlockType0*) Mi;
      for (Long s=0; s<blocks; s++, MJ++, MI++) {
	block_compressed x;
	x.x = *MJ & *MI;
	sum+= (table[x.b[0]] + table[x.b[1]] + table[x.b[2]] + table[x.b[3]]);
      }
      ans[j * individuals] = *res_ptr_r = (double) sum;
    }
  }
}


