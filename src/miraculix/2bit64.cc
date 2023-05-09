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
#define MY_VARIANT 64

//#include <stdio.h>

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "intrinsicsIntern.h"
#include "xport_import.h"
#include "options.h"
#include "MX.h"
#include "haplogeno.h"
#include "Files.h"
#include "bitBase.h"

#ifdef AVX
#error AVX defined
#endif


extern table_type *TABLE2AND, *TABLE2OR;


FileBinary(64)
 
void crossprod2v64(unit_t * M, Long snps, Long individuals, Long lda,
		   int VARIABLE_IS_NOT_USED cores, double *A){
  Long ldAns = individuals;
  //  printf("TABLE2_SIZE %ld\n", TABLE2_SIZE);
  assert(TABLE2_SIZE == 65536);
  assert(TABLE2AND != NULL);
  Long blocks = Blocks(snps);
  table_type tableAnd[TABLE2_SIZE], tableOr[TABLE2_SIZE];

  MEMCOPY(tableAnd, TABLE2AND, TABLE2_SIZE * sizeof(table_type));
  MEMCOPY(tableOr, TABLE2OR, TABLE2_SIZE * sizeof(table_type));
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic, 20)  
#endif	  
  for (Long i=0; i<individuals; i++) {   
    unit_t * Mi = M + i * lda,
      *Mj = Mi;
    double *ans = A + i,
      *res_ptr_r = ans + i * ldAns;
    for (Long j  =  i; j<individuals; j++, res_ptr_r++, Mj += lda) {
      Long sum = 0;
      BlockType *MJ = (BlockType0*) Mj,
	*MI = (BlockType0*) Mi;
      for (Long s=0; s<blocks; s++, MJ++, MI++) {
	block_compressed x;
	x.x = *MJ & *MI;
	sum += tableAnd[x.b[0]] + tableAnd[x.b[1]] +
	  tableAnd[x.b[2]] + tableAnd[x.b[3]];

 	x.x = *MJ | *MI;
	sum += tableOr[x.b[0]] + tableOr[x.b[1]]
	  + tableOr[x.b[2]] + tableOr[x.b[3]];
      }
      ans[j * individuals] = *res_ptr_r = (double) sum;
    }
  }
}

 	
