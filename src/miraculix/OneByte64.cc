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


#define BitsPerCode 8
#define MY_VARIANT 64
#define PlainInteger32 1

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "MXinfo.h"



void OneByteHaplo2geno64(unit_t *H,
			 Long VARIABLE_IS_NOT_USED snps, Long individuals,
			 Long lda,
			 int VARIABLE_IS_NOT_USED cores,
			 unit_t *Ans,
			 Long ldAns) {
  // works also for 128,256,512; to do
  assert((intptr_t) H % (MY_VARIANT / BitsPerByte) == 0 &&
	 (intptr_t) Ans % (MY_VARIANT / BitsPerByte) == 0);
  if (H == Ans && ldAns > lda) ERR0("1byteH to geno not possible"); // to do?
  const Long total = lda * individuals;
   const Long blocks = total / UnitsPerBlock; // OK
  BlockType
    *H1 = (BlockType*) H,
    *H2 = (BlockType*) (H + total),
    *G = (BlockType*) Ans;

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif  
  for (Long i=0; i<blocks; i++)
    STORE(G + i, ADD64(LOAD(H1 + i), LOAD(H2 + i)));
  for (Long i=blocks * UnitsPerBlock; i<total; i++) // OK
    Ans[i] = H[i] + H[i+total];
}

void FourHaplo2geno64(unit_t *H, Long snps, Long individuals, Long lda,
		      int cores, unit_t *Ans, Long ldAns) {
  // alignment not garanteed!
  if ((intptr_t) H % (MY_VARIANT / BitsPerByte) == 0 &&
      (intptr_t) Ans % (MY_VARIANT / BitsPerByte) == 0) {
    OneByteHaplo2geno64(H, snps, individuals, lda, cores, Ans, ldAns);
    return;
  }

if (H == Ans && ldAns > lda) ERR0("1byteH to geno not possible"); // to do?
  const Long total = lda * individuals;
  const Long blocks = total / UnitsPerBlock; // OK
  BlockType
    *H1 = (BlockType*) H,
    *H2 = (BlockType*) (H + total),
    *G = (BlockType*) Ans;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif  
  for (Long i=0; i<blocks; i++)
    STOREU(G + i, ADD64(LOADU(H1 + i), LOADU(H2 + i)));
  for (Long i=blocks * CodesPerBlock; i<total; i++)
    Ans[i] = H[i] + H[i+total];
}

