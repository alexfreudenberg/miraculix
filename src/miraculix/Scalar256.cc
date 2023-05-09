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

#define NO_AVX512 1
#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "Scalar.h"


#if defined AVX2
Long scalarIntAVX2(int * V1, int * V2, int N) { return scalarInt8by8(V1, V2, N); }
Ulong scalarUintAVX2XXX(Uint * V1, Uint * V2, Long N){
  BUG;
  // arg primitiv -- FEHLER?!: MULT32 liest auch nur 128 BIT !!
  Long steps = N / UnitsPerBlock; // OK
  BlockType0    
    *v1 = (BlockType0 *) V1,
    *endv1 = v1 + steps,
    *v2 = (BlockType0 *) V2;
  UnionType sum; // OK
  sum.vi=ZERO();

  Uint zaehler = 0;
  for (; v1 < endv1; v1 ++, v2 ++) {
    zaehler++;
    BlockType dummy, s1, s2;
    s1=LOADU(v1);
    s2=LOADU(v2);    
    dummy = MULT32(s1, s2);
    sum.vi = ADD32(sum.vi, dummy);
  }
  
  Ulong totalsum = sum.u32[0] + sum.u32[1] + sum.u32[2] + sum.u32[3] +
                  sum.u32[4] + sum.u32[5] + sum.u32[6] + sum.u32[7];
  Uint *end = V1 + N;  				
  // V1 += steps * Uni tsPerBlock;  wrong?!
  //  V2 += steps * Un itsPerBlock; wrongh?!
  for (; V1 < end; V1++, V2++) totalsum += V2[0] * V1[0];			
  return totalsum;								
}


#else// !defined AVX

Uint scalarUint8by8(Uint * V1, Uint * V2, Uint N);
//Ulong scalarUintAVX2(Uint * V1, Uint * V2, Uint N) { return scalarUint8by8(V1, V2, N); }
Long scalarIntAVX2(int * V1, int * V2, int N) { return scalarInt8by8(V1, V2, N); }
#endif
