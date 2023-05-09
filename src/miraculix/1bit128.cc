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
#define MY_VARIANT 128

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
//
//

#if defined SSSE3 

ASSERT_SIMD(1bit128, ssse3);

#include "MX.h"
#include "haplogeno.h"
#include "intrinsicsIntern.h"
#include "1bitIntern.h"



void trafo1Geno2Geno128(unit_t *OLD, Long snps, Long individuals,
			Long oldLDA,  unit_t *Ans, Long ldAns) {
  assert(2 * oldLDA == ldAns);
  
  Long blocks = Blocks(snps);

  if (OLD == Ans) {
    ERR0("not programed yet");
    // (Braucht nur 1/4 Speicherplatz extra)
    // (i)drittes Viertel rauskopieren + 1-2 Spalten mehr vorne freilasen
    // (ii) haelfte durcharbeitne
    // (iii) drittel Vierel aus speicher + in die Spalte(n) davor im
    //       Speicher viertel Viertel spaltenweise reinkopieren
    // (iv) letztes Viertel aus dem Speicher mit den 1-2 Spalten frueher
    //      abarbeiten.
  }	
  
   
  const BlockType 
    O2F2 = SET32(0x00FF00FF),
    F2O2 = SET32(0xFF00FF00),
    
    O4F4 = SET32(0x0000FFFF),
    F4O4 = SET32(0xFFFF0000);

  
#if ! defined PlainInteger32
  const BlockType O8F8 = SET64(0x00000000FFFFFFFF);
#endif	

#if defined AVX512F
  const BlockType cross_idx =
    _mm512_set_epi32(TAKE_B + 14, TAKE_B + 12, TAKE_B + 10, TAKE_B + 8,
		     TAKE_B + 6, TAKE_B + 4, TAKE_B + 2, TAKE_B + 0,
		     14, 12, 10, 8, 6, 4, 2, 0);
#elif defined AVX2
  const BlockType cross_idxA = _mm256_set_epi32(3, 3, 2, 2, 1, 1, 0, 0);
  const BlockType cross_idxC = _mm256_set_epi32(7, 7, 6, 6, 5, 5, 4, 4);
  const int halfreg_idx = 0 + (2 << 4);
#elif defined SSE2
  const int cross_idxA =  0 + (1 << 4);
  const int cross_idxC = 2 + (3 << 4);
#elif defined PlainInteger64 || defined PlainIteger32
  #define cross_idxA 0
  #define cross_idxC 1
#endif

 

  for (Long j=0; j<individuals; j++) {
    unit_t *New = Ans + j * ldAns,
      *halfmatrix = New + individuals * ldAns,
    *Old = OLD + j * oldLDA;
    
    for (Long i=0; i<blocks; i++) {	
      const BlockType a0 = LOAD((BlockType0*) Old);
      const BlockType b0 = LOAD((BlockType0*) halfmatrix);
      
      const BlockType a1 = CROSSLINE_PERMUTE(a0, cross_idxA);
      const BlockType b1 = CROSSLINE_PERMUTE(b0, cross_idxA);
      const BlockType c1 = CROSSLINE_PERMUTE(a0, cross_idxC);
      const BlockType d1 = CROSSLINE_PERMUTE(b0, cross_idxC);

    
#if defined PlainInteger32
      const BlockType a5 = a1,
	b5 = b1,
	c5 = c1,
	d5 = d1;
#else
      AND4_C(2, 1, O8F8);
      AND4_C(3, 2, F4O4);
      SHL32_4(4, 3, 16);
      OR4(5, 4, 2);    
#endif
    
      AND4_C(6, 5, O4F4);
      AND4_C(7, 6, F2O2);
      SHL32_4(8, 7, 8);
      OR4(9, 8, 6);
      
      AND4_C(10, 9, O2F2);
      AND4_C(11, 10, F1O1);
      SHL32_4(12, 11, 4);
      OR4(13, 12, 10);
      
      AND4_C(14, 13, O1F1);
      AND4_C(15, 14, E2N2);
      SHL32_4(16, 15, 2);
      OR4(17, 16, 14);
     
      AND4_C(18, 17, N2E2);
      AND4_C(19, 18, E1N1);
      SHL32_4(20, 19, 1);
      OR4(21, 20, 18);
            
      AND4_C(22, 21, N1E1);
      const BlockType c23 = SHL32(c22, 1);
      const BlockType d23 = SHL32(d22, 1);
      const BlockType a23 = OR(a22, c23);
      const BlockType b23 = OR(b22, d23);
       
      *((BlockType0*) New) = a23;
      //      *((BlockType0*) (New + UnitsPerBlock)) = c23; // OK
     *((BlockType0*) (New + UnitsPerBlock)) = b23; // OK

      Old += UnitsPerBlock;// OK
      halfmatrix += UnitsPerBlock;// OK
      New += 2 * UnitsPerBlock;// OK
    } // for i
  } // for j
}

SCALAR012(128)	




#else // !defined AVX
#include "Template.h"
#include "1bit.h"
#include "avx_miss.h"


  DeclareVersion1(128,Su,Sv)
  SIMD_MISS(1bit128, ssse3);

#endif
