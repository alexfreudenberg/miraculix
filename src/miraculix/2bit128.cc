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
#define MY_VARIANT 128

#define NO_SSE41 


#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "Files.h"
#include "Template.h"
#include "2bit.h"

#if defined SSSE3 

ASSERT_SIMD(2bit128, ssse3);

#include "MX.h"
#include "haplogeno.h"
#include "intrinsicsIntern.h"
#include "bitBase.h"
#include "2bitIntern.h"



sumGeno_start(2v128)
  return sumGenoIntern(code, snps, individuals, lda, cores, sums);
}

	
void TwoBithaplo2geno128(unit_t *M, Long snps, Long individuals, Long lda, 
			 int cores, unit_t *Ans,Long ldAns) {
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  TwoBithaplo2genoIntern(M, snps, individuals, lda, cores, Ans, ldAns);
}


/*

   for (int j=0; j<x15; j++) {
      const BlockType a0 = LOAD(xx++);
      const BlockType b0 = LOAD(yy++);

      BlockType z, S, T, y0, z0, y1;  // + and2scalar + or2scalar + a + b (+ sum)
      z0 = AND0(a0, b0);
      y0 = XOR0(a0, b0);
      z = AND0(y0, E1N1); // 4 + 1
      z = SHR32_0(z, 1);
      z = AND0(z, y0); // 3+ 1
      y1 = SHL32_0(z, 1); // 4+ 1
      z = OR0(z, y1);
      z = OR0(z, z0); // 2+ 1
      
      y1 = AND0(z, LH); 
      z = SHR32_0(z, 4);
      z = AND0(z, LH); // 3+ 1
	      
      S = SHUFFLE8_0(to_scalar, y1); // 4+ 1
      T = SHUFFLE8_0(to_scalar, z); // 4+ 1
      sum = ADD8_0(sum, S); 
      sum = ADD8_0(sum, T); 
    }
    ans += ADDUPVECTOR(sum); 

*/



/*

static Uint scalar(BlockType0 *x, BlockType0 *y, Long blocks) {
  Uint ans = 0;
  //#define x15 (255L / 16L)  // 16L = 4 * 2^2 (a byte full of 1's)
#define x15 (255L / 16L)  // warum ist dies schneller ?????
  Long
    bundles = blocks  / x15;

  for (Long i=0; i<bundles; i++) { 
    BlockType sum,
      *xx = x + i * x15,
      *yy = y + i * x15;
    sum = ZERO();
    for (int j=0; j<x15; j++) {
      const BlockType a0 = LOAD(xx++);
      const BlockType b0 = LOAD(yy++);
      const BlockType z0 = AND0(a0, b0);
      const BlockType z1 = XOR0(a0, b0);
      
      const BlockType z2 = AND0(z1, E1N1); // 4 + 1
      const BlockType z3 = SHR32_0(z2, 1);
      const BlockType z4 = AND0(z3, z1); // 3+ 1
      const BlockType z5 = SHL32_0(z4, 1); // 4+ 1
      const BlockType z6 = OR0(z4, z5);
      const BlockType z7 = OR0(z6, z0); // 2+ 1
      
      const BlockType z8 = AND0(z7, LH); 
      const BlockType z9 = SHR32_0(z7, 4);
      const BlockType z10 = AND0(z9, LH); // 3+ 1
	      
      const BlockType z11 = SHUFFLE8_0(to_scalar, z8); // 4+ 1
      const BlockType z12 = SHUFFLE8_0(to_scalar, z10); // 4+ 1
      sum = ADD8_0(sum, z11); 
      sum = ADD8_0(sum, z12); 
    }
    ans += ADDUPVECTOR(sum); 
  }

  BlockType sum ,
    *endx = x + blocks,
    *xx = x + bundles * x15,
    *yy = y + bundles * x15;
  sum = ZERO();
  for (; xx < endx; xx++, yy++) {
    BlockType  L1 = LOAD(xx),
    L2 = LOAD(yy);
    sum = ADD8(sum, VECTORADD(L1, L2));
  }
  ans += ADDUPVECTOR(sum);
  return ans;
}



 */




void trafo2Geno1Geno128(unit_t *OLD, Long snps, Long individuals,
			Long oldLDA, unit_t *ANS, Long ldAns) {
  // NOTE: Old and Ans could be the same destination

  assert(oldLDA == 2 * ldAns);
  
  const Long
    blocks = Blocks(snps),
    halfblocks = 1L + (blocks - 1L) / 2L;
   
  if (OLD == ANS) {
    ERR0("not programed yet");
    // (Braucht nur 1/4 Speicherplatz extra)
    // (i)drittes Viertel rauskopieren + 1-2 Spalten mehr vorne freilasen
    // (ii) haelfte durcharbeitne
    // (iii) drittel Vierel aus speicher + in die Spalte(n) davor im
    //       Speicher viertel Viertel spaltenweise reinkopieren
    // (iv) letztes Viertel aus dem Speicher mit den 1-2 Spalten frueher
    //      abarbeiten.
  }	
  
#if defined AVX512F
  const BlockType cross_idx =
    _mm512_set_epi32(TAKE_B + 14, TAKE_B + 12, TAKE_B + 10, TAKE_B + 8,
		     TAKE_B + 6, TAKE_B + 4, TAKE_B + 2, TAKE_B + 0,
		     14, 12, 10, 8, 6, 4, 2, 0);
#elif defined AVX2
  const BlockType cross_idx = _mm256_set_epi32(0, 0, 0, 0, 6, 4, 2, 0);
  const int halfreg_idx = 0 + (2 << 4);
#elif defined SSE2
  const int cross_idx = 0 + (2 << 2);
  const int halfreg_idx = 0 + (0 << 1);   
#endif
   
  const BlockType 
    E1N3 = SET32(0x88888888),
    N2E2N4=SET32(0x30303030),
    O1F1O2=SET32(0x0F000F00),  
    O2F2 = SET32(0x00FF00FF),
    O2F2O4=SET32(0x00FF0000),  
    O4F4 = SET32(0x0000FFFF);
#if ! defined PlainInteger32
  const BlockType O4F4O8 = SET64(0x0000FFFF00000000);
#endif	

  for (Long j=0; j<individuals; j++) {
    unit_t *Ans = ANS + j * ldAns,
      *halfmatrix = Ans + individuals * ldAns,
      *Old = OLD + j * oldLDA;
    
    for (Long i=0; i<halfblocks; i++) {	
      
      /*
	..    01010101 
	.. (& 01000100) >> 1, then "or" to previous
	Result:
	..  & 00110011 (clean up)
	.. (& 00110000) >> 2, then "or to previous
	 etc

	 
      */
      
      const BlockType a0 = LOAD((BlockType0*) Old);
      const BlockType b0 = LOAD((BlockType0*) (Old + UnitsPerBlock)); // OK
      const BlockType a1 = AND(a0, E1N1);
      const BlockType b1 = AND(b0, N1E1);
      const BlockType c1 = SHR32(AND(a0, N1E1), 1);
      const BlockType d1 = SHR32(AND(b0, E1N1), 1);
      AND4_C(2,1,E1N3);
      SHR32_4(3,2,1);
      OR4(4,3,1);
      
      AND4_C(5, 4, N2E2);
      AND4_C(6, 5, N2E2N4);
      SHR32_4(7, 6, 2);
      OR4 (8, 7, 5);
      
      AND4_C(9, 8, O1F1);
      AND4_C(10, 9, O1F1O2);
      SHR32_4(11, 10, 4);
      OR4 (12, 11, 9);
       
      AND4_C(13, 12, O2F2);
      AND4_C(14, 13, O2F2O4);
      SHR32_4(15, 14, 8);
      OR4 (16, 15, 13);
	 
#if ! defined PlainInteger32
       AND4_C(17, 16, O4F4);
       AND4_C(18, 17, O4F4O8);
       SHR64_4(19, 18, 16);
       OR4(20, 19, 17);
#endif	

#if defined AVX512F
#define TAKE_B 16
       const BlockType a22 = DOUBLE_CROSSLINE_PERMUTE(a20, cross_idx, b20);
       const BlockType c22 = DOUBLE_CROSSLINE_PERMUTE(c20, cross_idx, c20);
#elif defined AVX2 || defined SSE2
       PERMUTE4(21, 20, cross_idx);
       const BlockType a22 = HALFREGISTER_PERMUTE(a21, b21, halfreg_idx);
       const BlockType c22 = HALFREGISTER_PERMUTE(c21, d21, halfreg_idx);
#elif defined PlainInteger64
       BlockType a22, a22;
       ((uint32_t*) a22)[0] = ((uint32_t*) a20)[0];
       ((uint32_t*) a22)[1] = ((uint32_t*) b20)[0];
       ((uint32_t*) c22)[0] = ((uint32_t*) c20)[0];
       ((uint32_t*) c22)[1] = ((uint32_t*) d20)[0];
#elif defined PlainInteger32
       const BlockType a22;
       const BlockType c22; 
       ((uint_16_t*) a22)[0] = ((uint_16_t*) a20)[0];
       ((uint_16_t*) a22)[1] = ((uint_16_t*) b20)[0];
       ((uint_16_t*) c22)[0] = ((uint_16_t*) c20)[0];
       ((uint_16_t*) c22)[1] = ((uint_16_t*) d20)[0];
#else
#error ERROR.
#endif	
	  
       *((BlockType0*) Ans) = a22;
       *((BlockType0*) halfmatrix) = c22;
       Ans += UnitsPerBlock; // OK
       halfmatrix += UnitsPerBlock; // OK
       Old += 2 * UnitsPerBlock; // OK
     } // for i
   } // for j
}


FileBinary(128)
SCALAR012(128)


bool coding_2v128(unit_t *M, Long block16, Long rest, uint32_t* ans) {
  __m128i *pM = (__m128i*) M;
  BlockType two = SET32(2),
    wronga=ZERO(), 
    wrongb=ZERO(), 
    wrongc=ZERO(), 
    wrongd=ZERO();
  for (Long jj=0; jj<block16; jj++) {
     const BlockType
      a = LOAD(pM++),
      b = LOAD(pM++),
      c = LOAD(pM++),
      d = LOAD(pM++);
    wronga = OR(wronga, GREATER_THAN(a, two));
    wrongb = OR(wrongb, GREATER_THAN(b, two));
    wrongc = OR(wrongc, GREATER_THAN(c, two));
    wrongd = OR(wrongd, GREATER_THAN(d, two));
    const BlockType
      b1 = SHL32(b, BitsPerByte),
      c1 = SHL32(c, 2 * BitsPerByte),
      d1 = SHL32(d, 3 * BitsPerByte);
    const BlockType
      ab = OR(a, b1),
      cd = OR(c1, d1);
    const BlockType
      e = OR(ab, cd),
      f = SHR64(e, 30),
      g = OR(e, f),
      h = SHUFFLE32(g, 2),
      i = SHL32(h, 4),
      j = OR(i, g);
    ans[jj] = (uint32_t) _mm_cvtsi128_si32(j);
  }
 	
  wronga = OR(wronga, wrongb);
  wrongc = OR(wrongc, wrongd);
  UnionType wrongx; // OK
  wrongx.vi = OR(wronga, wrongc);
  bool wrong = wrongx.u64[0] | wrongx.u64[1];
  
  if (rest > 0) {
    rest *= BitsPerCode;
    _4Byte *mm = (_4Byte*) pM;
    _4Byte C = 0;
    for (int shft = 0; shft < rest; mm++, shft+=BitsPerCode) {
      wrong |= *mm > 2;
      C |= *mm << shft;
    }
    ans[block16] = (uint32_t) C;
  }
  return wrong;
}


void coding_2v128_4Byte(unit_t *M,Long ldM,
		       Long start_snp, Long end_snp, 
		       Long start_individual, Long end_individual,
		       int VARIABLE_IS_NOT_USED cores,
		       double VARIABLE_IS_NOT_USED *G,
		       unit_t *ans, Long VARIABLE_IS_NOT_USED individuals,
		       Long ldAns) {
  if (start_snp % CodesPerUnit != 0) BUG; // OK guaranted by calling function  
  Long
    cur_snps = end_snp - start_snp,
    codes_atonce = CodesPerByte * UnitsPerBlock, // OK
    block16 = cur_snps / codes_atonce,
    rest = cur_snps % codes_atonce;
  // assert(codes_atonce == 64);
  
  //  printf("cur_snps=%d %d %d atonce=%d %d\n", cur_snps, block16, rest,CodesPerByte, codes_atonce);

  unit_t
    *a= ans + start_snp / CodesPerUnit + start_individual * (Long) ldAns;// OK
  Long deltaIndiv = end_individual - start_individual;

  for (Long i=0; i<deltaIndiv; i++) {
    if (coding_2v128(M + i * ldM, block16, rest, (uint32_t*)(a + i * ldAns))) {
      ERR0("2bit#128: values outside 0,1,2.");
    }
  }
}


void coding_2v128(unit_t *M,  Long ldM,
		  Long start_snp, Long end_snp,
		  Long start_individual, Long end_individual,     
		  basic_options *opt, double *G, SEXP Ans) {
  Long *info =  GetInfo(Ans);
  ASSERT_LITTLE_ENDIAN;
  coding_2v128_4Byte(M, ldM,
		    start_snp, end_snp, 
		    start_individual,  end_individual,	
		    GreaterZero(opt->cores),
		    G, Align(Ans, opt), info[INDIVIDUALS], info[LDA]);
}


void allele_sum_2v128(unit_t *Code, Long snps, Long individuals, Long lda,
		      Long ldabitalign,
		      basic_options VARIABLE_IS_NOT_USED *opt,
		      Ulong *ans) {
  allele_sum_I_2bit(Code, snps, individuals, lda, ldabitalign, NULL, NULL, ans);
}


#else // !defined AVX
#include "avx_miss.h"
#include "MXinfo.h"


bool coding_2v128(_4Byte *M, Long block16, Long rest, uint32_t* ans) Su

DeclareVersion2(128,Su,Sv,Sm)
SIMD_MISS(2bit128, ssse3);


#endif  // defined SSSE3 || SSE2
