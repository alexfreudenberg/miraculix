 
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2019  Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#define MY_METHOD Packed

#include "IntrinsicsBase.h"

#include "error.h"
#include "MX.h"
#include "options.h"

#if defined SSE2

#if defined AVX2
#undef AVX2
#endif

#if defined AVX512
#undef AVX512
#endif

#if defined AVX512
#undef AVX512
#endif


/* 
__m128i _mm_mullo_epi16 (__m128i a, __m128i b) // sehr teuer

gleich teuer aber besser:
__m128i _mm_madd_epi16 (__m128i a, __m128i b)


*/


#include "packedDef.h"


// [1 + (X - 1) / UnitsPerBlock] * UnitsPerBlock + (UnitsPerBlock - 1)
// last part: to make finally sure that first used element of R-vector
// is aligned (maximum distance to 0-elemtn: UnitsPerBlock - 1)
// first part: then memory is read in blocks; so next number of Bytes that is
// is divisible by Unitsperblock * Bytesperunit is taken:
// 1 + (X - 1) / UnitsPerBlock]

Uint BytesPerBlockPacked() { return BytesPerBlock; }
Uint CodesPerBlockPacked() { return CodesPerBlock; }
//Uint UnitsPerIndivPacked(Uint snps) { return UnitsPerIndiv256(snps); }
#define UnitsPerIndivPacked UnitsPerIndiv256
#define UPI  UnitsPerIndivPacked




Uint INLINE ADDUPVECTOR(BlockType0 x) {
  //  BlockType zero;
  //  ZERO(zero);
  // BlockType vsum;
  // SAD8(vsum, x, zero); //
  //  Uint L1, L2;
  // EXTRACT16(L1, vsum, 0);
  //EXTRACT16(L2, vsum, 4);

  BlockUnitType *y = (BlockUnitType0*) &x;
  Uint ans = 0;
  for (Uint i=0; i<16; i++) ans += y->u8[i];
  return ans;

}


#include "packedIntern.h"
PackedInitIntern(128)

Ulong sumGeno128(Uint *S, Uint snps, Uint individuals) {
  return sumGenoIntern(S, snps, individuals);
}

	
void haplo2geno128(Uint *X, Uint snps, Uint individuals, Uint unitsPerIndiv,
		       Uint *Ans) {
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  haplo2genoIntern(X, snps, individuals, unitsPerIndiv, Ans);
}
 
 
SEXP matrix_start_packed( Uint snps, Uint individuals,SEXP file) {  
  return matrix_start_Intern( snps, individuals, file);
}


void crossprod_packed(Uint  *CGM, Uint snps, Uint individuals, double  *ans) {
  crossprodIntern(CGM, snps, individuals, ans);
}


SEXP allele_freq128(SEXP GM) {
  return allele_freqIntern(GM);
}
  
bool usePacked(snpcoding method) {
   option_type *global = &(KEYT()->global);
  return method == Packed ? true : global->genetics.efficient;  
}


bool use128() { return sse2; }


#else // !defined SSE2
void static SSE2missing() { ERR("'packed' needs the availablity of 'SSE2'"); }
#define Sm { SSE2missing(); return R_NilValue; }
#define Su { SSE2missing(); return 0; }
#define Sv { SSE2missing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
SEXP allele_freq128(SEXP V GM) Sm
void crossprod_packed(Uint V* CGM, Uint V snps, Uint V individuals,
			double V *ans) Sv
SEXP matrix_start_packed(Uint V snps, Uint V individuals, SEXP V G) Sm

Uint BytesPerBlockPacked() Su
Uint CodesPerBlockPacked() Su
Uint UnitsPerIndivPacked(Uint V snps) Su
void haplo2geno128(Uint *X, Uint snps, Uint individuals, Uint unitsPerIndiv,
		       Uint *Ans) Sv
Ulong sumGeno128(Uint V *S, Uint V snps, Uint V individuals) Su

//#include "miraculix.h"
bool usePacked(snpcoding method) {
  option_type *global = &(KEYT()->global);
  if (global->genetics.efficient || method != Packed) return false;
  ERR("'SSE2' is needed in 'Packed'. Set 'RFoptions(efficient=TRUE)'.");
}

bool use128() { return false; }

#endif  // defined SSE2 || SSE2
