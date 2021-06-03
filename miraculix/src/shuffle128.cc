
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

#define MY_METHOD Shuffle

#include "IntrinsicsBase.h"

#include "error.h"
#include "MX.h"
#include "options.h"

#if defined SSSE3 

#if defined AVX2
#undef AVX2
#endif

#if defined AVX512
#undef AVX512
#endif


#include "shuffleDef.h"


// [1 + (X - 1) / UnitsPerBlock] * UnitsPerBlock + (UnitsPerBlock - 1)
// last part: to make finally sure that first used element of R-vector
// is aligned (maximum distance to 0-elemtn: UnitsPerBlock - 1)
// first part: then memory is read in blocks; so next number of Bytes that is
// is divisible by Unitsperblock * Bytesperunit is taken:
// 1 + (X - 1) / UnitsPerBlock]

Uint BytesPerBlockShuffle() { return BytesPerBlock; }
Uint CodesPerBlockShuffle() { return CodesPerBlock; }
#define UPI UnitsPerIndiv256

Uint INLINE ADDUPVECTOR(BlockType0 x) {
  BlockType zero;
  ZERO(zero);
  BlockType vsum;
  SAD8(vsum, x, zero); //
  Uint L1, L2;
  EXTRACT16(L1, vsum, 0);
  EXTRACT16(L2, vsum, 4);
  return L1 + L2;
}


#include "shuffleIntern.h"
ShuffleInitIntern(128)


void crossprod_shuffle(Uint *CGM, Uint snps, Uint individuals, double  *ans) {
  crossprodIntern(CGM, snps, individuals, ans);
}
SEXP matrix_start_shuffle( Uint snps,Uint individuals, SEXP file) {  
  return matrix_start_Intern(snps, individuals, file);
}

bool useShuffle(snpcoding method) {
  option_type *global = &(KEYT()->global);
  return method == Shuffle ? true : global->genetics.efficient;  
}


#else // !defined SSSE3
void static SSSE3missing() {ERR("'shuffle' needs the availablity of 'SSSE3'"); }
#define Sm { SSSE3missing(); return R_NilValue; }
#define Su { SSSE3missing(); return 0; }
#define Sv { SSSE3missing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
Uint BytesPerBlockShuffle() Su
Uint CodesPerBlockShuffle() Su
void crossprod_shuffle(Uint V* CGM, Uint V snps, Uint V individuals,
			double V *ans) Sv
SEXP matrix_start_shuffle(Uint V snps, Uint V individuals, SEXP V G) Sm
 

bool useShuffle(snpcoding method) {
  option_type *global = &(KEYT()->global);
  if (global->genetics.efficient || method != Shuffle) return false;
  ERR("'SSSE3' is needed for vector multiplication in 'Shuffle'. Set 'RFoptions(efficient=TRUE)'.");
}
void ShuffleInit128(){ return; }

#endif  // defined SSSE3 || SSE2
