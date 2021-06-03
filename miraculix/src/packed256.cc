
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

#define MY_METHOD Packed256

#include "IntrinsicsBase.h"
#include "error.h"
#include "MX.h"
#include "options.h"

#if defined AVX2

#include "packedDef.h"


// [1 + (X - 1) / UnitsPerBlock] * UnitsPerBlock + (UnitsPerBlock - 1)
// last part: to make finally sure that first used element of R-vector
// is aligned (maximum distance to 0-elemtn: UnitsPerBlock - 1)
// first part: then memory is read in blocks; so next number of Bytes that is
// is divisible by Unitsperblock * Bytesperunit is taken:
// 1 + (X - 1) / UnitsPerBlock]

Uint BytesPerBlockPacked256() { return BytesPerBlock; }
Uint CodesPerBlockPacked256() { return CodesPerBlock; }
// Uint UnitsPerIndivPacked256(Uint snps)  { return UnitsPerIndiv256(snps); }
#define UnitsPerIndivPacked256 UnitsPerIndiv256
#define UPI  UnitsPerIndivPacked256



Uint INLINE ADDUPVECTOR(BlockType0 x) {
  BlockType zero;
  ZERO(zero);
  BlockUnitType vsum;
  SAD8(vsum VI, x, zero);  
  return vsum.u32[0]+ vsum.u32[2] + vsum.u32[4] + vsum.u32[6];
}


#include "packedIntern.h"
PackedInitIntern(256)


Ulong sumGeno256(Uint *S, Uint snps, Uint individuals) {
  return sumGenoIntern(S, snps, individuals);
}

void haplo2geno256(Uint *X, Uint snps, Uint individuals,
			  Uint unitsPerIndiv,Uint *Ans) {
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  haplo2genoIntern(X, snps, individuals, unitsPerIndiv, Ans);
}

     
SEXP matrix_start_packed256( Uint snps,Uint individuals, SEXP file) {  
  return matrix_start_Intern(snps, individuals, file);
}
 

void crossprod_packed256(Uint *CGM, Uint snps, Uint individuals, double *ans){
  crossprodIntern(CGM, snps, individuals, ans);
}


      
SEXP allele_freq256(SEXP GM) {
  return allele_freqIntern(GM);
}

bool usePacked256(snpcoding method) {
  option_type *global = &(KEYT()->global);
  return method == Packed256 ? true : global->genetics.efficient;  
}

bool use256() { return avx2; }


#else // !defined AVX2
void static AVX2missing() { ERR("'packed256' needs the availablity of 'AVX2'"); }
#define Sm { AVX2missing(); return R_NilValue; }
#define Su { AVX2missing(); return 0; }
#define Sv { AVX2missing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
SEXP allele_freq256(SEXP V GM) Sm
void crossprod_packed256(Uint V* CGM, Uint V snps, Uint V individuals,
			double V *ans) Sv
SEXP matrix_start_packed256(Uint V snps, Uint V individuals, SEXP V G) Sm

Uint BytesPerBlockPacked256() Su
Uint CodesPerBlockPacked256() Su
void haplo2geno256(Uint V * SNPxIndiv, Uint V snps, Uint V individuals,
			  Uint upi, Uint V *A) Sv
Ulong sumGeno256(Uint V *S, Uint V snps, Uint V individuals) Su

#include "miraculix.h"

bool usePacked256(snpcoding method) {
  option_type *global = &(KEYT()->global);
  if (global->genetics.efficient || method != Packed256) return false;
  ERR("'Packed256' is not available. Set 'RFoptions(efficient=TRUE)'.");
}


bool use256() { return false; }
void PackedInit256(){ return; }

#endif  // defined AVX2

