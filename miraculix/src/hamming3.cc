

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


Copyright (C) 2014 -- 2019 Martin Schlather
Copyright (C) 2014 -- 2015 Florian Skene: SSE2+SSSE3

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

#define BitsPerCode 4L
#define atonce 4L

#include "IntrinsicsBase.h"

#define MY_METHOD Hamming3

#if defined SSSE3

#if !defined SSE2
#define SSE2 1
#endif

#if defined AVX2
#undef AVX2
#endif

#if defined AVX512
#undef AVX512
#endif


#include "miraculix.h"
#include <General_utils.h>
#include "haplogeno.h"
#include "hamming.intern.h"
#include "xport_import.h"

static uint64_t code1[] = {0, 6, 15},    //0, 1, 2
  code2[] = {0, 3, 15};    //0, 1, 2
static Uint rev_coding[16] ={0, 0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 2};


Uint UnitsPerIndivH3 UnitsPerIndivH23

SEXP matrix_startH3 matrix_startH23
void codingH3 codingH23
void haplo2genoH3 h2gH23
SEXP get_matrixH3 get_H23
Ulong sumGenoH3 sumGeno23
SEXP allele_freqH3 allele_freq23
SEXP get_matrixN_H3 get_matrixN_H23



void crossprod_H3(Uint *M, Uint snps, Uint individuals, double *ergb) {
  Uint
    unitsPerIndiv = UnitsPerIndivH3(snps),
    *M_T =  M + unitsPerIndiv * individuals;

  Uint
    blocks = unitsPerIndiv / UnitsPerBlock,
    loops = blocks / atonce;
  
  const unsigned _LUT[] = {0x02010100, 0x03020201, 0x03020201, 0x04030302};
  BlockType LUT, mask, zero;


  LOADU(LUT, (BlockType0 *) _LUT); // _mm_set_epi32(,,,)
  SET8(mask, 0x0F);
  ZERO(zero);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 20)
#endif
  for (Ulong i = 0; i<individuals; i++) {
    BlockType
      *mpp = (BlockType0 *) (M + i * unitsPerIndiv);
     for (Ulong j = i; j<individuals; j++) {
     uint64_t tot = 0;
      BlockType *mptr = mpp,
	*m_tpp = (BlockType0*) (M_T + j * unitsPerIndiv);
      for (Uint k = 0; k < loops; k++) {    // '/4' hinzugefuegt
	uni128  intermed;
	BlockType v0, lo, count0, count1, L1, L2;
#define COUNT(count0)						\
	LOAD(L1, mptr);mptr++;					\
	LOAD(L2, m_tpp);m_tpp++;				\
	AND(v0, L1, L2);					\
	AND(lo, mask, v0);					\
	SHR16(L1, v0, 4);					\
	AND(v0, mask, L1); /* hi */				\
	SHUFFLE8(lo, LUT, lo);					\
	SHUFFLE8(v0, LUT, v0);					\
	ADD8(count0, lo, v0);
	
	COUNT(count0);
	COUNT(v0);
	ADD8(count0, count0, v0);
	COUNT(count1);
	COUNT(v0);
	ADD8(count1, count1, v0);
	ADD8(count0, count0, count1);
	
	//bitwise difference between 'tota' and all zeros
	SAD8(intermed VI, count0, zero);
	// adds the two bitcounts which are saved as unsigned 64bit integers
	// to continuous bitcount
	tot += (Uint) intermed.u64[0] + (Uint)intermed.u64[1];
      }
      ergb[j + i * individuals] = ergb[i + j * individuals] = (double) tot;
    }
  }
}


#else 
#include "error.h"
#include "MX.h"
void static SSSE3missing() {
  ERR("'Hamming' needs at least the availablity of 'SSSE32'");
}
#define Sm { SSSE3missing(); return R_NilValue; }
#define Su { SSSE3missing(); return 0; }
#define Sv { SSSE3missing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
Uint UnitsPerIndivH3(Uint snps) Su
SEXP matrix_startH3( Uint V snps, Uint V individuals, SEXP V file)Sm
void codingH3(Uint V *M, Uint V start_individual, Uint V end_individual, 
	     Uint V start_snp, Uint V end_snp, Uint V Mnrow, SEXP V Ans,
	     double V *G)Sv
void haplo2genoH3(Uint V *X, Uint V snps, Uint V individuals,
		    Uint V Mnrow, Uint V *code) Sv
void crossprod_H3(Uint V *M, Uint V snps, Uint V individuals, double V *ergb){
  ERR("Hamming3 snp coding needs at least SSSE3");
}
SEXP get_matrixH3(SEXP V SNPxIndiv) Sm
Ulong sumGenoH3(Uint V *M, Uint V snps, Uint V indiv) Su
SEXP allele_freqH3(SEXP V SNPxIndiv) Sm
SEXP get_matrixN_H3(SEXP V SNPxIndiv, SEXP V Snps) Sm

#endif // s

