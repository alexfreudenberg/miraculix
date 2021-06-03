

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
#define repet 10L
#define atonce (3L * repet)

#include "IntrinsicsBase.h"
#include "xport_import.h"

#define MY_METHOD Hamming2


#if defined SSE2

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

static uint64_t code1[] = {0, 5, 10},    //0, 1, 2
  code2[] = {0, 3, 15};    //0, 1, 2
static Uint rev_coding[16] ={0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 2, 0,  0, 0, 0, 0};


Uint UnitsPerIndivH2 UnitsPerIndivH23

SEXP matrix_startH2 matrix_startH23
void codingH2 codingH23
void haplo2genoH2 h2gH23
SEXP get_matrixH2 get_H23
Ulong sumGenoH2 sumGeno23
SEXP allele_freqH2 allele_freq23
SEXP get_matrixN_H2 get_matrixN_H23



void crossprod_H2(Uint *M, Uint snps, Uint individuals, double *ergb) {
  // cf line 817 onwards in plink_calc.c of plink-1.07-src.zip  on http://zzz.bwh.harvard.edu/plink/download.shtml

  
  Uint
    unitsPerIndiv = UnitsPerIndivH2(snps),
    *M_T =  M + unitsPerIndiv * individuals;
  
  // M*M^T:
  // outer for-loop represents rows of M (only within partition, if > 1 cores
  // inner for-loop represents columns of M^T starting at diagonal element
  Uint
    blocks = unitsPerIndiv / UnitsPerBlock,
    loops = blocks / atonce;
  BlockType m1, m2, m4, m8;
  uint64_t OOO1 = UINT64_C(0x0001000100010001);

  SET8(m1, 0x55);
  SET8(m2, 0x33);
  SET8(m4, 0x0f);
  SET16(m8,0x00ff);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 20)
#endif 
  for (Ulong i = 0; i<individuals; i++) {
    BlockType *mpp = (BlockType0 *) (M + i * unitsPerIndiv);
    for (Ulong j = i; j<individuals; j++) {
      BlockType
	*m_tpp = (BlockType0*) (M_T + j * unitsPerIndiv),
	*mptr = mpp; //mptr set back to first element of row in 	
      uint64_t tot = 0;
      for (Uint s=0; s<loops; s++) {
	uni128 acc;
	ZERO(acc VI);
	// do-while runs 10x before continuing
	for (Uint t=0; t<repet; t++) {
	  BlockType L1, L2, half1, half2;
	  LOAD(L1, mptr); mptr++;
	  LOAD(L2, m_tpp); m_tpp++;
	  AND(half1, L1, L2);
	  SHR64(L1, half1, 1);	  
	  AND(half2, L1, m1);
	  AND(half1, half1, m1);
	  OR(half1, half1, half2);
	  LOAD(L1, mptr); mptr++;
	  LOAD(L2, m_tpp); m_tpp++;
	    BlockType count1, count2;
	  AND(count1, L1, L2);
	  ADD64(count1, count1, half1);
	  LOAD(L1, mptr); mptr++;
	  LOAD(L2, m_tpp);m_tpp++;
	  AND(count2, L1, L2);
	  ADD64(count2, count2, half2);
#define dummy1 half1
#define dummy2 half2
	  AND(dummy1, count1, m2);
	  SHR64(L1, count1, 2);
	  AND(dummy2, L1, m2);
	  ADD64(count1,dummy1, dummy2);
	  AND(dummy1, count2, m2);
	  SHR64(L1, count2, 2);
	  AND(dummy2, L1, m2);
	  ADD64(dummy1, dummy1, dummy2);
	  ADD64(count1, count1, dummy1);
	  AND(dummy1, count1, m4);
	  SHR64(L1, count1, 4);
	  AND(dummy2, L1, m4);
	  ADD64(dummy1, dummy1, dummy2);
	  ADD64(acc VI, acc VI, dummy1);
	}
	BlockType L1, L2;
	AND(L1, acc VI, m8);
	SHR64(L2, acc VI, 8);
	AND(L2, L2, m8);
	ADD64(acc VI, L1, L2);
	tot += ((acc.u64[0] + acc.u64[1]) * OOO1) >> 48;
      }
      ergb[j + i * individuals] = ergb[i + j * individuals] = (double) tot;
    }
  }
}



#else 
#include "error.h"
#include "MX.h"
void static SSEmissing() {
  ERR("'Hamming' needs at least the availablity of 'SSE2'");
}
#define Sm { SSEmissing(); return R_NilValue; }
#define Su { SSEmissing(); return 0; }
#define Sv { SSEmissing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
Uint UnitsPerIndivH2(Uint snps) Su
SEXP matrix_startH2( Uint V snps, Uint V individuals,SEXP V file)Sm
void codingH2(Uint V *M, Uint V start_individual, Uint V end_individual, 
	     Uint V start_snp, Uint V end_snp, Uint V Mnrow, SEXP V Ans,
	     double V *G)Sv
void haplo2genoH2(Uint V *X, Uint V snps, Uint V individuals,
		    Uint V Mnrow, Uint V *code) Sv
void crossprod_H2(Uint V *M, Uint V snps, Uint V individuals, double V *ergb){
  ERR("Hamming snp coding needs at least SSE2");
}
SEXP get_matrixH2(SEXP V SNPxIndiv) Sm
Ulong sumGenoH2(Uint V *M, Uint V snps, Uint V indiv) Su
SEXP allele_freqH2(SEXP V SNPxIndiv) Sm
SEXP get_matrixN_H2(SEXP V SNPxIndiv, SEXP V Snps) Sm

#endif // s

