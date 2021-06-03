
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


#ifndef miraculix_Bit23intern_H
#define miraculix_Bit23intern_H 1

#include <stdio.h>
#ifdef DO_PARALLEL
#include <omp.h>
#endif

#include "intrinsics.h"
//#include "haplogeno.h"
#include "options.h"
#include "error.h"
#include "MX.h"
#include "Haplo.h"

#define BytesPerMiniblock 2L
#define BitsPerMiniblock (BitsPerByte * BytesPerMiniblock)
#define CodesPerMiniblock (BitsPerMiniblock / BitsPerCode)
#define genuineBitsPerMiniblock (CodesPerMiniblock * BitsPerCode)
#define two_genuineBitsPerMiniblock (1L << genuineBitsPerMiniblock)
#define deltaMiniblockBits (BitsPerMiniblock - genuineBitsPerMiniblock)


//  MiniblocksPerBlock guaranteed to be 4
#define MiniblocksPerBlock (BytesPerBlock / BytesPerMiniblock) 
#define genuineCodesPerBlock (MiniblocksPerBlock * CodesPerMiniblock)
#define genuineCodesPerUnit (genuineCodesPerBlock / UnitsPerBlock)
typedef uint16_t blockarray[MiniblocksPerBlock];

 
#define nr_genotypes 3

union block_compressed {
  blockarray b;
  uint32_t u32[2];
  BlockType x;
};

//
typedef  char table_type; //  10.585 
//typedef unsigned char table_type; //  10.612
//typedef unsigned int table_type; //   10.325 
//typedef int table_type; //  10.217 
void initiate_tableI(table_type **TABLE, Uint TABLESIZE,
		     Uint codesPerMiniblock, Uint bits, BlockType0 *result_code,
		     Uint *result_value, Uint NrResults);

SEXP get_matrix23_start(Uint snps, Uint individuals, SEXP G);


void printbits(BlockType0 x, Uint size, Uint bits);

SEXP matrix_start2(Uint snps, Uint individuals, SEXP file);

SEXP matrix_start3(Uint snps, Uint individuals,  SEXP file);


/* matrix_start */
#define START23(Init23) (Uint snps, Uint individuals, SEXP file) {	\
  assert(TWO_GENUINEBITSPERMINIBLOCK == two_genuineBitsPerMiniblock);	\
  assert(genuineCodesPerBlock == CodesPerBlock);			\
  SEXP Code;								\
  PROTECT(Code =  CreateEmptyCodeVector(snps, individuals, MY_METHOD));	\
  start_info(Code, file);				\
  UNPROTECT(1);								\
  return(Code);								\
  }


#define Coding23(start_individual, end_individual, start_snp, end_snp, Mnrow, FROM) \
  for (Ulong i=start_individual; i<end_individual; i++) {		\
    Uint *pX = X + (i - start_individual) * Mnrow,			\
      shift = 0,							\
      counter = 0;							\
    BlockType compressed = 0,						\
      *pAns = (BlockType0 *) (ans + i * unitsPerIndiv);			\
    for (Uint s=start_snp; s<end_snp; s++) {				\
      compressed |= geno_code[FROM] << shift;				\
      shift += BitsPerCode;						\
      if (++counter >= CodesPerMiniblock) {				\
	shift += deltaMiniblockBits; /* 3bit */				\
	counter = 0;							\
      }									\
      if (shift >= BitsPerBlock) {					\
	*pAns = compressed;						\
	pAns++;								\
	shift = counter = 0;						\
	compressed = 0;							\
      }									\
    }									\
    if (shift > 0) {							\
      *pAns = compressed;						\
      pAns++;								\
    }									\
  }


#define GET_MATRIX_N_23(GET, UPI)					\
  Uint									\
  *info = GetInfo(SNPxIndiv),						\
    individuals = info[INDIVIDUALS],					\
    len = length(Snps),							\
    *which = (Uint*) INTEGER(Snps),					\
    snps = info[SNPS],							\
    unitsPerIndiv = UPI(snps);					\
  SEXP Ans;								\
  PROTECT(Ans = get_matrix23_start(len, individuals, R_NilValue));	\
  Uint *ans = (Uint *) INTEGER(Ans);					\
  for (Uint i=0; i<individuals; i++) {					\
    BlockType *Ma = (BlockType0*) (M + i * unitsPerIndiv);		\
    for (Uint s=0; s<len; s++) {					\
      Rint S = which[s];						\
      *(ans++) = GET(Ma, S);						\
    }									\
  }									\
  UNPROTECT(1);								\
  return Ans

#define ALLELE_FREQ23(GET, UPI)						\
  Uint									\
    *info = GetInfo(SNPxIndiv),						\
    individuals = info[INDIVIDUALS],					\
    snps = info[SNPS],							\
    unitsPerIndiv = UPI(snps);					\
  SEXP Ans;								\
  PROTECT(Ans=allocVector(REALSXP, snps));				\
  double *ans = REAL(Ans);						\
  for (Uint i=0; i<snps; ans[i++] = 0.0);				\
  for (Uint i=0; i<individuals; i++) {					\
    BlockType *Ma = (BlockType0*) (M + i * unitsPerIndiv);		\
    for (Uint s=0; s<snps; s++) ans[s] += (double) GET(Ma, s);		\
  }									\
  double factor = 0.5 / (double) individuals;				\
  for (Uint i=0; i<snps; i++) ans[i] *= factor;				\
  UNPROTECT(1);								\
  return Ans




#define ZERONTH(UPI)							\
  Uint									\
    *info = GetInfo(SNPxIndiv),						\
    individuals = info[INDIVIDUALS],					\
    snps = info[SNPS],							\
    len = length(Snps),							\
    *which = (Uint*) INTEGER(Snps),					\
    unitsPerIndiv = UPI(snps);					\
  for (Uint i=0; i<individuals; i++) {					\
    BlockType *Ma = (BlockType0*) (M + i * unitsPerIndiv);		\
    for (Uint s=0; s<len; s++) {					\
      Rint S = which[s];						\
      Uint  j = S / genuineCodesPerBlock,				\
	idx = S % genuineCodesPerBlock;					\
      block_compressed x;						\
      x.x = Ma[j];							\
      x.b[idx / CodesPerMiniblock] &=					\
	~(CodeMask << ( (idx % CodesPerMiniblock) * BitsPerCode));	\
      Ma[j] = x.x;							\
    }									\
  }								       
  


#endif
