
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


#ifndef miraculix_internH_H
#define miraculix_internH_H 1

void codingH(Uint *X, Uint start_individual, Uint end_individual, 
	     Uint start_snp, Uint end_snp, Uint Mnrow,
	     SEXP Code, uint64_t *code1, uint64_t *code2, snpcoding method);
void haplo2genoH(Uint *X, Uint snps, Uint individuals, Uint Mnrow, Uint *code,
		 uint64_t *code1, uint64_t *code2, snpcoding method);
SEXP get_matrixH(SEXP SNPxIndiv, Uint *rev_coding);
SEXP allele_freqH(SEXP SNPxIndiv, Uint *rev_coding);		
Ulong sumGenoH(Uint *M, Uint snps, Uint individuals, snpcoding method,
	       Uint *rev_coding);
SEXP get_matrixN_H(SEXP SNPxIndiv, SEXP Snps, Uint *rev_coding);


#define  MatrixCoding(start_snp, cur_snps, start_individual, cur_indiv, FROM) \
  /*  Ulong unitsPerIndiv = GetUPI(snps, MY_METHOD),	*/		\
  /*  memPerMatrix = unitsPerIndiv * individuals; */			\
  Uint  mini = sizeof(BlockType0) / sizeof(uint64_t),			\
    CodesPerMiniblock = CodesPerBlock / mini;				\
  Uint									\
  *M = code + start_snp / CodesPerUnit + start_individual * unitsPerIndiv, \
    *M_T = M + unitsPerIndiv * individuals;				\
  for (Ulong i = 0; i<cur_indiv; i++) {					\
    Uint count=0,							\
      *pX = X + i * Mnrow;						\
    uint64_t *Mptr = (uint64_t*) (M + i * unitsPerIndiv),		\
      *M_Tptr =  (uint64_t*) (M_T + i * unitsPerIndiv);			\
    for (Uint s = 0; s<cur_snps; s++) {					\
      Uint snp = FROM;							\
      *Mptr = *Mptr << BitsPerCode;					\
      *Mptr |= code1[snp];						\
      *M_Tptr = *M_Tptr << BitsPerCode;					\
      *M_Tptr |= code2[snp];						\
      if (++count == CodesPerMiniblock) {				\
	count = 0;							\
	Mptr++;								\
	M_Tptr++;							\
      }									\
    }									\
    if (count > 0) {							\
      Uint shift = BitsPerCode * (CodesPerMiniblock - count);		\
      *Mptr = *Mptr << shift;						\
      *M_Tptr = *M_Tptr << shift;					\
    }									\
  }


#define DOGET(DO)							\
  Ulong mini = sizeof(BlockType0) / sizeof(uint64_t),			\
    BitsPerMiniblock = BitsPerBlock / mini,				\
    CodesPerMiniblock = CodesPerBlock / mini;				\
  for (Uint i=0; i<individuals; i++) {					\
    uint64_t *Ma = (uint64_t*) (M + i * unitsPerIndiv);			\
    for (Uint s=0; s<snps; s++) {					\
      DO rev_coding[(Ma[s / CodesPerMiniblock]				\
		   >> (BitsPerMiniblock - BitsPerCode -			\
		       BitsPerCode  * (s % CodesPerMiniblock))) & CodeMask]; \
    }									\
  }


#define DONGET(DO)							\
  Ulong mini = sizeof(BlockType0) / sizeof(uint64_t),			\
    BitsPerMiniblock = BitsPerBlock / mini,				\
    CodesPerMiniblock = CodesPerBlock / mini;				\
  for (Uint i=0; i<individuals; i++) {					\
    uint64_t *Ma = (uint64_t*) (M + i * unitsPerIndiv);			\
    for (Uint s=0; s<len; s++) {					\
      Rint S = which[s];						\
      DO rev_coding[(Ma[S / CodesPerMiniblock]				\
		     >> (BitsPerMiniblock - BitsPerCode -		\
			 BitsPerCode  * (S % CodesPerMiniblock))) & CodeMask]; \
    }									\
  }


#define UnitsPerIndivH23				\
  (Uint snps) {						\
    Uint blocks = (1L + (snps - 1L) / (atonce * CodesPerBlock)) * atonce; \
    return UnitsPerBlock * blocks;					\
  }


#define matrix_startH23							\
  (Uint snps,Uint individuals,  SEXP file) {				\
    SEXP Code;								\
    assert(BytesPerBlock == sizeof(BlockType0));			\
    PROTECT(Code = CreateEmptyCodeVector(snps, individuals, MY_METHOD)); \
    start_info(Code, file);		\
    UNPROTECT(1);							\
    return(Code);							\
}

#define codingH23							\
  (Uint *X, Uint start_individual, Uint end_individual,			\
   Uint start_snp, Uint end_snp, Uint Mnrow,				\
   SEXP Code,  double VARIABLE_IS_NOT_USED *G) { /* from file_get */	\
    codingH(X, start_individual, end_individual, start_snp, end_snp,	\
	    Mnrow, Code, code1, code2, MY_METHOD);	\
  }

#define h2gH23\
  (Uint *X, Uint snps, Uint individuals, Uint Mnrow, Uint *code){	\
  haplo2genoH(X, snps, individuals, Mnrow, code, code1, code2, MY_METHOD);\
}				

#define sumGeno23\
  (Uint *M, Uint snps, Uint individuals) {	\
    return sumGenoH(M, snps, individuals, MY_METHOD, rev_coding);	\
  }

#define allele_freq23					\
  (SEXP SNPxIndiv) {					\
    return allele_freqH(SNPxIndiv, rev_coding);	\
}

#define get_codingH23						\
  (SEXP SNPxIndiv) {						\
    return get_codingH(SNPxIndiv, MY_METHOD);			\
  }


#define allele_freqH23						\
  (SEXP SNPxIndiv) {						\
    return allele_freqH(SNPxIndiv, MY_METHOD);			\
  }

#define get_H23							\
  (SEXP SNPxIndiv) { return get_matrixH(SNPxIndiv, rev_coding); }

#define  get_matrixN_H23				\
  (SEXP SNPxIndiv, SEXP Snps) {				\
    return get_matrixN_H(SNPxIndiv, Snps, rev_coding);	\
  }


#endif
