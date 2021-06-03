
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

#include "IntrinsicsBase.h"

#if defined SSSE3 || defined SSE2

#if defined SSSE3 && !defined SSE2
#define SSE2 1
#endif

#if defined AVX2
#undef AVX2
#endif

#if defined AVX512
#undef AVX512
#endif



#define BitsPerCode 4L
#include <General_utils.h>
#include "haplogeno.h"
#include "Haplo.h"
#include "align.h"
#include "hamming.intern.h"

Uint CodesPerBlockH() { return CodesPerBlock; }
Uint BitsPerCodeH() { return BitsPerCode; }
Uint BytesPerBlockH() { return BytesPerBlock; }
Uint *AlignH(SEXP Code, Uint nr, bool test) { return AlignTest(Code, nr, test);}

void codingH(Uint *X, Uint start_individual, Uint end_individual,		
	     Uint start_snp, Uint end_snp, Uint Mnrow,			
	     SEXP Code, uint64_t *code1, uint64_t *code2,
	     snpcoding method) {
  if (start_snp % CodesPerBlock != 0) BUG;				
  Uint
    *info = GetInfo(Code),						
    individuals = info[INDIVIDUALS],					
    snps = info[SNPS],						
    unitsPerIndiv = GetUPI(snps, method),
    cur_snps = end_snp - start_snp,					
    cur_indiv = end_individual-start_individual,			
    *code = Align(Code, ALIGN_SSE);				
  MatrixCoding(start_snp, cur_snps, start_individual, cur_indiv, FROMINPUT); 
}

void haplo2genoH(Uint *X, Uint snps, Uint individuals, Uint Mnrow, Uint *code,
		 uint64_t *code1, uint64_t *code2, snpcoding method) {
  Uint unitsPerIndiv = GetUPI(snps, method);
  MatrixCoding(0, snps, 0, individuals, FROMHAPLO);
}


#define ANSPLUS *(ans++) =
SEXP get_matrixH(SEXP SNPxIndiv, Uint *rev_coding) {
  Uint *info = GetInfo(SNPxIndiv),				
    individuals = info[INDIVIDUALS],				
    snps = info[SNPS],
    method = info[METHOD],
    unitsPerIndiv = GetUPI(snps, (snpcoding) method),
    *M = Align(SNPxIndiv, ALIGN_SSE);	

  SEXP Ans;							
  PROTECT(Ans=allocMatrix(INTSXP, snps, individuals));	
  Uint *ans = (Uint *) INTEGER(Ans);				
  DOGET(ANSPLUS);						
  UNPROTECT(1);						
  return Ans;							
}


#define SUMUP sum +=
Ulong sumGenoH(Uint *M, Uint snps, Uint individuals, snpcoding method,
	       Uint *rev_coding) {
  Ulong
    unitsPerIndiv = GetUPI(snps, method),
    sum = 0L;					 
  DOGET(SUMUP);						
  return sum;							
}


#define FREQ ans[s] += (double)
SEXP allele_freqH(SEXP SNPxIndiv, Uint *rev_coding) {		
  Uint *info = GetInfo(SNPxIndiv),				
    individuals = info[INDIVIDUALS],				
    snps = info[SNPS],						
    method = info[METHOD],					
    unitsPerIndiv = GetUPI(snps, (snpcoding) method),
    *M = Align(SNPxIndiv, ALIGN_SSE);	
  SEXP Ans;							
   PROTECT(Ans=allocVector(REALSXP, snps));			
  double *ans = REAL(Ans);					
  for (Uint i=0; i<snps; ans[i++] = 0.0);			
  DOGET(FREQ);							
  double factor = 0.5 / (double) individuals;			
  for (Uint i=0; i<snps; i++) ans[i] *= factor;			
  UNPROTECT(1);							
  return Ans;							
}


SEXP get_matrixN_H(SEXP SNPxIndiv, SEXP Snps, Uint *rev_coding) {
  Uint *info = GetInfo(SNPxIndiv),				
    individuals = info[INDIVIDUALS],				
    len = length(Snps),						
    *which = (Uint*) INTEGER(Snps),				
    snps = info[SNPS],						
    method = info[METHOD],
    unitsPerIndiv = GetUPI(snps, (snpcoding) method),
    *M = Align(SNPxIndiv, ALIGN_SSE);
  SEXP Ans;							
  PROTECT(Ans=allocMatrix(INTSXP, len, individuals));		
  Uint *ans = (Uint *) INTEGER(Ans);				
  DONGET(ANSPLUS);						
  UNPROTECT(1);							
  return Ans;							
}



SEXP zeroNthGenoH(SEXP SNPxIndiv, SEXP Snps, snpcoding method) {		
  Uint *info = GetInfo(SNPxIndiv),					
    individuals = info[INDIVIDUALS],					
    len = length(Snps),							
    *which = (Uint*) INTEGER(Snps),					
    snps = info[SNPS],
    unitsPerIndiv = GetUPI(snps, method),
    *M = Align(SNPxIndiv, ALIGN_SSE);		
  Ulong									
    mini = sizeof(BlockType0) / sizeof(uint64_t),			
    BitsPerMiniblock = BitsPerBlock / mini,				
    CodesPerMiniblock = CodesPerBlock / mini;				
  for (Uint i=0; i<individuals; i++) {		
    uint64_t *Ma = (uint64_t*) (M + i * unitsPerIndiv);   
    for (Uint s=0; s<len; s++) {					
      Rint S = which[s];						
      Ma[S / CodesPerMiniblock]	&=					
	~ (CodeMask << (BitsPerMiniblock - BitsPerCode - 		
			BitsPerCode  * (S % CodesPerMiniblock)));	
    }									
  }									
  return SNPxIndiv;							
}




#else // no sse2
#include "error.h"
#include "MX.h"
void static SSEmissing() { ERR("'Hamming' needs the availablity of 'SSE2'");}
#define Sm { SSEmissing(); return R_NilValue; }
#define Su { SSEmissing(); return 0; }
#define Sv { SSEmissing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif


Uint CodesPerBlockH() Su
Uint BitsPerCodeH() Su
Uint BytesPerBlockH() Su
SEXP zeroNthGenoH(SEXP V SNPxIndiv, SEXP V Snps, snpcoding V method) Sm
Uint *AlignH(SEXP Code, Uint nr, bool test) Su


#endif // sse2

