
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2019 -- 2019  Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, writne to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#ifndef miraculix_shuffle_H
#define miraculix_shuffle_H 1

#include "MX.h"


void matrixshuffle_mult(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_shuffle(Uint individuals, Uint snps, SEXP G);
void matrix_shuffle(Uint *M, Uint start_individual, Uint end_individual, 
		    Uint start_snp, Uint end_snp, Uint Mnrow,
		    SEXP Ans, double * G);

SEXP get_matrixshuffle(SEXP SNPxIndiv);
SEXP matrix_coding_shuffle(Uint *M, Uint snps, Uint individuals);

SEXP create_codevector(Uint snps, Uint individuals);
Uint *AlignShuffle(SEXP Code, Uint nr, bool test);
void InitShuffle();
Uint CodesPerBlockShuffle();
Uint UnitsPerIndivShuffle(Uint snps);
void haplo2genoShuffle(Uint * SNPxIndiv, Uint snps, Uint individuals, Uint *A);
Ulong sumGenoShuffle(Uint *S, Uint snps, Uint individuals);

void zeroNthGenoShuffle(SEXP CM, SEXP NN);
SEXP allele_freqShuffle(SEXP GM);
SEXP get_matrixN_shuffle(SEXP CM, SEXP NN);
void ReUseAsShuffle(SEXP Code);



#endif
