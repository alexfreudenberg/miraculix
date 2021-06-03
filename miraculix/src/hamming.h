
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


#ifndef miraculix_sse_H
#define miraculix_sse_H 1

Uint BytesPerBlockH();
Uint CodesPerBlockH();
Uint BitsPerCodeH();

SEXP zeroNthGenoH(SEXP SNPxIndiv, SEXP Snps, snpcoding method);
Uint *AlignH(SEXP Code, Uint nr, bool test);


Uint UnitsPerIndivH2(Uint snps);
SEXP matrix_startH2( Uint snps,Uint individuals, SEXP file);
void codingH2(Uint *M, Uint start_individual, Uint end_individual, 
	     Uint start_snp, Uint end_snp, Uint Mnrow, SEXP Ans,
	     double *G);
void crossprod_H2(Uint *SNPxIndiv, Uint snps, Uint individuals, double *A);
void haplo2genoH2(Uint *SNPxIndiv, Uint snps, Uint individuals, Uint U,
		  Uint *A);
SEXP get_matrixH2(SEXP SNPxIndiv);
Ulong sumGenoH2(Uint *M, Uint snps, Uint individuals);
SEXP allele_freqH2(SEXP SNPxIndiv);
SEXP get_matrixN_H2(SEXP SNPxIndiv, SEXP Snps);


Uint UnitsPerIndivH3(Uint snps);
SEXP matrix_startH3(Uint snps, Uint individuals, SEXP file);
void codingH3(Uint *M, Uint start_individual, Uint end_individual, 
	     Uint start_snp, Uint end_snp, Uint Mnrow, SEXP Ans,
	     double *G);
void crossprod_H3(Uint *SNPxIndiv, Uint snps, Uint individuals, double *A);
void haplo2genoH3(Uint *SNPxIndiv, Uint snps, Uint individuals, Uint U,
		  Uint *A);
SEXP get_matrixH3(SEXP SNPxIndiv);
Ulong sumGenoH3(Uint *M, Uint snps, Uint individuals);
SEXP allele_freqH3(SEXP SNPxIndiv);
SEXP get_matrixN_H3(SEXP SNPxIndiv, SEXP Snps);





#endif
