
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


#ifndef miraculix_packed_H
#define miraculix_packed_H 1

#include "MX.h"


Uint BytesPerBlockPacked();
Uint CodesPerBlockPacked();
void crossprod_packed(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_packed(Uint snps,Uint individuals,  SEXP G);
bool usePacked(snpcoding method);

//Uint UnitsPerIndivPacked(Uint snps);
#define CodesPerBlock128 CodesPerBlockPacked
void haplo2geno128(Uint * SNPxIndiv, Uint snps, Uint individuals,
		   Uint upiH, Uint *A);
Ulong sumGeno128(Uint *S, Uint snps, Uint individuals);
SEXP allele_freq128(SEXP GM);
bool use128();


Uint BytesPerBlockPacked256();
Uint CodesPerBlockPacked256();
void crossprod_packed256(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_packed256(Uint snps, Uint individuals, SEXP G);

void haplo2geno256(Uint * SNPxIndiv, Uint snps, Uint individuals,
			  Uint upiH, Uint *A);
Ulong sumGeno256(Uint *S, Uint snps, Uint individuals);
SEXP allele_freq256(SEXP GM);
bool use256();

#endif
