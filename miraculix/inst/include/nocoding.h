
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


#ifndef miraculix_auto_H
#define miraculix_auto_H 1

SEXP matrix_start_plain(Uint snps, Uint individuals,  SEXP G);
void coding_plain(Uint *M, Uint start_individual, Uint end_individual, 
		  Uint start_snp, Uint end_snp,  Uint Mnrow, SEXP Ans,
		  double *G);
SEXP get_matrixPlain(SEXP SNPxIndiv);
Uint *AlignPlain(SEXP Code, Uint nr, bool test);
Uint BytesPerBlockPlain();
Uint CodesPerBlockPlain();
Uint UnitsPerIndivPlain(Uint snps);
Uint BitsPerCodePlain();
Ulong sumGenoPlain(Uint *S, Uint snps, Uint individuals);
void haplo2genoPlain(Uint * SNPxIndiv, Uint snps, Uint individuals,
		     Uint unitsPerIndiv, Uint *A);
Ulong sumGenoPlain(Uint *S, Uint snps, Uint individuals);
void crossprod_Plain(Uint * SNPxIndiv, Uint snps, Uint individuals, double *A);

#endif

