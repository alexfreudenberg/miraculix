
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

Uint BytesPerBlockShuffle();
Uint CodesPerBlockShuffle();
void crossprod_shuffle(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_shuffle(Uint snps, Uint individuals, SEXP G);
bool useShuffle(snpcoding method);

Uint BytesPerBlockShuffle256();
Uint CodesPerBlockShuffle256();
void crossprod_shuffle256(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_shuffle256(Uint snps, Uint individuals, SEXP G);
bool useShuffle256(snpcoding method);


#endif
