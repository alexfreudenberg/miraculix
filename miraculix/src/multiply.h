
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


#ifndef miraculix_multiply_H
#define miraculix_multiply_H 1

#include "MX.h"

Uint BytesPerBlockMultiply();
Uint CodesPerBlockMultiply();
void crossprod_multiply(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_multiply(Uint snps, Uint individuals, SEXP G);


Uint BytesPerBlockMultiply256();
Uint CodesPerBlockMultiply256();
void crossprod_multiply256(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_multiply256(Uint snps, Uint individuals, SEXP G);


#endif
