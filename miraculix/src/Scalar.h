


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2019 -- 2019 Martin Schlather

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


#ifndef miraculix_ScalarInt_H
#define miraculix_ScalarInt_H 1

#include "MX.h"
#include "IntrinsicsBase.h"
#include "xport_import.h"

#define SCALAR_INT_8 0
#define SCALAR_INT_16 1
#define SCALAR_INT_AVX2 2

#if defined AVX2
//
#define SCALAR_INT_DEFAULT SCALAR_INT_AVX2
//#define SCALAR_INT_DEFAULT SCALAR_INT_16
#else
#define SCALAR_INT_DEFAULT SCALAR_INT_8
#endif

#define SCALARUINT(x, y, len) scalarUint(x, y, len, SCALAR_INT_DEFAULT);
Uint scalarUint(Uint *x, Uint *y, Uint len, Uint n);
void matmulttransposedUint(Uint *A, Uint *B, double *c, Uint m, Uint l, Uint n);

#endif
