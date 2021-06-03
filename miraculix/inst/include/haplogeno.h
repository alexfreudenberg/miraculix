


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 -- 2019 Martin Schlather

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


#ifndef miraculix_haplogeno_H
#define miraculix_haplogeno_H 1

#include <inttypes.h>
#include "miraculix.h"
#include "AutoMiraculix.h"
// #include "IntrinsicsBase.h" --- darf nicht rein !!!
#include "intrinsics.h"

Uint *AlignBase(SEXP CM, Uint nr, Uint bytesperblock, bool test);

//Uint* DoAlign(SEXP SNPxIndiv, Uint nr, snpcoding method);
//Uint* DoAlignWithoutTest(SEXP CM, Uint nr, snpcoding method);
SEXP createSNPmatrix(Uint snps, Uint individuals, snpcoding method);
Ulong sumGeno(Uint * SNPxIndiv, Uint snps, Uint individuals, snpcoding method);

Uint haplo2geno(Uint * SNPxIndiv, Uint snps, Uint individuals,
		snpcoding method, Uint unitsPerIndivH, Uint *A);
void haplo2geno(SEXP H, Uint *code, snpcoding method);
double *crossprod(Uint * SNPxIndiv, Uint snps, Uint individuals,
		    snpcoding method,  bool centred, bool normalized,
		    Ulong SumGeno);

Uint GetBytesPerBlock(snpcoding method);
Uint GetCodesPerBlock(snpcoding method);
Uint GetBitsPerCode(snpcoding method);
Uint GetUPI(Uint snps, snpcoding method);


void ReUseAs(SEXP Code, snpcoding method);
SEXP CreateEmptyCodeVector(Uint snps, Uint individuals, snpcoding method);
void start_info(SEXP Code, SEXP file);
void allInfo(Uint * info);
void allInfo(SEXP M);

Ulong calculateAlignedMem(Ulong memInUnits, snpcoding method,
			  Uint bytesPerBlock);


#endif

