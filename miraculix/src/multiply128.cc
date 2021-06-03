
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

#define MY_METHOD Multiply

#include "IntrinsicsBase.h"
#include "error.h"
#include "MX.h"
#include "options.h"

#if defined SSE2


#if defined AVX2
#undef AVX2
#endif

#if defined AVX512
#undef AVX512
#endif



#include "2bit.h"
#include "intrinsics.h"
#include  "options.h"
#include "xport_import.h"
#include "align.h"
#include "haplogeno.h"


#define UPI  UnitsPerIndiv256
#include "multiplyIntern.h"
MultiplyInitIntern(128)

Uint BytesPerBlockMultiply() { return BytesPerBlock; }
Uint CodesPerBlockMultiply() { return CodesPerBlock; }
Uint BitsPerCodeMultiply() { return BitsPerCode; }



SEXP matrix_start_multiply( Uint snps,Uint individuals, SEXP file) {  
  return matrix_start_Intern(snps, individuals, file);
}
 

void crossprod_multiply(Uint *CGM, Uint snps, Uint individuals, double *ans){
  crossprodIntern(CGM, snps, individuals, ans);
}



#else // !defined AVX2
void static AVX2missing() { ERR("'multiply' needs the availablity of 'AVX2'"); }
#define Sm { AVX2missing(); return R_NilValue; }
#define Su { AVX2missing(); return 0; }
#define Sv { AVX2missing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
void crossprod_multiply(Uint V* CGM, Uint V snps, Uint V individuals,
			double V *ans) Sv
SEXP matrix_start_multiply(Uint V snps, Uint V individuals, SEXP V G) Sm


Uint BytesPerBlockMultiply() Su
Uint CodesPerBlockMultiply() Su
Uint UnitsPerIndivMultiply(Uint V snps) Su


#endif  // defined AVX2
