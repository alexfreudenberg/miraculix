
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

#define MY_METHOD Multiply256

#include "IntrinsicsBase.h"
#include "error.h"
#include "MX.h"
#include "options.h"

#if defined AVX2



#include "2bit.h"
#include "intrinsics.h"
#include  "options.h"
#include "xport_import.h"
#include "align.h"
#include "haplogeno.h"

#define UPI  UnitsPerIndiv256
#include "multiplyIntern.h"
MultiplyInitIntern(256)


Uint BytesPerBlockMultiply256() { return BytesPerBlock; }
Uint CodesPerBlockMultiply256() { return CodesPerBlock; }
Uint BitsPerCodeMultiply256() { return BitsPerCode; }


SEXP matrix_start_multiply256( Uint snps,Uint individuals, SEXP file) {  
  return matrix_start_Intern(snps, individuals, file);
}
 

void crossprod_multiply256(Uint *CGM, Uint snps, Uint individuals, double *ans){
  crossprodIntern(CGM, snps, individuals, ans);
}


bool useMultiply256(snpcoding method) {
  option_type *global = &(KEYT()->global);
  return method == Multiply ? true : global->genetics.efficient;  
}


#else // !defined AVX2
void static AVX2missing() { ERR("'multiply256' needs the availablity of 'AVX2'"); }
#define Sm { AVX2missing(); return R_NilValue; }
#define Su { AVX2missing(); return 0; }
#define Sv { AVX2missing(); }
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif
void crossprod_multiply256(Uint V* CGM, Uint V snps, Uint V individuals,
			double V *ans) Sv
SEXP matrix_start_multiply256(Uint V snps, Uint V individuals, SEXP V G) Sm

Uint BytesPerBlockMultiply256() Su
Uint CodesPerBlockMultiply256() Su
Uint UnitsPerIndivMultiply256(Uint V snps) Su
Uint BitsPerCodeMultiply256() Su
  void MultiplyInit256() { return; }

#endif  // defined AVX2
