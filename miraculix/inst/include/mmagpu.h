
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2020 -- 2021  Martin Schlather, Alexander Freudenberg

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


#ifndef miraculix_mmagpu_H
#define miraculix_mmagpu_H 1

#include "MX.h"


bool useMMAGPU(snpcoding method);
void crossprod_mmagpu(Uint* CGM, Uint snps, Uint individuals, double *ans);
SEXP matrix_start_mmagpu(Uint snps,Uint individuals,  SEXP G);

void coding2gpu(Uint *M, Uint start_individual, Uint end_individual, 
	       Uint start_snp, Uint end_snp, Uint Mnrow, SEXP Ans,
	       double VARIABLE_IS_NOT_USED *G);
bool check_7_5();



#endif
