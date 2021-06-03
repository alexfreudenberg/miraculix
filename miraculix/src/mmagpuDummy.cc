
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


#include "error.h"
#include "mmagpu.h"
#include "options.h"
#include "xport_import.h"

#ifndef USEGPU
#define GPUmissing ERR("This installation of miraculix hasn't been configured to use GPUs.")
 
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif


bool check_7_5() {
  if (PL > 1) PRINTF("Note that the package has not been compiled adequately to use the GPU. See the start-up messages for details.");
  return false;
}

SEXP matrix_start_mmagpu(Uint V snps, Uint V individuals,  SEXP V G) {
  GPUmissing; return R_NilValue;}

void crossprod_mmagpu(Uint V * CGM,
		      Uint V snps,
		      Uint V individuals,
		      double V *ans) { GPUmissing }

void coding2gpu(Uint V *M,
		Uint V start_individual,
		Uint V end_individual, 
		Uint V start_snp,
		Uint V end_snp,
		Uint V Mnrow,
		SEXP V Ans,
		double V *G){ GPUmissing }

bool useMMAGPU(snpcoding method) {
  option_type *global = &(KEYT()->global);
  if (global->genetics.efficient || method != MMAGPU) return false;
  ERR("'CUDA' is needed for vector multiplication in 'MMAGPU'. Set 'RFoptions(efficient=TRUE)'.")
}

#endif
