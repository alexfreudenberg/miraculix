/*
 Authors 
 Martin Schlather, Alexander Freudenberg, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Martin Schlather, Alexander Freudenberg

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/


#include "Basic_miraculix.h"
#include "xport_import.h"
#include "mmagpu.h"
#include "options.h"
#include "extern.h"

#define GPUmissing ERR0("This installation of miraculix hasn't been configured to use GPUs.")
 
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif


#if ! defined USEGPU 
void crossprod_PlainGPU(Uint V *M, Long V snps, Long V individuals,
			Uint V Vlda, int V cores, double V *A ){ GPUmissing }
SIMD_MISS(mma_61, gpu);
#endif


#if ! defined USEGPU || ! defined DO_PARALLEL
void check_7_5() {
  basic_options opt;
  Ext_get_utils_basic(&opt, false);
  if (opt.Rprintlevel > 1) PRINTF("Note that the package has not been compiled adequately to use the GPU. See the start-up messages for details.");
}
void crossprod_mmagpu(Uint V *code,
		      Long V snps,
		      Long V individuals,
		      Uint V lda,
		      int V cores,
		      double V *ans) { GPUmissing }

SIMD_MISS(mma_75, gpu);

#else
void check_7_5() { };

#endif
