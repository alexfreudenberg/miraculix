/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Martin Schlather

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



#define BitsPerCode 2
#define MY_VARIANT 512


#include "Basic_miraculix.h"
#include"intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "Template.h"
#include "Files.h"
#include "MX.h"
#include "2bit.h"

#if defined AVX512

ASSERT_SIMD(2bit512, avx512f);

#include "haplogeno.h"
#include "intrinsicsIntern.h"
#include "2bitIntern.h"


void coding_2v512(unit_t *M, Long ldM,
		  Long start_individual, Long end_individual, 
		  Long start_snp, Long end_snp, 
		  basic_options VARIABLE_IS_NOT_USED  *opt,
		  double VARIABLE_IS_NOT_USED *G, SEXP Ans) {
  /*
    use 
__m256i _mm256_maskz_compress_epi8 (__mmask32 k, __m256i a)
void _mm512_mask_compressstoreu_epi8 (void* base_addr, __mmask64 k, __m512i a)
  
 */
  
  BUG;
}

FileBinary(512);
SCALAR012(512)

#else // !defined AVX
#include "avx_miss.h"


DeclareVersion2(512,Su,Sv,Sm)
SIMD_MISS(2bit512, avx512f);


#endif  // defined AVX2
// trafo2Geno1Geno128
