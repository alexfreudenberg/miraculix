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

   
  
#include "5codesDef.h"
#define MY_VARIANT 512


#include "Basic_miraculix.h"
#include"intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "Template.h"
#include "5codesDef.h" 
#include "5codes.h"


#if defined AVX2    


ASSERT_SIMD(5codes512, avx2);

#include "MX.h"
#include "haplogeno.h" 
#include "intrinsicsIntern.h" 
#include "5codesIntern.h"
 
 


#else
#include "avx_miss.h"
  
   
gV5_header(512, floatD, LongDouble, float, double) Sv
gV5_header(512, double, LongDouble, double, double) Sv


SIMD_MISS(5codes512, avx2);

#endif  // defined AVX
// trafo2Geno1Geno128
 
