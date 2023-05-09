/*
 Authors 
 Guido Moerkotto
 maintained and modified by Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Guido Moerkotte, Martin Schlather

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

#define BitsPerCode 1
#define MY_VARIANT 256

#include "Basic_miraculix.h"
#include"intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "bitBase.h"


#if defined AVX2

ASSERT_SIMD(1bit256, avx2);

#include "MX.h"
#include "haplogeno.h"
#include "intrinsicsIntern.h"
#include "1bitIntern.h"


SCALAR012(256)



#else // !defined AVX
#include "Template.h"
#include "1bit.h"
#include "avx_miss.h"


  DeclareVersion1(256,Su,Sv)
  SIMD_MISS(1bit256, avx2);

#endif  // defined AVX2
