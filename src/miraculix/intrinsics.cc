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
 

/*
#if defined MINGWCPUID
#include <cpuid.h>
#if defined( __MINGW64__ )
    #include <xmmintrin.h>
#endif
#if defined WIN32 || defined _WIN32 || defined __WIN32__
#include <limits.h>
#include <intrin.h>
typedef unsigned __int32  uint32_t;
#else
#include <stdint.h>
#endif

*/

#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "MXinfo.h"
#include "extern.h"


AVAILABLE_SIMD

#define check_variant_error -1
int // NICHT Uint!
check_variant(coding_type coding, int variant,
	      bool efficient, bool stopOnError) {
  // check whether compiled code is run on a compatible kernel
  // if (coding != Haplo) printf("coding %d", coding);
  //  code runs on a incompatible kernel !!

 	 
   bool
    ok = true,
    switch_ok = (efficient && stopOnError) || variant == VARIANT_DEFAULT;


   //   printf("check variant %sv%d effi=%d switch=%d\n", CODING_NAMES[coding], variant, efficient, switch_ok);


#if defined USEGPU
  if (variant == VARIANT_GPU) {
    if (!avx2) {
      if (!stopOnError) return check_variant_error;
      ERR1("modern graphic card, but too old kernel. If this is an issue:%200s",
	   CONTACT);
   }
    return variant;
  }
#else
  if (variant == VARIANT_GPU) {
    if (!switch_ok) {
      if (!stopOnError) return check_variant_error;
      ERR0("graphic card not available");
    }
    variant = 0;
  }  
#endif

  //  if (coding > 2) printf("%s %d\n", CODING_NAMES[coding], variant);
  assert(FirstUserHaploCoding == 20 && LastUserHaploCoding == 24);

  if (variant == VARIANT_DEFAULT || //variant == VARIANT_R ||
      efficient) {
    switch(coding) {
    case OneBitGeno : variant = 512 + VARIANT_A; break;
    case TwoBitGeno :
    case TwoBitGenoTransposed :
      variant =
#if defined USEGPU
	VARIANT_GPU;
#else
      256;
#endif
      break;
    case ThreeBit : variant = 64;break;
    case FourBit : BUG;
    case OneByteGeno : variant = 32; break;
    case TwoByte : BUG;
    case FourByteGeno :
      variant =
#if defined USEGPU
	VARIANT_GPU;
#else
	32;
#endif
      break;
    case Plink : variant = 32; break;
    case PlinkTransposed : variant = 256; break;
    case OrigPlink : variant = 257; break;
    case OrigPlinkTransposed : variant = 257; break;
     
    case FiveCodes: 
    case FiveCodesTransposed: variant = 256; break;
    case FourByte_Tr: variant = 32; break;
    case OneBitHaplo: variant = 32; break;
    case TwoBitHaplo: variant = 32; break;
    case OneByteHaplo: variant = 32; break;
    case FourByteHaplo: variant = 32; break;
    case FourByteSingleBit: variant = 32; break;
    case EightByteHaplo: variant = 32; break;
    case FourByteSingleBit_Tr: variant = 32; break;
    case EightByteHaplo_Tr: variant = 32; break;
    default :
      //printf("intrinics %d\n", coding);    printf("intrinics %s\n", CODING_NAMES[coding]);
      BUG;
    }
  }
  
#if defined USEGPU
  if (variant == VARIANT_GPU) return variant;
#endif
  if (coding == FourByteGeno) return variant;
  
#if defined AVX512
  if (((ok = sse41Avail) && variant < 256)) 
 if ((ok = avx512fAvail)) {
    if (coding != OneBitGeno) return variant;
  #if defined AVX512POPCNT
    if ((ok = avx512popcnt)) return variant;
  #else   
    if (variant <= 512) return variant;
  #endif
  }
#elif defined AVX2
  if (((ok = avx2Avail) && variant < 512) )
    return variant;
#elif defined SSE41
  if (((ok = sse41Avail) && variant < 256))
    return variant; 
 #elif defined SSSE3
  if (((ok = ssse3Avail) && variant < 256))
    return variant;  
#elif defined SSE2
  if (((ok = sse2Avail) && variant < 128))
    return variant; 
#else
   if (variant < 128) return variant; 
#endif

  
  if (switch_ok) {
    switch(coding) {
    case OneBitGeno :
#if defined AVX512F  
      if (avx512fAvail) return 512;
#elif defined AVX2
      if (avx2Avail) return 256;
#elif defined SSSE3
      if (ssse3Avail) return 128;
#endif
      return 64;
    case TwoBitGeno :
    case TwoBitGenoTransposed :
#if defined AVX2
      if (avx2Avail) return 256; 
#elif defined SSSE3
      if (ssse3Avail) return 128; 
#endif
      return 64;
    case FiveCodes :
    case FiveCodesTransposed: 
#if defined AVX2
      if (avx2Avail) return 256; 
#elif defined SSE41
      if (sse41Avail) return 128;
#endif
      return 32;
    default:
      BUG; // ThreeBit or 
    }
  }
   
  if (!stopOnError) return check_variant_error;
  if (switch_ok) BUG;
  if (ok) {
    ERR3("variant '%d' for '%.50s' not compiled. %.50s",
	 variant,
	  CODING_NAMES[coding],
	 efficient ? "Bug?" : "Set basic.efficient=TRUE?"); // OK
  }
  if (false) {
    ERR6("The machine was compiled with different flags than the currant kernel understands. Current RFoptions are basic.efficient=%d, snpcoding=%.50s, variant=%d (%d; %d %d)",
	 efficient, CODING_NAMES[coding],
	 variant, sse2Avail, switch_ok, stopOnError);
  } else {
    ERR3("The machine was compiled with different flags than the currant kernel understands. Current RFoptions are basic.efficient=%d, snpcoding=%.50s, variant=%d",
	 efficient, CODING_NAMES[coding], variant);
  }
  return TwoBitGeno;
}

int check_variant(coding_type coding, int variant, bool efficient) {
  return check_variant(coding, variant, efficient, true);
}


#if defined SSE2 
bool any128(__m128i A) {
  uint64_t *x = (uint64_t*) (&A);
  return x[0] || x[1];
}

uint64_t ADD8plain64(uint64_t B, uint64_t C) { // safe ADD8, but slow
  uint64_t ans;
  unsigned char
    *x=(unsigned char*) &B,
    *y=(unsigned char*) &C,
    *z=(unsigned char*) &ans;
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
  z[3] = x[3] + y[3];
  z[4] = x[4] + y[4];
  z[5] = x[5] + y[5];
  z[6] = x[6] + y[6];
  z[7] = x[7] + y[7];
  return ans;
}

uint32_t ADD8plain32(uint32_t B, uint32_t C) { // safe ADD8, but slow
  uint32_t ans;
  unsigned char
    *x=(unsigned char*) &B,
    *y=(unsigned char*) &C,
    *z=(unsigned char*) &ans;
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
  z[3] = x[3] + y[3];
  return ans;
}


#endif
