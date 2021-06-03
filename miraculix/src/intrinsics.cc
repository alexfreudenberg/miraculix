
 
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

#include <cpuid.h>
#include "IntrinsicsBase.h"
#include "intrinsics.h"
#include "AutoMiraculix.h"
//#include <Basic_utils.h>7
#include "options.h"
#include "MX.h"
#include "xport_import.h"
#include <zzz_RandomFieldsUtils.h>


AVAILABLE

#ifdef _WIN32
#include <limits.h>
#include <intrin.h>
typedef unsigned __int32  uint32_t;
#else
#include <stdint.h>
#endif


#if __GNUC__ >= 7
#define FALLTHROUGH_OK __attribute__ ((fallthrough));
#else
#define FALLTHROUGH_OK   
#endif


// https://github.com/m-j-w/CpuId.jl/blob/master/src/cpufeature.jl
// https://stackoverflow.com/questions/1666093/cpuid-implementations-in-c
// https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?view=msvc-160
//https://cpp.hotexamples.com/examples/-/-/__get_cpuid/cpp-__get_cpuid-function-examples.html


snpcoding check_method(snpcoding method) {
  // check whether compiled code is run on a compatible kernel
  snpcoding unmet_method = UnknownSNPcoding;
  // if (method != Haplo) printf("method %d", method);
#if defined USEGP
  if (!avx2)
    ERR1("modern graphic card, but too old kernel. If this is an issue: %200s",
       CONTRACT);
  return method;
#elif defined AVX2
  if (avx2 || method == Haplo) return method;
  unmet_method = Shuffle256;
#elif defined SSSE3
  if (ssse 3|| method == Haplo) return method; 
  unmet_method = Shuffle;
#elif defined SSE2
  if (sse2|| method == Haplo) return method; 
  unmet_method = Packed;
#else
 return method; 
#endif

  //  code runs on a incompatible kernel !!
  option_type *global = &(KEYT()->global);
  bool switch_ok = !global->genetics.efficient && is256(method);

  //  printf("switchok = %d %d %d unmet=%d\n", switch_ok, ssse3, sse2, unmet_method);
  
  switch(unmet_method) {
  case Shuffle256 :
    if (ssse3) {
      if ((method >= FirstMoBPSmethod && method <= Shuffle)
	  || method <= Hamming3 ) return method;
      if (switch_ok && method >= FirstMoBPSmethod && method <= Shuffle256) {
	if (PL > 2)
	  PRINTF("kernel change: '%20s' not recognized. Switch to %20s\n",
		 SNPCODING_NAMES[method], SNPCODING_NAMES[Shuffle]);  
	return Shuffle;
      }
    }
    FALLTHROUGH_OK;
  case Shuffle :
    if (sse2) {
      if ((method >= FirstMoBPSmethod && method <= Packed)
	  || method <= Hamming2) return method;
      if (switch_ok && method >= FirstMoBPSmethod && method <= Shuffle256) {
	if (PL > 2)
	  PRINTF("kernel change: '%20s' not recognized. Switch to %20s\n",
		 SNPCODING_NAMES[method], SNPCODING_NAMES[Packed]);  
	return Packed;
      }
    }
    FALLTHROUGH_OK;
 case Packed :
   if ((method >= FirstMoBPSmethod && method <= TwoBit)
       || method <= ThreeBit) return method;
   if (switch_ok && method >= FirstMoBPSmethod && method <= Shuffle256) {
     if (PL > 2)
       PRINTF("kernel change: '%20s' not recognized. Switch to %20s\n",
	      SNPCODING_NAMES[method], SNPCODING_NAMES[Packed]);  
     return TwoBit;
   }
   break;
  default : BUG;
  }

  ERR("The machine was compiled with different flags than the currant kernel understands. Set RFoptions(genetics.efficient=FALSE, genetics.snpcoding=Shuffle256).");
  
  return TwoBit;
}
