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


#ifndef parallel_omp_base_H
#define parallel_omp_base_H 1

// NEVER_** is used by parallel_**, so no conflict with NO_, which is used
// in the programmes. Except that there is no differences between
// NEVER_ and NO_

// #define NEVER_OMP 1
// #define NEVER_AVX 1
// #define NEVER_SSE 1

#if defined WIN32 || defined _WIN32 || defined __WIN32__
#define MSDOS_WINDOWS 1
#elif defined MSDOS_WINDOWS
#undef MSDOS_WINDOWS
#endif


#if defined __x86_64 || defined __x86_64__ || defined __amd64__ || defined __amd64 || defined _M_X64
#define X86_64 1
#elif defined X86_64
#undef X86_64
#endif


#if defined __arm64__ || defined __arm64 || defined __aarch64__
#define ARM64 1
#elif defined ARM64
#undef ARM64
#endif

#if defined __arm32__ || defined __arm__ || defined ARM64
#define ARM32 1
#define NO_AVX 1
#endif

#if defined __ARM_NEON__ || defined __aarch64__ || defined _M_ARM || defined _M_ARM64
#define NEON 1
#elif defined NEON
#undef NEON
#endif

#if defined MSDOS_WINDOWS || defined (__APPLE__) || defined(__sun)
  #if defined TIME_AVAILABLE
    #undef TIME_AVAILABLE
  #endif
#else
#define TIME_AVAILABLE 1
#endif


#if defined _OPENMP && ! defined NO_OMP && ! defined NEVER_OMP && ! defined ARM32 && ! defined __APPLE__ // 15 Jan 2022 
  #if defined SCHLATHERS_MACHINE
    #define DO_PARALLEL 1 // may change value when debugging
  #else
    #define DO_PARALLEL 1// never changes value 
  #endif
#elif defined DO_PARALLEL
  #undef DO_PARALLEL
#endif


#if defined NEVER_SSE
  #ifndef NO_SSE2
     #define NO_SSE2 1
  #endif
#elif defined NEVER_AVX
  #ifndef NO_AVX
     #define NO_AVX 1
  #endif
#elif defined NEVER_AVX512
  #ifndef NO_AVX512
     #define NO_AVX512 1
  #endif
#endif


#if defined NO_SSE2 && ! defined NO_SSE3
#define NO_SSE3 1
#endif
#if defined NO_SSE3 && ! defined NO_SSSE3
#define NO_SSSE3 1
#endif
#if defined NO_SSSE3 && ! defined NO_SSE41
#define NO_SSE41 1
#endif
#if defined NO_SSE41 && ! defined NO_AVX
#define NO_AVX 1
#endif
#if defined NO_AVX && ! defined NO_AVX2
#define NO_AVX2 1
#endif
#if defined NO_AVX2 && ! defined NO_AVX512
#define NO_AVX512 1
#endif


#if ! defined NO_AVX512

#if ! defined DO_AVX512BITALG && ! defined DO_AVX512BW && ! defined DO_AVX512CD && ! defined DO_AVX512DQ && ! defined DO_AVX512ER && ! defined DO_AVX512F && ! defined DO_AVX512IFMA && ! defined DO_AVX512PF && ! defined DO_AVX512VBMI && ! defined DO_AVX512VL && ! defined DO_AVX512VPOPCNTDQ  && ! defined DO_AVX5124FMAPS && ! defined DO_AVX5124VNNIW
#define DO_AVX512BITALG 1
#define DO_AVX512BW 1 
#define DO_AVX512CD 1
#define DO_AVX512DQ 1
#define DO_AVX512ER 1
#define DO_AVX512F 1
#define DO_AVX512IFMA 1
#define DO_AVX512PF 1
#define DO_AVX512VBMI 1
#define DO_AVX512VL 1
#define DO_AVX512VPOPCNTDQ 1
#define DO_AVX5124FMAPS 1
#define DO_AVX5124VNNIW 1
#endif

#if defined __AVX512BITALG__ && defined DO_AVX512BITALG
#define AVX512BITALG 1
#endif
#if defined __AVX512BW__ && defined DO_AVX512BW
#define AVX512BW 1
#endif
#if defined __AVX512CD__ && defined DO_AVX512CD
#define AVX512CD 1
#endif
#if defined __AVX512DQ__ && defined DO_AVX512DQ
#define AVX512DQ 1
#endif
#if defined __AVX512ER__ && defined DO_AVX512ER
#define AVX512ER 1
#endif
#if defined __AVX512F__ && defined DO_AVX512F
#define AVX512F 1
#define AVX512 1
#endif
#if defined __AVX512IFMA__ && defined DO_AVX512IFMA
#define AVX512IFMA 1
#endif
#if defined __AVX512PF__ && defined DO_AVX512PF
#define AVX512PF 1
#endif
#if defined __AVX512VBMI__ && defined DO_AVX512VBMI
#define AVX512VBMI 1
#endif
#if defined __AVX512VL__ && defined DO_AVX512VL
#define AVX512VL 1 //
#endif
#if defined __AVX512VPOPCNTDQ__ && defined DO_AVX512VPOPCNTDQ
#define AVX512VPOPCNTDQ 1 //
#endif
#if defined __AVX5124FMAPS__ && defined DO_AVX5124FMAPS
#define AVX5124FMAPS 1
#endif
#if defined __AVX5124VNNIW__ && defined DO_AVX5124VNNIW
#define AVX5124VNNIW 1
#endif

#endif // end ! no 512


#if defined __AVX2__ && ! defined NO_AVX2
#define AVX2 1
#elif defined AVX2
#undef AVX2
#endif

#if defined __AVX__  && ! defined NO_AVX
#define AVX 1
#elif defined AVX
#undef AVX
#endif

#if (defined __SSE4_1__ || defined NEON) && ! defined NO_SSE41
#define SSE41 1
#elif defined SSE41
#undef SSE41
#endif

#if (defined __SSSE3__ || defined NEON)  && ! defined NO_SSSE3
#define SSSE3 1
#elif defined SSSE3
#undef SSSE3
#endif

#if (defined __SSE3__ || defined NEON) && ! defined NO_SSE3
#define SSE3 1
#elif defined SSE3
#undef SSE3
#endif

#if (defined __SSE2__ || defined NEON) && ! defined NO_SSE2
#define SSE2 1 
#elif defined SSE2
#undef SSE2
#endif



#endif
