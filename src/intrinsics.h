/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2018  Martin Schlather

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


#ifndef miraculix_initrinsics_H
#define miraculix_initrinsics_H 1


// PKG_CXXFLAGS =  $(SHLIB_OPENMP_CXXFLAGS) -mavx -msse -msse2 -msse3 -msse4 -msse4a -march=core-avx2

#if defined __AVX2__
//#define AVX2 __AVX2__
#endif
#if defined __AVX__
#define AVX 1
#endif
#if defined __SSE41__ 
#define SSE412 1
#endif
#if defined  __SSE4A__
#define SSE4A 1
#endif
#if defined  __SSSE3__ 
#define SSSE3 1
#endif
#if defined  __SSE3__ 
#define SSE3 1
#endif
#if defined  __SSE2__ 
#define SSE2 1
#endif
#if defined __SSE__ 
#define SSE 1
#endif
#if defined __MMX__ 
#define MMX 1
#endif
//

#if defined MMX
//#include <mmintrin.h>
#endif

#if defined SSE
//#include <xmmintrin.h>
#endif

#if defined SSE2
#include <emmintrin.h>
#endif

#if defined SSE3
//#include <pmmintrin.h>
#endif

#if defined SSSE3
#include <tmmintrin.h>
#endif

#if defined SSE4A
//#include <ammintrin.h>
#endif

#if defined SSE412
//#include <smmintrin.h>
#endif

#if defined AVX or defined AVX2
#include <x86intrin.h>
#endif

#if defined AVX512
//#include <immintrin.h>
#endif

#if defined (AVX512)
#define SSEBITS 512
#define SSEMODE 30
#elif defined (AVX)
#define SSEBITS 256
#define SSEMODE 20
#elif defined (SSE)
#define SSEBITS 128
#define SSEMODE 10
#else
#define SSEBITS 64
#define SSEMODE 0
#endif

#define ALLIGNED __declspec(align(SSEBITS/8))


#define bitsperbyte 8
#define Uint unsigned int
#define uint64 uint64_t

#ifndef INT64_C
#define INT64_C(X) X
#endif


#endif




