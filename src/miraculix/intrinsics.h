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



#ifndef rfutils_intrinsics_H
#define rfutils_intrinsics_H 1

#ifdef PRINTF
#error "intrinsics.h not very first"
#endif

#include<inttypes.h> // uintptr_t
#include "parallel_simd.h"


#if defined MINGWCPUID
#include <cpuid.h>
#elif defined WINCPUID
//#warning loading intrin.h the first time
#include <intrin.h>
#endif

//#if defined _ _ARM_NEON
//#include <arm_neon.h>
//#if defined(_ _LP64_ _) && _ _LP64_ _
//#endif
//#endif

#if ! defined MEMisALIGNED
#define MEMisALIGNED Nan
#endif

#if defined ARM32 && defined SSE2
#include "sse2neon.h"
#elif defined AVX || defined SSE2 //|| defined AVX2 || defined 
#include <immintrin.h>
#endif

#if __GNUC__ > 4 ||							\
  (__GNUC__ == 4 && (__GNUC_MINOR__ > 9 ||				\
		     (__GNUC_MINOR__ == 9 &&  __GNUC_PATCHLEVEL__ >= 1)))
//#define OpenMP4 1
#endif


union uni32{
  uint32_t vi;
  float f4[1];
  uint32_t u32[1];
  uint16_t u16[2];
  uint8_t u8[4];
};

union uni64{
  uint64_t vi;
  uint64_t u64[1];
  double d8[1];
  float f4[2];
  uint32_t u32[2];
  uint16_t u16[4];
  uint8_t u8[8];
};

union uni128{
#if defined SSE2
  __m128i vi;
  __m128d d;
  __m128 f;
  __m128d d128[1];
#endif
  uint64_t u64[2];
  uint32_t u32[4];
  uint8_t u8[16];
  //  __m64 m64[2];
  double halfd[2], d8[2];
  float halff[4], f4[4];
};


union uni256 {
#if defined AVX2
  __m256i vi; 
  __m256d d;
  __m256 f;
#endif
#if defined SSE2 || defined AVX2
  __m128i i128[2];
  __m128d d128[2];
  __m128d halfd[2];
  __m128 halff[2];
#endif
  uint64_t u64[4];
  uint32_t u32[8];
  uint8_t u8[32]; 
  //  __m64 m64[4];
  double d8[4];
  float f4[8];
};

union uni512 {
#if defined AVX512
  __m512i vi; 
  __m512d d;
  __m512 f;
#endif
#if defined AVX2 || defined AVX512
  __m256i i256[2];
  __m256d d256[2];
  __m256d halfd[2];
  __m256 halff[2];
#endif
#if defined SSE2 || defined AVX2 || defined AVX512
  __m128i i128[4];
  __m128d d128[4];
  __m128 f128[4];
#endif
  uint64_t u64[8];
  uint32_t u32[16];
  uint8_t u8[64]; 
  //  __m64 m64[4];
  double d8[8];
  float f4[16];
};


#define BitsPerByte 8U


#if defined AVX512
#define SIMD_AVAILABILITY avx512f
#define SSEBITS 512U
#define SSEMODE 30U

#define BlockType0 __m512i 
#define BlockType __m512i ALIGNED
#define UnionType0 uni512
#define Doubles __m512d
#define Floats __m512

#define LOADuDOUBLE _mm512_loadu_pd
#define LOADuFLOAT _mm512_loadu_ps
#define LOADU  _mm512_loadu_si512 // _mm512_lddqu_si512
#if defined MEM_IS_ALIGNED
  #define LOADDOUBLE _mm512_load_pd
  #define LOADFLOAT _mm512_load_ps
  #define LOAD  _mm512_load_si512
#else
  #define LOAD LOADU
  #define LOADFLOAT LOADuFLOAT
  #define LOADDOUBLE LOADuDOUBLE
#endif
#define MAXDOUBLE _mm512_max_pd
#define ADDDOUBLE  _mm512_add_pd
#define SUBDOUBLE  _mm512_sub_pd
#define MULTDOUBLE _mm512_mul_pd 
#define STOREuDOUBLE _mm512_storeu_pd
#define ZERODOUBLE _mm512_setzero_pd
#define MULTFLOAT  _mm512_mul_ps 
#define ADDFLOAT  _mm512_add_ps 
#define SUBFLOAT   _mm512_sub_ps 
#define ZEROFLOAT _mm512_setzero_ps
#define STOREuFLOAT _mm512_storeu_ps
//#define BLENDFLOAT  _mm256_blend_ps
//#define DUPLICATEFLOAT  _mm512_moveldup_ps
#define MASK0ADDDOUBLE(A,M,B)  _mm512_maskz_add_pd(A, M, A, B)
//
#define BLENDvDOUBLE  _mm512_blendv_pd
#define DUPLICATEDOUBLE  _mm512_movedup_pd

#define MAXINTEGER _mm512_max_epi32
#define AND  _mm512_and_si512
#define OR  _mm512_or_si512
#define XOR  _mm512_xor_si512
#define ANY(A) (! _mm512_kortestz(_mm512_test_epi32_mask(A, A), _mm512_test_epi32_mask(A, A)))
#define SHR32  _mm512_srli_epi32 // see also _mm512512_rol_epi64,
#define SHL32  _mm512_slli_epi32
#define SHR16  _mm512_srli_epi16
#define SHR64  _mm512_srli_epi64
#define SHL64  _mm512_slli_epi64


#define SET16  _mm512_set1_epi16
#define SET32  _mm512_set1_epi32
#define SET64  _mm512_set1_epi64 // oder _m512d _mm512_set1_pd (double a)
#define ZERO   _mm512_setzero_si512
#define STORE_DOUBLE _mm512_store_pd
//#define EXTRACT16  _mm512_extract_epi16

#define ADD32  _mm512_add_epi32
#define MADD16  _mm512_madd_epi16 
#define ADD64  _mm512_add_epi64
#define MULT32  _mm512_mullo_epi32

#define SET8  _mm512_set1_epi8 // nicht! BW
#define SETREV8(  B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0) \
  _mm512_set_epi8(B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0, \
		  B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0, \
		  B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0, \
		  B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0)

#if defined AVX512BW
  #define ADD8 _mm512_add_epi8
  #define SAD8  _mm512_sad_epu8
  #define SHUFFLE8 _mm512_shuffle_epi8
#else
  #define LOWER256(A) (__m256i) _mm512_extractf64x4_pd((__m512d) (A), 0)
  #define UPPER256(A) (__m256i) _mm512_extractf64x4_pd((__m512d) (A), 1)
  #define DO_256(X, A, B) \
    _mm512_inserti64x4(_mm512_zextsi256_si512(X(LOWER256(A), LOWER256(B))), \
		       X(UPPER256(A), UPPER256(B)), 1)

  #define ADD8(A, B) DO_256(_mm256_add_epi8, A, B)		
  #define SAD8(A, B) DO_256(_mm256_sad_epu8, A, B)
  #define SHUFFLE8(A, B) DO_256(_mm256_shuffle_epi8, A, B) 
#endif


#elif defined AVX
#define SSEBITS 256U
#define SSEMODE 20U

#define BlockType0 __m256i 
#define BlockType __m256i ALIGNED
#define UnionType0 uni256
#define Doubles __m256d
#define Floats __m256


#define LOADuDOUBLE _mm256_loadu_pd
#define LOADuFLOAT _mm256_loadu_ps
#if defined MEM_IS_ALIGNED
  #define LOADDOUBLE _mm256_load_pd
  #define LOADFLOAT _mm256_load_ps
#else
  #define LOADDOUBLE LOADuDOUBLE
  #define LOADFLOAT LOADuFLOAT
#endif
#define MAXDOUBLE _mm256_max_pd
#define ADDDOUBLE  _mm256_add_pd
#define SUBDOUBLE  _mm256_sub_pd
#define MULTDOUBLE _mm256_mul_pd 
#define STOREuDOUBLE _mm256_storeu_pd
#define ZERODOUBLE _mm256_setzero_pd

#define MULTFLOAT  _mm256_mul_ps 
#define ADDFLOAT  _mm256_add_ps 
#define SUBFLOAT   _mm256_sub_ps 
#define ZEROFLOAT _mm256_setzero_ps
#define BLENDFLOAT  _mm256_blend_ps
#define STOREuFLOAT _mm256_storeu_ps
#define STORE_FLOAT _mm256_store_ps
#define DUPLICATEFLOAT  _mm256_moveldup_ps
#define MASK0ADDDOUBLE(A,M,B)  _mm256_maskz_add_pd(A, M, A, B)
#define BLENDDOUBLE  _mm256_blend_pd
#define BLENDvDOUBLE  _mm256_blendv_pd
#define DUPLICATEDOUBLE  _mm256_movedup_pd
#define BROADCAST _mm256_broadcast_sd 

#if defined AVX2
#define LOADU _mm256_loadu_si256 // _mm256_lddqu_si256
#if defined MEM_IS_ALIGNED
  #define LOAD _mm256_load_si256 // _mm256_lddqu_si256
#else
  #define LOAD LOADU
#endif

#define MAXINTEGER _mm256_max_epi32

#define AND  _mm256_and_si256
#define OR  _mm256_or_si256
#define XOR  _mm256_xor_si256
#define ANY(A) (!_mm256_testz_si256(A, A))
#define SHR32  _mm256_srli_epi32 // see also _mm256512_rol_epi64,
#define SHL32  _mm256_slli_epi32
#define SHR16  _mm256_srli_epi16
#define SHR64  _mm256_srli_epi64
#define SHL64  _mm256_slli_epi64
#define SHUFFLE8 _mm256_shuffle_epi8
#define GATHER_FLOAT _mm256_i32gather_ps
#define GATHER_DOUBLE _mm256_i64gather_pd

#define SET8  _mm256_set1_epi8
#define SETREV8(   B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0) \
  _mm256_setr_epi8(B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0, \
		   B15,B14,B13,B12,B11,B10,B9,B8,B7,B6,B5,B4,B3,B2,B1,B0)

#define SET16  _mm256_set1_epi16
#define SET32  _mm256_set1_epi32
#define SET64  _mm256_set1_epi64x // oder _m256d _mm256_set1_pd (double a)
#define ZERO   _mm256_setzero_si256
#define STORE_DOUBLE _mm256_store_pd
#define EXTRACT16  _mm256_extract_epi16

#define ADD8 _mm256_add_epi8
#define ADD32  _mm256_add_epi32
#define MADD16  _mm256_madd_epi16 
#define ADD64  _mm256_add_epi64
#define SAD8   _mm256_sad_epu8
#define MULT32  _mm256_mullo_epi32
#define SIMD_AVAILABILITY avx2

#else
  #define SIMD_AVAILABILITY avx
  #define MAXINTEGER _mm_max_epi32
#endif



#elif defined SSE2
#define SSEBITS 128U
#define SSEMODE 10U
#define BlockType0 __m128i
#define BlockType __m128i ALIGNED
#define UnionType0 uni128
#define Doubles __m128d
#define Floats __m128

#define LOADU  _mm_loadu_si128
#define LOADuDOUBLE _mm_loadu_pd
#define LOADuFLOAT _mm_loadu_ps
#if defined MEM_IS_ALIGNED
  #define LOADDOUBLE _mm_load_pd
  #define LOADFLOAT _mm_load_ps
  #define LOAD  _mm_load_si128
#else
  #define LOAD LOADU
  #define LOADDOUBLE LOADuDOUBLE
  #define LOADFLOAT LOADuFLOAT
#endif

#define MAXDOUBLE _mm_max_pd
#define MAXINTEGER _mm_max_epi32
#define ADDDOUBLE  _mm_add_pd
#define SUBDOUBLE  _mm_sub_pd
#define MULTDOUBLE _mm_mul_pd 
#define STOREuDOUBLE _mm_storeu_pd
#define ZERODOUBLE _mm_setzero_pd


#define MULTFLOAT  _mm_mul_ps 
#define ADDFLOAT  _mm_add_ps 
#define SUBFLOAT   _mm_sub_ps 
#define ZEROFLOAT _mm_setzero_ps
#define BLENDFLOAT  _mm_blend_ps
#define DUPLICATEFLOAT  _mm_moveldup_ps
#define STOREuFLOAT _mm_storeu_ps


#define AND  _mm_and_si128
#define OR  _mm_or_si128
#define XOR  _mm_xor_si128
bool any128(__m128i A);
#define ANY(A) any128(A)
#define SHR32  _mm_srli_epi32 // see also _mm512_rol_epi64,
#define SHL32  _mm_slli_epi32
#define SHR16  _mm_srli_epi16
#define SHR64  _mm_srli_epi64
#define SHL64  _mm_slli_epi64

#define SET8  _mm_set1_epi8
#define SETREV8 _mm_setr_epi8
#define SET16  _mm_set1_epi16
#define SET32 _mm_set1_epi32
#define SET64  _mm_set1_epi64x
#define ZERO   _mm_setzero_si128
#define STORE_DOUBLE _mm_store_pd
#define STORE_FLOAT _mm_store_ps
#define EXTRACT16  _mm_extract_epi16

#define ADD8  _mm_add_epi8
#define ADD32  _mm_add_epi32
#define ADD64  _mm_add_epi64
#define MADD16  _mm_madd_epi16 
#define SAD8   _mm_sad_epu8 // _pu8?
#define INT2FLOAT  _mm_cvtepi32_ps
#define INT2DOUBLE _mm_cvtpi32_pd // very expensive

#define BLENDDOUBLE  _mm_blend_pd
#define DUPLICATEDOUBLE  _mm_movedup_pd
//#define MOVEMASK _mm_movemask_ps
//#define BLEND _mm_blend_pd //see also _mm512_mask_inserti64x4_mm_insert_epi64


#if defined SSSE3 // within SSE2
#define SIMD_AVAILABILITY sse2
#define SHUFFLE8  _mm_shuffle_epi8
#else
#define SIMD_AVAILABILITY ssse3
#endif


#elif defined MMX || defined PlainInteger64 // particularly Bit23
#define SIMD_AVAILABILITY no_sse
#define SSEBITS 64U
#define SSEMODE 0U
#define BlockType0 uint64_t
#define BlockType BlockType0
#define UnionType0 uni64

#define AND(B,C)  (B) & (C)
#define  OR(B,C)  (B) | (C)
#define XOR(B,C)  (B) xor (C) 
#define SHR64(B,C)  (B) >> (C)
#define SHR32 SHR64 // unsafe
#define SHR16 SHR64 // unsafe
#define SHL64(B,C)  (B) << (C)
#define SHL32 SHL64 // unsafe
#define SHL16 SHL64 // unsafe
#define SET32 (Ulong) 0x0000000100000001UL * (Ulong)
#define ADD64(B,C)  (B) + (C)
#define ADD32 ADD64 // unsafe
#define ZERO()  0UL
#define LOADU(A) *(A)
#define LOAD LOADU

#define SET8(A) (((BlockType0) (A)) * ((BlockType0) 0x0101010101010101UL))
#if defined MMX
  #define ADD8(B,C)  (BlockType0) _mm_add_pi8((__m64) B, (__m64) C)
#else
  #define ADD8(B,C)  (((BlockType0) (B)) + ((BlockType0) (C))) // unsafe
#endif
#define ANY

  
#else  
#define SIMD_AVAILABILITY no_sse
#define SSEBITS 32U
#define SSEMODE 0U
#define BlockType0 uint32_t
#define BlockType BlockType0
#define UnionType0 uni32
#if defined PlainInteger32
  #define AND(B,C)  (B) & (C)
  #define  OR(B,C)  (B) | (C)
  #define XOR(B,C)  (B) xor (C) 
  #define SHR32(B,C)  (B) >> (C)
  #define SHR16 SHR32 // unsafe
  #define SHL32(B,C)  (B) << (C)
  #define SHL16 SHL32  // unsafe
  #define SHL64 SHL32 // unsafe
  #define SET32 
  #define ADD64(B,C) (B) + (C)
  #define ZERO()  0U
  #define ANY
  #define LOADU(A) *(A)
  #define LOAD LOADU
  #define SET8(A) (((BlockType0) (A)) * ((BlockType0) 0x01010101U))
  #define ADD8(B,C)  (((BlockType0) (B)) + ((BlockType0) (C))) // unsafe !
#else
  #if defined __GNUC__ && defined SCHLATHERS_MACHINE
    #warning No specification of any SIMD.
  #endif
#endif


#endif // AVX512 .. PlaintInteger32


#if defined AVX
#define SCALAR_DEFAULT SCALAR_NEARFMA
#else
#define SCALAR_DEFAULT SCALAR_BASE
#endif



#define BytesPerBlock (SSEBITS / BitsPerByte)
#define ALIGNED __attribute__ ((aligned (BytesPerBlock)))
#define SizeOfDouble 8U
#define SizeOfFloat 4U
#define SizeOfInt 4U
#define doubles (BytesPerBlock / SizeOfDouble)
#define floats (BytesPerBlock / SizeOfFloat)
#define integers (BytesPerBlock / SizeOfInt)
#define BitsPerBlock (BytesPerBlock * BitsPerByte) 



///////////////////////////////////////////////////////////////////////
// checks whether current hardware matches the compilation
//  * mainly intel (and amd) cores
//  * but also GPU
///////////////////////////////////////////////////////////////////////


#define noMISS 0U
#define noUSE 0U
#define anyrelevantUSE 0U
#define gpuUSE 1U
#define avx2USE 2U
#define avxUSE 3U
#define ssse3USE 4U
#define sse2USE 5U
#define avx512fUSE 6U
#define sse41USE 7U

#define USEnMISS 10U
#define gpuMISS 11U
#define avx2MISS 12U
#define avxMISS 13U 
#define ssse3MISS 14U
#define sse2MISS 15U
#define avx512fMISS 16U
#define sse41MISS 17U

#define anyMISS (1 << gpuMISS) | (1 << avx2MISS) | (1 << avxMISS) | \
  (1 << ssse3MISS) | (1 << sse2MISS) | (1 << avx512fMISS)
#define SIMD_INFO					\
  allmiss | alluse | (HAS_PARALLEL || alluse != 0) * (1 << anyrelevantUSE) | \
  ((HAS_PARALLEL || alluse != noUSE) && !(HAS_PARALLEL && allmiss==noMISS)) * \
  (1 << USEnMISS)

#if defined EAX
#if EAX != 0
#define EXX_REDEFINED 1
#endif
#else
#define EAX 0
#endif

#if defined EBX
#if EBX != 1U
#define EXX_REDEFINED 1
#endif
#else
#define EBX 1
#endif

#if defined ECX
#if ECX != 2
#define EXX_REDEFINED 1
#endif
#else
#define ECX 2
#endif

#if defined EDX
#if EDX != 3
#define EXX_REDEFINED 1
#endif
#else
#define EDX 3
#endif


//#define sse3 Available(1, ECX,0)
#define no_sseAvail true
#define no_sseMISS 999U
#define no_sseUSE  999U

#define ssse3Avail Available(1, ECX,9)
#define sse41Avail Available(1, ECX,19)
#define avxAvail Available(1, ECX,28)

#define sseAvail Available(1, EDX,25)
#define sse2Avail Available(1, EDX,26)

#define avx2Avail Available(7, EBX,5)

#define avx512fAvail Available(7, EBX,16)
#define avx512dqAvail Available(7, EBX, 17)
#define avx512pfAvail Available(7, EBX,26)
#define avx512erAvail Available(7, EBX,27)
#define avx512cdAvail Available(7, EBX,28)
#define avx512bwAvail Available(7, EBX,30)
#define avx512vlAvail Available(7, EBX,31)

#define avx512vbmiAvail Available(7, ECX, 1)
#define avx512vmbi2Avail Available(7, ECX, 6)
#define avx512vnniAvail Available(7, ECX, 11)
#define avx512bitalgAvail Available(7, ECX, 12)
#define avx512popcntAvail Available(7, ECX, 14)

#define avx512intersectAvail Available(7, EDX, 8)
#define avx512fp16Avail Available(7, EDX, 23)

#define avx512bf16Avail Available(7, EAX, 5)

// intel Advanced Matrix Calculations
#define amxbf16Avail Available(7, EDX, 22)
#define amxtileAvail Available(7, EDX, 24)
#define amxint8Avail Available(7, EDX, 25)

/*
  PRINTF("blatt %d: %u %u %u %u\n", Blatt, s[0], s[1], s[2], s[3]);	\
      uint32_t a = s[Register];\
      for (int i=31; i>=0; i--){if (i == Bit) PRINTF(" :");PRINTF("%s", (a >> i) & 1 ? "1" : "0");if (i%4 == 0) PRINTF(" ");} PRINTF(" register=%d bit=%d %d: %d %d\n", Register, Bit, bit_SSE, s[Register] & (1 << (Bit)), (s[Register] >> Bit) & 1); \
*/


#define  AVAILABLE_SIMD_OK static inline bool				\
    Available(unsigned VARIABLE_IS_NOT_USED B, int VARIABLE_IS_NOT_USED R, \
	      int VARIABLE_IS_NOT_USED Bit) {	return true; }

#if defined EXX_REDEFINED // unknown system -- don't perform checks
  #define INSTALL_DEFAULT Inone
  #define AVAILABLE_SIMD AVAILABLE_SIMD_OK 
#elif defined ARM32
  #define INSTALL_DEFAULT Iask
  #if defined CROSS_CAPACITY 
    #error "ARM allows only CROSS=noflags and CROSS=FALSE"
  #elif defined REQUIRED_SIMD && REQUIRED_SIMD <= 2
    #error "ARM allows CROSS=noflags and CROSS=FALSE, only."
  #endif
  #define AVAILABLE_SIMD  AVAILABLE_SIMD_OK 
#elif defined __APPLE__ // i.e. apple but isn't arm
  #define INSTALL_DEFAULT Inone
  #if defined CROSS_CAPACITY 
    #error "old MAC-OS allows only CROSS=noflags and CROSS=FALSE"
  #elif defined REQUIRED_SIMD && REQUIRED_SIMD != 3
    #error "old MAC-OS allows CROSS=noflags and CROSS=FALSE, only."
  #endif
  #if defined REQUIRED_SIMD
    #undef REQUIRED_SIMD
  #endif
  #define AVAILABLE_SIMD AVAILABLE_SIMD_OK 
#elif defined WINCPUID
  #define INSTALL_DEFAULT Iask
  #define AVAILABLE_SIMD static inline bool			\
    Available(unsigned Blatt, int Register, int Bit) {		\
      uint32_t s[4];							\
      __cpuid((int *)s, (int) Blatt);					\
      return s[Register] & (1 << (Bit));				\
    }
  }
  #if ! defined MSDOS_WINDOWS
    #error Puzzled about the underlying system. Please contact maintainer.
  #endif
#elif defined LINUXCPUID
  #define INSTALL_DEFAULT Iask
  #define AVAILABLE_SIMD static inline bool				\
    Available(unsigned Blatt, int Register, int Bit) {			\
      uint32_t s[4];							\
      asm volatile							\
      ("cpuid": "=a"(s[0]), "=b"(s[1]),"=c"(s[2]),			\
	"=d"(s[3]):"a"(Blatt),"c"(0));					\
      return s[Register] & (1 << (Bit));				\
      }
#elif defined MINGWCPUID
  #define INSTALL_DEFAULT Iask
// vgl https://github.com/luzexi/MinGW/blob/master/x64/lib/gcc/x86_64-w64-mingw32/4.8.0/include/cpuid.h
  #if defined SCHLATHERS_MACHINE
    #define REACT_ON_DIFFERENT_CPUID_RESULTS				\
         uint32_t u[4];							\
	 asm volatile							\
	 ("cpuid": "=a"(u[0]), "=b"(u[1]),"=c"(u[2]),			\
	  "=d"(u[3]):"a"(Blatt),"c"(0));				\
	 PRINTF("%u %u %u %u\n%u %u %u %u\n%u %u %u %u\n",	       	\
		u[0],u[1],u[2],u[3],					\
		t[0],t[1],t[2],t[3],					\
		s[0],s[1],s[2],s[3]);					\
	 if ((s[0] != t[0] || s[1] != t[1] || s[2] != t[2] || s[3] !=t[3])) BUG 
  #else
    #define REACT_ON_DIFFERENT_CPUID_RESULTS return false
  #endif
  #define AVAILABLE_SIMD static inline bool				\
    Available(unsigned Blatt, int Register, int Bit) {			\
    unsigned int t[4];							\
      if (!__get_cpuid(Blatt, t, t+1, t+2, t+3))	       		\
       ERR1("unallowed cpuid access. %.80s", CONTACT);			\
     unsigned int s[4];							\
     __cpuid(Blatt, s[0], s[1], s[2], s[3]); 				\
    if ((s[0] != t[0] || s[1] != t[1] || s[2] != t[2] || s[3] != t[3])) { \
      /* __get_cpuid does not seem to work for certain registers */	\
      /* indeed results may differ (14 Jan 2022) ! */			\
      REACT_ON_DIFFERENT_CPUID_RESULTS;	 }	      			\
    return s[Register] & (1 << (Bit));		\
   }
#else
  #define INSTALL_DEFAULT Inone
  #define AVAILABLE_SIMD static inline bool				\
      Available(unsigned VARIABLE_IS_NOT_USED B, int VARIABLE_IS_NOT_USED R, \
		int VARIABLE_IS_NOT_USED Bit) {			\
         RFERROR("SIMD checks are not available on your system (on MS systems only under Visual Studio). Use 'CROSS' on Linux systems and alike.");  \
      return false;							\
      }					
  #if defined REQUIRED_SIMD
    #undef REQUIRED_SIMD
  #endif
#endif


#if defined CROSS_CAPACITY
  #if defined REQUIRED_SIMD
   #define ASSERT_TEXT							\
"But the CPU doesn't know about it. As 'CROSS=TRUE' has been chosen as compilation option, it was assumed that the programme was compiled on the most unskilled CPU." // ok
  #else
    #define ASSERT_TEXT \
   "But the CPU doesn't know about it. As 'CROSS' has been chosen as compilation option, it was assumed that each CPU has at least the CROSS skills."
   #endif

#elif defined REQUIRED_SIMD // ! CROSS_CAPACITY
  #if REQUIRED_SIMD == 0 // CROSS = nosimd without -mno-sse2
    #define ASSERT_TEXT \
    "This means 'without SIMD', but the compiler includes SIMD. ('CROSS=nosimd' has been chosen.)"
  #elif REQUIRED_SIMD == 1 // CROSS = nosimd and -mno-sse2
    #define ASSERT_TEXT\
    "This means'without SIMD', but the CPU requires SIMD at a higher level.  Please contact the maintainer."
  #elif REQUIRED_SIMD == 2 // CROSS=NA
    #define ASSERT_TEXT \
    "This means 'without SIMD'), but the compiler includes SIMD (at a higher level). ('CROSS=NA' had been chosen.)"
  #elif REQUIRED_SIMD == 3 // CROSS=F ALSE
    #if defined AVAILABLE_SIMD
      #undef  AVAILABLE_SIMD
    #endif
    #define AVAILABLE_SIMD AVAILABLE_SIMD_OK
    #define ASSERT_TEXT\
    "This situation is unexpected for a PC. Please contact the maintainer."
  #elif REQUIRED_SIMD == 4
   #define ASSERT_TEXT\
    "This situation is unexpected on ARM. Please contact the maintainer."
  #else
    #define ASSERT_TEXT\
    "This leads to an unexpected situation. Please contact the maintainer."
  #endif

#else // ! CROSS_CAPACITY && ! REQUIRED_SIMD 
   #define ASSERT_TEXT\
    "This situation is unexpected for a server. Please contact the maintainer."
  #if defined AVAILABLE_SIMD
    #undef AVAILABLE_SIMD
  #endif
  #define AVAILABLE_SIMD AVAILABLE_SIMD_OK
#endif

    
#define ASSERT_AVAILABILITY(V,W) if ((V##Avail)) {} else {char msg[300]; SPRINTF(msg, "The program was compiled for '%.10s%.5s%.10s. %.200s'", #V, STRCMP(#V, #W) ? " && " : "", STRCMP(#V, #W) ? #W : "", ASSERT_TEXT); RFERROR(msg);}
#define ASSERT_AVAILABILITY_AUX(V,W) ASSERT_AVAILABILITY(V,W) // expands V

#define ASSERT_SIMD(FILE, WHAT)					\
  AVAILABLE_SIMD						\
  Uint check_simd_##FILE() {					\
    ASSERT_AVAILABILITY_AUX(SIMD_AVAILABILITY,WHAT);	 return noMISS;}\
  Uint simd_use_##FILE = WHAT##USE;		         	\
  Uint simd_miss_##FILE = WHAT##MISS

#define ASSERT_SIMD_AUX(FILE, WHAT) ASSERT_SIMD(FILE, WHAT)// expands WHAT	

#define THIS_FILE_ANYSIMD ASSERT_SIMD_AUX(this_file, SIMD_AVAILABILITY)

#define SIMD_MISS(FILE, WHAT)				        \
  Uint check_simd_##FILE() { return 1<<WHAT##MISS; }		\
  Uint simd_use_##FILE = WHAT##USE;		         	\
  Uint simd_miss_##FILE = WHAT##MISS		               

#define EXTERN_SIMD_CHECK(FILE)			\
  extern Uint simd_use_##FILE;			\
  extern Uint check_simd_##FILE()

#define CHECK_FILE(FILE)				\
  allmiss |= (miss = check_simd_##FILE());		\
  if (miss) {} else alluse |= 1 << simd_use_##FILE

#define CHECK_THIS_FILE_ANYSIMD				\
  Uint alluse = noUSE,					\
    miss = noMISS,					\
    allmiss = noMISS;					\
  CHECK_FILE(this_file)

#endif


/*  LocalWords:  endif
 */
