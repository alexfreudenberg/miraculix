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


#ifndef miraculix_initrinsics_H
#define miraculix_initrinsics_H 1



#ifndef UINT64_C
#define UINT64_C(X) X
#endif


#define BytesPerUnit 4 // OK
#define BitsPerUnit (BytesPerUnit * BitsPerByte) // OK
#define UnitsPerBlock (BytesPerBlock / BytesPerUnit) 


#if defined BytesPerBlock
#if ! defined BytesPerPartUnit
  #define BytesPerPartUnit BytesPerUnit
#endif
#define BitsPerPartUnit (BitsPerByte * BytesPerPartUnit)
#define PartUnitsPerBlock (BytesPerBlock / BytesPerPartUnit)

#define PartUnitsPerUnit (BytesPerUnit / BytesPerPartUnit)
#if defined PartUnitsPerUnit && PartUnitsPerUnit < 1
#undef PartUnitsPerUnit
#endif

union block_compressed{
  uint16_t b[PartUnitsPerBlock]; // OK
  // uint32_t u32[2]; falsch
  BlockType x;
};
#endif


#if defined BitsPerCode  /// real valued in general!!

  // none of the divisions here are truncated, except, maybe, the following:
#define CodesPerPartUnit ((Uint)((BitsPerPartUnit / BitsPerCode)))//truncated in general
#define codingBitsPerPartUnit ((Uint) (CodesPerPartUnit * BitsPerCode + 0.99))
#define two_codingBitsPerPartUnit (1UL << codingBitsPerPartUnit)
#define deltaBitsPartUnit (BitsPerPartUnit - codingBitsPerPartUnit)
#define CodesPerBlock (PartUnitsPerBlock * CodesPerPartUnit)
#define CodesPerByte ((Uint)(BitsPerByte / BitsPerCode + 0.1))
#define CodeMask ((1U << ((Uint) BitsPerCode)) - 1U)

#if defined PartUnitsPerUnit
#define CodesPerUnit (PartUnitsPerUnit * CodesPerPartUnit)
#endif

#endif // defined BitsPerCode


//#define ALIGNED __attribute__ (aligned BytesPerBlock)))
#define UnionType UnionType0 ALIGNED // OK


// #if defined AVX || defined SSE2 || defined AVX2 || defined AVX512
#if defined SSE2 || defined AVX2
#define VI .vi
#define U128(X) .u128[X]
#endif


// needed for immitations of sse/ssse3 functions
#define X128F64none(Z,A) Z(A,0); Z(A,1); 
#define X128F64(Z,A,B) Z(A,B,0); Z(A,B,1); 
#define X128F32(Z,A,B) X128F64(Z,A,B); Z(A,B,2); Z(A,B,3); 
#define X128F16(Z,A,B) X128F32(Z,A,B); Z(A,B,4); Z(A,B,5); Z(A,B,6); Z(A,B,7);
#define X128F8(Z,A,B) X128F16(Z,A,B); Z(A,B,8); Z(A,B,9);Z(A,B,10); Z(A,B,11);Z(A,B,12); Z(A,B,13); Z(A,B,14); Z(A,B,15);
#define X128F64bi(Z,A,B,C) Z(A,B,C,0); Z(A,B,C,1); 
#define X128F32bi(Z,A,B,C) X128F64bi(Z,A,B,C); Z(A,B,C,2); Z(A,B,C,3);
#define X128F16bi(Z,A,B,C) X128F32bi(Z,A,B,C); Z(A,B,C,4); Z(A,B,C,5); Z(A,B,C,6); Z(A,B,C,7);
#define X128F8bi(Z,A,B,C) X128F16bi(Z,A,B,C); Z(A,B,C,8); Z(A,B,C,9);Z(A,B,C,10); Z(A,B,C,11);Z(A,B,C,12); Z(A,B,C,13); Z(A,B,C,14); Z(A,B,C,15);

/*
const BlockType F1O1 = SET32(0xF0F0F0F0);
const BlockType O1F1 = SET32(0x0F0F0F0F);

const BlockType N1E1 = SET32(0x55555555);
const BlockType E1N1 = SET32(0xAAAAAAAA);
const BlockType N2E2 = SET32(0x33333333);
const BlockType E2N2 = SET32(0xCCCCCCCC); 
*/


#if defined AVX512

#define DOUBLE_CROSSLINE_PERMUTE _mm512_permutex2var_epi32
#define CROSSLINE_PERMUTE
#if defined AVX512VPOPCNTDQ || defined AVX512VL
  #define POPCNT64 _mm512_popcnt_epi64
#else
  #define POPCNT64 POPCNT_EMUL64
#endif
#define STORE _mm512_store_si512
#define STOREU _mm512_storeu_si512

#elif defined AVX2
#define INT2FLOAT  _mm_cvtepi32_ps // 128 bit only !!
#define INT2DOUBLE(B)  _mm_cvtpi32_pd(B) // 128 bit only !! latency 4 thr 1
#define CROSSLINE_PERMUTE _mm256_permutevar8x32_epi32
#define HALFREGISTER_PERMUTE _mm256_permute2x128_si256
#define POPCNT64 POPCNT_EMUL64
#define SHUFFLE32 _mm256_shuffle_epi32
#define GREATER_THAN _mm256_cmpgt_epi32
#define STORE _mm256_store_si256 
#define STOREU _mm256_storeu_si256 


#elif defined SSE2
#if defined SSSE3 // within SSE2
#else
  #define SHUFFLE8X(A,B,C,X) ((uni128*) &A)->u8[X] = ((uni128*) &B)->u8[((uni128*) &C)->u8[X]] // ACHTUNG A != B !!
  #define SHUFFLE8(A,B,C) { X128F8bi(SHUFFLE8X,A,B,C) }
#endif

#define SHUFFLE32 _mm_shuffle_epi32

#define CROSSLINE_PERMUTE _mm_shuffle_epi32
#define HALFREGISTER_PERMUTE(A,B,C)\
  (__m128i) _mm_shuffle_pd((__m128d) A, (__m128d) B, C)
#define POPCNT64 POPCNT_EMUL64
#define GREATER_THAN _mm_cmpgt_epi32
#define STORE _mm_store_si128 
#define STOREU _mm_storeu_si128 

#else

#define POPCNT64 POPCNT_EMUL64
#define STORE(A,B) *(A)=B
#define STOREU STORE

#if defined PlainInteger64
  #define VI
  #define CROSSLINE_PERMUTE(A, B) (((A) >> ((B) << 5)) AND 0x00000000FFFFFFFF)

#elif defined PlainInteger32
  #define VI
  #define CROSSLINE_PERMUTE(A, B) (((A) >> ((B) << 5)) AND 0x0000FFFF)
#endif

#endif // AVX512 ... PlainInteger


#define AND4(To, From, From2)				\
  const BlockType a##To = AND(a##From, a##From2);		\
  const BlockType b##To = AND(b##From, b##From2);		\
  const BlockType c##To = AND(c##From, c##From2);		\
  const BlockType d##To = AND(d##From, d##From2);
#define OR4(To, From, From2)				\
  const BlockType a##To = OR(a##From, a##From2);		\
  const BlockType b##To = OR(b##From, b##From2);		\
  const BlockType c##To = OR(c##From, c##From2);		\
  const BlockType d##To = OR(d##From, d##From2);


#define AND4_C(To, From, N)				\
  const BlockType a##To = AND(a##From, N);		\
  const BlockType b##To = AND(b##From, N);		\
  const BlockType c##To = AND(c##From, N);		\
  const BlockType d##To = AND(d##From, N);
#define SHR32_4(To, From, N)				\
  const BlockType a##To = SHR32(a##From, N);		\
  const BlockType b##To = SHR32(b##From, N);		\
  const BlockType c##To = SHR32(c##From, N);		\
  const BlockType d##To = SHR32(d##From, N);
#define SHL32_4(To, From, N)				\
  const BlockType a##To = SHL32(a##From, N);		\
  const BlockType b##To = SHL32(b##From, N);		\
  const BlockType c##To = SHL32(c##From, N);		\
  const BlockType d##To = SHL32(d##From, N);
#define SHR64_4(To, From, N)			\
  const BlockType a##To = SHR64(a##From, N);		\
  const BlockType b##To = SHR64(b##From, N);		\
  const BlockType c##To = SHR64(c##From, N);		\
  const BlockType d##To = SHR64(d##From, N);
#define PERMUTE4(To, From, N)					\
  const BlockType a##To = CROSSLINE_PERMUTE(a##From, N);	\
  const BlockType b##To = CROSSLINE_PERMUTE(b##From, N);	\
  const BlockType c##To = CROSSLINE_PERMUTE(c##From, N);	\
  const BlockType d##To = CROSSLINE_PERMUTE(d##From, N);
#define SHUFFLE8_4(To, N, From)				\
  const BlockType a##To = SHUFFLE8(N, a##From);		\
  const BlockType b##To = SHUFFLE8(N, b##From);		\
  const BlockType c##To = SHUFFLE8(N, c##From);		\
  const BlockType d##To = SHUFFLE8(N, d##From);
#define SAD8_4(To, From, N)			\
  const BlockType a##To = SAD8(a##From, N);		\
  const BlockType b##To = SAD8(b##From, N);		\
  const BlockType c##To = SAD8(c##From, N);		\
  const BlockType d##To = SAD8(d##From, N);


#define ADDUP8_4(Sum, From)				\
  a##Sum = ADD8(a##Sum, a##From);			\
  b##Sum = ADD8(b##Sum, b##From);			\
  c##Sum = ADD8(c##Sum, c##From);			\
  d##Sum = ADD8(d##Sum, d##From);
#define HORIZ_ADD64_4(From)				\
  a += HORIZ_ADD64(a##From);		\
  b += HORIZ_ADD64(b##From);		\
  c += HORIZ_ADD64(c##From);		\
  d += HORIZ_ADD64(d##From);



// #define Uint uint32 // Window




#endif




