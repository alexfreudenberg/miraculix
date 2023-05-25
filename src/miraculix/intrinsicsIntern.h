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

#ifndef miraculix_inl_H
#define miraculix_inl_H 1


/*

in particular for scalar of 1bit and 2bit

 */



#define DEF_POPCNT_EMUL64						\
  static const unsigned char PopCountTable0[16] =			\
    {/*    0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F */		\
      /**/ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};		\
  static unsigned int PopCountTable1[256 + 0] = {255};			\
  static inline BlockType						\
  POPCNT_EMUL64(BlockType0 v) {						\
    if (PopCountTable1[0] == 255) {					\
      for (unsigned short i=0; i<256; i++) {				\
	unsigned short lo = i % 16,					\
	  hi = i / 16;							\
	PopCountTable1[i] = PopCountTable0[lo] + PopCountTable0[hi];	\
      }									\
    }									\
    BlockType sum = 0;							\
    unsigned char *p = (unsigned char*) &v;				\
    assert(sizeof(BlockType0) == BytesPerBlock);			\
    for (unsigned short i=0; i<BytesPerBlock; i++) {			\
      sum += (BlockType0) PopCountTable1[MAX(0, MIN(p[i], 255))  ];	\
    }									\
  return sum;								\
  }


const BlockType F1O1 = SET32(0xF0F0F0F0);
const BlockType O1F1 = SET32(0x0F0F0F0F);
const BlockType N1E1 = SET32(0x55555555);
const BlockType E1N1 = SET32(0xAAAAAAAA);
const BlockType N2E2 = SET32(0x33333333);
const BlockType E2N2 = SET32(0xCCCCCCCC); 


#if (defined AVX512 && ! (defined AVX512VPOPCNTDQ || defined AVX512VL)) || defined AVX2 || defined SSSE3

const BlockType PopCountTable = SETREV8(//   0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F
			  /**/ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
static inline BlockType
POPCNT_EMUL64(BlockType0 v) { // sizeof(Long)=64
  const BlockType vlo = AND(v, O1F1);
  const BlockType vhi = AND(SHR32(v, 4), O1F1);
  const BlockType plo = SHUFFLE8(PopCountTable, vlo);
  const BlockType phi = SHUFFLE8(PopCountTable, vhi);
  const BlockType sum = ADD8(plo, phi);
  const BlockType res = SAD8(sum, ZERO());
  return res;
}
#endif


#if defined AVX512

static inline Ulong HORIZ_ADD64(BlockType0 A) {
 __m256i z = _mm256_add_epi64(_mm512_extracti64x4_epi64(A, 1),
			      _mm512_extracti64x4_epi64(A, 0));
  uni128 x;
  x.vi = _mm_add_epi64(_mm256_extracti128_si256(z, 1),
		       _mm256_extracti128_si256(z, 0));
  return (Ulong) x.u64[0] + (Ulong) x.u64[1];
}


#elif defined AVX2

static inline Ulong HORIZ_ADD64(BlockType0 A) {		      
  uni128 x;
  x.vi = _mm_add_epi64(_mm256_extracti128_si256(A, 1),
		       _mm256_extracti128_si256(A, 0));
  return (Ulong) x.u64[0] + (Ulong) x.u64[1];
}


#elif defined SSE2
static inline Ulong HORIZ_ADD64(BlockType0 A) {
  UnionType x; // OK
  x.vi = A;
  return (Ulong) x.u64[0] + (Ulong) x.u64[1];
}

#if ! defined SSSE3
DEF_POPCNT_EMUL64
#endif



#else

#define HORIZ_ADD64
DEF_POPCNT_EMUL64


#endif


#endif
