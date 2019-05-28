
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




/*

__m128d _mm_blend_pd(__m128d v1, __m128d v2, const int mask) 



_mm256_movemask_ps


sse.m128_f32[0], sse.m128_f32[1], sse.m128_f32[2], sse.m128_f32[3]

typedef union __declspec(intrin_type) _CRT_ALIGN(16) __m128 {
float m128_f32[4];
unsigned __uint64 m128_u64[2];
__int8 m128_i8[16];
__int16 m128_i16[8];
__int32 m128_i32[4];
__uint64 m128_i64[2];
unsigned __int8 m128_u8[16];
unsigned __int16 m128_u16[8];
unsigned __int32 m128_u32[4];
} __m128;

 */

#include <smmintrin.h>
#include <inttypes.h> // uintptr_t
#include "intrinsics.h"
#include "kleinkram.h"
#include "miraculix.h"
#include <General_utils.h>
#include "xport_import.h"
#include "AutoMiraculix.h"
#include "options.h"
#include "general.relmatrix.h"

//#define DO_FLOAT 1

// #define DO_PARALLEL 1 // done by basic_utils.h


bool debugging = false;
SEXP Debug() { debugging = true; return R_NilValue; }
SEXP StopDebug() { debugging = false;  return R_NilValue; }

//#define show printf
//
#define show(X)
// void show(const char *s) {  if (debugging) PRINTF("%.50s", s); }


#define checking true


 
// #define Uint uint32 // Window

// #define INLINE
#define INLINE inline

#if defined AVX2
#define BlockType0 __m256i
#define ALIGNED __attribute__ ((aligned (32)))
#define BlockType __m256i ALIGNED
#define BytesPerBlock 32
#define Double __m256d
#define Float __m256
#define AND _mm256_and_si256
#define OR  _mm256_or_si256
//#define XOR(X, Y) X xor Y // _mm256_xor_si256
#define XOR _mm256_xor_si256
#define SHR _mm256_srli_epi32
#define SHL _mm256_slli_epi32
#define ZERO _mm256_setzero_si256()
#define ZEROFLOAT _mm256_setzero_ps()

BlockType set_epi32(Uint X) {
  union { BlockType b; __m128i s[2]; } Y;
  Y.s[0] = Y.s[1] = _mm_set1_epi32(X);
  return(Y.b);
}
#define SET_EPI32(X) set_epi32(X) // _mm256_set1_epi32(X) funktioniert nicht!?
//#define SET_EPI32(X) _mm256_set_epi32(X,X,X,X,X,X,X,X)
#define ADD32 _mm256_add_epi32
#define INT2FLOAT _mm256_cvtepi32_ps
#define MULTFLOAT _mm256_mul_ps 
#define ADDFLOAT  _mm256_add_ps 
#define SUBFLOAT  _mm256_sub_ps 
#define ZEROFLOAT _mm256_setzero_ps()
#define INT2DOUBLE _mm256_cvtepi32_pd  // very expensive
#define MULTDOUBLE _mm256_mul_pd 
#define ADDDOUBLE  _mm256_add_pd 
#define SUBDOUBLE  _mm256_sub_pd
#define ZERODOUBLE _mm256_setzero_pd()
#define LOAD _mm256_load_si256
#define LOAD_DOUBLE _mm256_load_pd
#define LOAD1_DOUBLE _mm256_load1_pd
#define STORE_DOUBLE _mm256_store_pd
#define SHUFFLE8 _mm256_shuffle_epi8 //https://software.intel.com/en-us/node/524017
#define ADD8 _mm256_add_epi8 //https://software.intel.com/en-us/node/523879
#define SAD8 _mm256_sad_epu8
#define HalfBlockType __m128i
#define MOVEMASK _mm256_movemask_ps

#elif defined SSE2
#define BlockType0 __m128i 
#define ALIGNED __attribute__ ((aligned (16)))
#define BlockType __m128i ALIGNED
#define BytesPerBlock 16
#define Double __m128d
#define Float __m128
#define AND _mm_and_si128
#define OR  _mm_or_si128
#define XOR  _mm_xor_si128
#define SHR _mm_srli_epi32 // see also _mm512_rol_epi64,
#define SHL _mm_slli_epi32
#define ZERO _mm_setzero_si128()
#define ZEROFLOAT _mm_setzero_ps()
#define SET_EPI32(X) _mm_set1_epi32(X)
#define SET_EPI64 _mm_set1_epi64
//#define SET_EPI32(X) _mm_set_epi32(X,X,X,X)
#define ADD32 _mm_add_epi32
#define INT2FLOAT _mm_cvtepi32_ps
#define MULTFLOAT _mm_mul_ps 
#define ADDFLOAT  _mm_add_ps 
#define SUBFLOAT  _mm_sub_ps 
#define ZEROFLOAT _mm_setzero_ps()
#define INT2DOUBLE _mm_cvtpi32_pd // very expensive
#define MULTDOUBLE _mm_mul_pd 
#define ADDDOUBLE  _mm_add_pd 
#define SUBDOUBLE  _mm_sub_pd
#define ZERODOUBLE _mm_setzero_pd()
#define LOAD _mm_load_si128
#define LOADU _mm_loadu_si128
#define LOAD_DOUBLE _mm_load_pd
#define LOAD1_DOUBLE _mm_load1_pd
#define STORE_DOUBLE _mm_store_pd
#define SHUFFLE8 _mm_shuffle_epi8
#define ADD8 _mm_add_epi8
#define SAD8 _mm_sad_epu8
#define HalfBlockType __m64
#define MOVEMASK _mm_movemask_ps
#define BLEND _mm_blend_pd //see also _mm512_mask_inserti64x4_mm_insert_epi64


#else
#define BlockType Programm_nneds_at_least_SSE_commands
#endif



#ifdef DO_FLOAT
#define real float
#define Real Float
#define ZEROREAL ZEROFLOAT 
#define INT2REAL INT2FLOAT 
#define MULTREAL MULTFLOAT 
#define ADDREAL ADDFLOAT  
#define SUBREAL SUBFLOAT  
#else
#define real double
#define Real Double
#define ZEROREAL ZERODOUBLE 
#define INT2REAL INT2DOUBLE 
#define MULTREAL MULTDOUBLE 
#define ADDREAL ADDDOUBLE  
#define SUBREAL SUBDOUBLE  
#endif

					   
#define BytesPerDouble 8L
#define BitsPerByte 8L
#define BitsPerCode 2L 
#define BitsPerGenoCode 2L
#define BytesPerUnit 4L
#define UnitsPerBlock (BytesPerBlock / BytesPerUnit)
#define BitsPerUnit (BytesPerUnit * BitsPerByte)
#define BitsPerBlock (BytesPerBlock * BitsPerByte)
#define CodesPerUnit (BitsPerUnit / BitsPerCode)
#define CodesPerBlock (CodesPerUnit * UnitsPerBlock)
#define CodesPer32 (32 / BitsPerGenoCode)
#define CodesPerByte (BitsPerByte / BitsPerCode)
#define DoublesPerBlock (BytesPerBlock / BytesPerDouble)
#define Bytes128 16

// #define Blocks(X) (1L + (X - 1L) / CodesPerBlock)
#define algn_general(X)  ((1L + (uintptr_t) (((uintptr_t) X - 1L) / BytesPerBlock)) * BytesPerBlock)
//#define UnitsAlign1(S) ( (1L + Blocks(S)) * UnitsPerBlock - 1L)
//#define UnitsAlign(S,I) (UnitsAlign1(S) + Blocks(S) * (I-1) * UnitsPerBlock)
//#define RealAlign(X) (BytesPerBlock / sizeof(real) + (1L +  (X - 1L) / CodesPerBlock) * CodesPerBlock)

Uint inline *algn(int *X) {assert(algn_general(X)>=(uintptr_t)X); return (Uint *) algn_general(X); }
real inline *algnrl(real *X) {assert(algn_general(X)>=(uintptr_t)X); return (real *)  algn_general(X); }
Uint inline Blocks(Uint X) { return 1L + (X - 1L) / CodesPerBlock; }
Uint inline UnitsAlign1(Uint S) { return (1L + Blocks(S)) * UnitsPerBlock - 1L; }
Uint inline UnitsAlign(Uint S, Uint I) {
  return UnitsAlign1(S) + Blocks(S) * (I-1) * UnitsPerBlock;
}
Uint inline RealAlign(Uint X) {
  return BytesPerBlock / sizeof(real) + (1L +  (X - 1L) / CodesPerBlock) * CodesPerBlock;
}


Uint zaehler = 0;
#define WHAT_NAMES {"haplo vector", "geno vector", "geno matrix"}

#define BITS_GENE_INPUT 6
#define BITS_SEX 1
#define BITS_INDIVIDUALS 22
#define BITS_HAPLO 3
//////// total sum 32 bit

#define MAX_GENE_INPUT (1 << BITS_GENE_INPUT)
#define MAX_SEX (1 << BITS_SEX)
#define MAX_INDIVIDUALS (1 << BITS_INDIVIDUALS)
#define MAX_HAPLO (1 << BITS_HAPLO)

#define PATTERN_GENE_INPUT  (MAX_GENE_INPUT - 1)
#define PATTERN_SEX (MAX_SEX - 1)
#define PATTERN_INDIVIDUALS (MAX_INDIVIDUALS - 1)
#define PATTERN_HAPLO (MAX_HAPLO - 1)



//93843

// [1 + (X - 1) / UnitsPerBlock] * UnitsPerBlock + (UnitsPerBlock - 1)
// last part: to make finally sure that first used element of R-vector
// is aligned (maximum distance to 0-elemtn: UnitsPerBlock - 1)
// first part: then memory is read in blocks; so next number of Bytes that is
// is divisible by Unitsperblock * Bytesperunit is taken:
// 1 + (X - 1) / UnitsPerBlock]



					   //                 0  1  2  3  4  5  6  7   8  9 10 11 12 13 14 15
#define andshuffle 0, 1, 4, 0, 1, 2, 5, 0,  4, 5, 8, 0, 0, 0, 0, 0
#define xorshuffle 0, 0, 0, 2, 0, 0, 0, 2,  0, 0, 0, 2, 2, 2, 2, 4

Uint permut_finv_g_h[CodesPerBlock];

BlockType and2scalar, xor2scalar, MASK1, MASK2,
      LH, BitMaskStart[CodesPerBlock], BitMaskEnd[CodesPerBlock];
bool do_warn_change = true, do_warn_changeOK = true;

#define ALIGN_HAPLO 0
#define ALIGN_HAPLOGENO 1
#define ALIGN_VECTOR 2
#define ALIGN_RELATION 3
#define ALIGN_ALLELE 4
#define ALIGN_COMPUTE 5
#define ALIGN_CODED 6
#define ALIGN_LAST ALIGN_CODED
bool PRESERVE = false;
int *Memory[ALIGN_LAST + 1] = { NULL };
Uint
     nMem[ALIGN_LAST + 1] = { 0 };
real *Pd0 = NULL; 					   
 
typedef union addr_Uint{
  uintptr_t a;
  Uint u[sizeof(uintptr_t) / BytesPerUnit];
} addr_Uint;

typedef union {
    BlockType b;
    Uint u[UnitsPerBlock];
  } BlockUnitType ALIGNED;



Uint ToInti(SEXP X, int i) {
  switch(TYPEOF(X)) {
  case INTSXP : return (Uint) INTEGER(X)[i];
  case LGLSXP : return (Uint) LOGICAL(X)[i];
  case REALSXP : return (Uint) REAL(X)[i];
  default : ERR("not of numerical type (in ToInti)");
  }
}
#define ToInt0(X) ToInti(X, 0)

SEXP dolocking(SEXP Do) {
 if (Do != R_NilValue) {
   PRESERVE = LOGICAL(Do)[0];
   if (PRESERVE) BUG;
   for (Uint i=0; i<= ALIGN_LAST; i++) {
     FREE(Memory[i]);
     nMem[i] = 0;
   }
 }

 SEXP Ans;
 PROTECT(Ans = allocVector(LGLSXP, 1));
 LOGICAL(Ans)[0] = PRESERVE;
 UNPROTECT(1);
 return Ans;
}

void showblock(Uint *X, Uint snps) { //
  Uint *x = X,
    bits = snps * BitsPerCode,
    endfor = 1 + (Uint) ((bits - 1) / BitsPerBlock);
  
  for (Uint d=0; d<endfor; d++, x+=UnitsPerBlock) {
    for (Uint b=0; b<UnitsPerBlock; b++, x++) {
      Uint pattern = (Uint) 0x80000000;
      for (Uint s=0; s<BitsPerUnit; s++) {
        PRINTF("%d", (*x & pattern) != 0);
	pattern >>= 1;
	if (s==15) { PRINTF("."); }
      }
      PRINTF(" ");
    }
  }

  PRINTF("\n");
}

void showblock(BlockType0 X) { showblock((Uint *) &X, CodesPerBlock); }
void showblock(BlockType0 X, int blocks) {
  showblock((Uint *) &X, blocks * CodesPerBlock);
}

void showblock(Uint *X, Uint snps, Uint individuals) { 
  for (Uint i=0; i<individuals; i++) {
#ifdef SCHLATHERS_MACHINE    
    PRINTF("%ld\n", (uintptr_t) (X + i * Blocks(snps) * UnitsPerBlock));
#endif 
    showblock(X + i * Blocks(snps) * UnitsPerBlock, snps);
  }
}


SEXP Information = R_NilValue,
  Method = R_NilValue;
static void Init() {
  assert(BytesPerBlock == sizeof(BlockType0));
  assert(BytesPerUnit == sizeof(Uint));
  assert(sizeof(int) == sizeof(Uint));
  PRINTF("The following basic compiler options have been set:\nparallel computing: %.50s\nfloating point single precision: %.50s\nSIMD: %.50s\n",
#ifdef DO_PARALLEL
	 "yes",
#else
	 "no",
#endif
#ifdef DO_FLOAT
	 "yes",
#else
	 "no",
#endif	 
#ifdef AVX
	 "AVX"
#elif SSSE3
	 "SSSE3"
#elif SSE2
	 "SSE2"
#else
	 "none"
#endif
	 );


  if (sizeof(uintptr_t) > 4 * BytesPerUnit) ERR("address cannot be stored. Please contract author.");

  Information = install("information");
  Method = install("method");
 
  BlockType Einsen = SET_EPI32(0xFFFFFFFF);
  BlockUnitType *BME = (BlockUnitType*) BitMaskEnd;
  BME->u[0] = 0x03;

  for (Uint k=1; k<CodesPerBlock; k++) {
    BitMaskStart[k] = XOR(BitMaskEnd[k-1], Einsen); // not
    BitMaskEnd[k] = BitMaskEnd[k-1];
    BME = (BlockUnitType*) BitMaskEnd + k;
    BME->u[k % UnitsPerBlock] |= 0x03 << (BitsPerCode * (k / UnitsPerBlock));
  }
  BitMaskStart[0] = Einsen;
   
  __m128i *as = (__m128i *) &and2scalar,
    *xs = (__m128i *) &xor2scalar;
  Uint endfor = BytesPerBlock / sizeof(__m128i);
  for (Uint i=0; i<endfor; i++) {
    _mm_setr_epi8(andshuffle);
    as[i] = _mm_setr_epi8(andshuffle);
    xs[i] = _mm_setr_epi8(xorshuffle);
  }	 

  LH = SET_EPI32(0x0F0F0F0F);
  MASK1 = SET_EPI32(0x01010101);
  MASK2 = SET_EPI32(0x02020202);

  // see permutation.tex !
  Uint f[CodesPerBlock], g[CodesPerBlock], h[CodesPerBlock],
    finv[CodesPerBlock];
    
  for (Uint i=0, z=0; i<CodesPerUnit; i++) {
    for (Uint j=0; j<UnitsPerBlock; j++) {
      f[z++] = j * CodesPerUnit + i;
    }
  }
  Ext_orderingInt((int *) f, CodesPerBlock, 1, (int *) finv);

  for (Uint i=0, z=0; i<CodesPerByte; i++) {
    for (Uint j=0; j<BytesPerBlock; j++) {
      g[z++] = j * CodesPerByte + i;
   }
  }
  
  for (Uint k=0, z=0; k<CodesPerByte; k++) {
    for (Uint j=0; j<BytesPerUnit; j++) {
      for (Uint i=0; i<UnitsPerBlock; i++) {
	h[z++] = i * BytesPerUnit + j + k * BytesPerBlock;
     }
    }
  }

  for (Uint k=0; k<CodesPerBlock; k++) {
    permut_finv_g_h[k] = finv[g[h[k]]];
  }
  
}

SEXP init_parallelcompute() { Init(); return  R_NilValue; }



#define ORIGIN_GENERATION 0
#define ORIGIN_SEX 1
#define ORIGIN_NR 2
#define ORIGIN_HAPLO 3

void INLINE decodeOrigin(Uint x, Uint *ans) {
  ans[ORIGIN_HAPLO] = (x & PATTERN_HAPLO);
  x >>= BITS_HAPLO;
  ans[ORIGIN_NR] = (x & PATTERN_INDIVIDUALS);
  x >>= BITS_INDIVIDUALS;  
  ans[ORIGIN_SEX] = (x & PATTERN_SEX);
  x >>= BITS_SEX;
  ans[ORIGIN_GENERATION] = x;
}

SEXP decodeOrigins(SEXP M, SEXP Line) {  
  if (Information == R_NilValue) Init();
  Uint line = TYPEOF(Line) == INTSXP ? INTEGER(Line)[0] : (int) REAL(Line)[0];
  SEXP Ans;
  PROTECT(Ans = allocVector(INTSXP, 4));
  Uint *ans = (Uint*) INTEGER(Ans);
  decodeOrigin(ToInti(M, line - 1), ans);
  for (Uint i=0; i<4; ans[i++]++);
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return Ans;
}



SEXP codeOrigins(SEXP M) {
  if (Information == R_NilValue) Init();
  Uint ncol = ncols(M),
    nrow = nrows(M);
  assert(ncol == 4);
  
  SEXP Ans;
  PROTECT(Ans = allocVector(INTSXP, nrow));
 
  Uint *ans = (Uint*) INTEGER(Ans);
  for (Uint i=0; i<nrow; i++) ans[i] = 0;  // unnoetig?!

#define DO1(TYPE, RTYPE)						\
    TYPE *gen = RTYPE(M),						\
      *sex = gen + nrow,						\
      *nr = sex + nrow,							\
      *haplo = nr + nrow					       
    
#define DO2								\
   for (Uint i=0; i<nrow; i++) {					\
      Uint g = (Uint) gen[i],					\
	s = (Uint) sex[i],					\
	n = (Uint) nr[i],					\
	H = (Uint) haplo[i]; 					\
      if (g < 1 || g > MAX_GENE_INPUT|| s < 1 || s > MAX_SEX ||		\
	  n < 1 || n > MAX_INDIVIDUALS ||	H < 1 || H > MAX_HAPLO) \
	ERR1("some value in row %d out of bound", i + 1);		\
      ans[i] = ((((((g - 1) << BITS_SEX) +				\
		   s - 1) << BITS_INDIVIDUALS) +			\
		 n - 1) << BITS_HAPLO) +				\
		   H - 1;						\
      									\
      Uint control[4];							\
      decodeOrigin(ans[i], control);					\
      if (control[0]+1 != g || control[1]+1 != s || control[2]+1 != n || \
	  control[3]+1 != H) {						\
	/*	p rintf("i=%d: generation:%d->%d sex:%d->%d nr:%d->%d haplo:%d->%d\n", i, g, control[0], s, control[1], n, control[2], H, control[3]); */ \
	BUG;								\
      }									\
   }

  
  if (TYPEOF(M) == REALSXP) {
    DO1(double, REAL);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
    DO2   
  } else if (TYPEOF(M) == INTSXP) {
    DO1(int, INTEGER);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)  
#endif
    DO2   
  } else BUG;
  
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return Ans;
}

static bool warnold = true;
Uint *Align(SEXP CM, Uint nr) {
  Uint *info;
  SEXP Infos = getAttrib(CM, Information);
  if (TYPEOF(Infos) == INTSXP) {
    if (warnold) {
      warnold = false;
      PRINTF("warning: stored data from obsolete miraculix coding used.\n");
    }
    info = (Uint *) INTEGER(Infos);
  } else {
    assert(Infos != R_NilValue);
    info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  }
  uintptr_t address = (uintptr_t) INTEGER(CM);
  Uint *algnaddress = algn(INTEGER(CM));
  addr_Uint au;
  for (Uint i=0; i<sizeof(uintptr_t)/BytesPerUnit; i++)
    au.u[i] = info[ADDR0 + i];
  if (au.a == address) {
    assert((uintptr_t) algnaddress % BytesPerBlock == 0);
    return algnaddress;
  } else if (au.a % BytesPerUnit ==  address % BytesPerUnit) {
#ifdef SCHLATHERS_MACHINE
  if (do_warn_changeOK) {
    PRINTF("Address has changed in a coded object (%d) by 'gc' or disk saving. But luckily got same modulus (%lu %lu)).\n", info[WHAT], (uintptr_t) au.a % (1024 * 1024), address % (1024 * 1024));
    do_warn_changeOK = debugging;
  }
#endif    
    assert((uintptr_t) algnaddress % BytesPerBlock == 0);
    return algnaddress;
  }
#ifdef SCHLATHERS_MACHINE
  if (do_warn_change) {
    PRINTF("Address has changed in a coded object (%d) by 'gc' or disk saving\n", info[WHAT]);
    do_warn_change = debugging;
  }
#endif    
  
  if (PRESERVE) ERR("severe error: has 'dolocking' be called in meanwhile? If 'dolocking' has not been called at all or only at the very beginning, please contaact maintainer.");
  if (au.a % BytesPerBlock == address % BytesPerBlock) return algnaddress;
  if (nMem[nr] < info[MEM]) {
    FREE(Memory[nr]);
    nMem[nr] = info[MEM];
    // printf("s/i=%d, %d; mem=%10g %d\n", info[SNPS], info[INDIVIDUALS], nMem[nr] / 1e9, BytesPerUnit);
    Memory[nr] = (int*) CALLOC(nMem[nr], BytesPerUnit);
  }
  Uint *algnMem = algn(Memory[nr]),
    bytes = Blocks(info[SNPS]) * info[INDIVIDUALS] * BytesPerUnit;
  assert((uintptr_t) algnMem % BytesPerBlock == 0);
  assert((uintptr_t)algnMem + bytes < (uintptr_t) Memory[nr] + nMem[nr]);
  assert(bytes < (Uint) length(CM) * BytesPerUnit);
  assert(bytes < nMem[nr] * BytesPerUnit);
  
  MEMCOPY(algnMem, algnaddress, bytes);
  return algnMem;
}


SEXP static create_codevector(Uint what, Uint snps, Uint individuals) {
  Uint mem = UnitsAlign(snps, individuals);  
  SEXP Code =  CreateEmptyCodeVector(snps, individuals, mem, INTSXP,
				     Information, false),
    Infos = getAttrib(Code, Information);
  Uint *pi  = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  pi[WHAT] = what; // first SNP in SNP sequence
  addr_Uint au;
  au.a = (uintptr_t) INTEGER(Code);
  for (Uint i=0; i<sizeof(uintptr_t)/BytesPerUnit; i++) pi[ADDR0 + i] = au.u[i];

  if (PRESERVE) {BUG;  R_PreserveObject(Code);}
  return(Code);
}



void codeInner2(Uint *M1, Uint *M2, Uint snps, Uint *code) {
  Uint j,
    fullBlocks = snps / CodesPerBlock;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
  for (j=0; j<fullBlocks; j++) {
    Uint *C = code + j * UnitsPerBlock,
      *mm1 = M1 + j * CodesPerBlock,
      *mm2 = M2 + j * CodesPerBlock;
     for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = u * BitsPerCode;
      for (Uint k=0; k<UnitsPerBlock; k++, mm1++, mm2++) {
	assert(*mm1 == 0 || *mm1 == 1);
	assert(*mm2 == 0 || *mm2 == 1);
	C[k] |= (*mm1 << shift) | (*mm2 << (shift + 1));
      }
    } 
  }
  if (fullBlocks * CodesPerBlock < snps) {
    Uint *end = M1 + snps,
      *C = code + fullBlocks * UnitsPerBlock;
    M1 += CodesPerBlock * fullBlocks;
    M2 += CodesPerBlock * fullBlocks;
    for (Uint u=0; u<CodesPerUnit; u++) {
     
      Uint shift = u * 2;
      for (Uint k=0; k<UnitsPerBlock && M1 < end; k++, M1++, M2++) {
	assert(*M1 == 0 || *M1 == 1);
	assert(*M2 == 0 || *M2 == 1);
	C[k] |= (*M1 << shift) | (*M2 << (shift + 1));
      }
      if (M1 == end) break;
    }
  }
}


void codeInner(Uint *M, Uint snps, Uint *code) {  
  Uint j,
    bits = snps * BitsPerCode,
    fullBlocks = bits / BitsPerBlock;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
  for (j=0; j<fullBlocks; j++) {
    Uint *C = code + j * UnitsPerBlock,
      *mm = M + j * BitsPerBlock;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = u * BitsPerCode;
      for (Uint k=0; k<UnitsPerBlock; k++, mm++) {
	assert(*mm == 0 || *mm == 1);
	C[k] |= *mm << shift;
	mm++;
	assert(*mm == 0 || *mm == 1);
	C[k] |= *mm << (shift + 1);
      }
    }
  }
  if (fullBlocks * CodesPerBlock < snps) {
    Uint *end = M + bits,
      *C = code + fullBlocks * UnitsPerBlock;
    M += BitsPerBlock * fullBlocks;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = u * 2;
      for (Uint k=0; k<UnitsPerBlock && M < end; k++, M++) {
	assert(*M == 0 || *M == 1);
	C[k] |= *M << shift;
	M++; if (M == end) break;
	assert(*M == 0 || *M == 1);
	C[k] |= *M << (shift + 1);
      }
      if (M == end) break;
    }
  }
}

void assert_shuffle() {
  int method = GLOBAL.relationship.method;
  if (method != Shuffle) // will work if aligned: ( && method != TwoBit)
    ERR("'Haplo' can be coded only for method 'Shuffle'");
  if (Information == R_NilValue) Init();
}

SEXP codeHaplo2(SEXP M1, SEXP M2) {// 2 x Haplo - Matrix, as vector gespeichert
  assert_shuffle();
  if (length(M1) == 0 || length(M1) != length(M2))
    ERR("'M' has length 0 or lengths are unequal.");
  if (!isVector(M1) || !isVector(M2))
    ERR("'M' or 'N' not a vector or of wrong size.");
  Uint snps = (Uint) length(M1);  
  SEXP Code = create_codevector(HAPLO, snps, 1);
  ToInt(M1);  
  ToInt(M2);  
  codeInner2(M1int, M2int, snps, Align(Code, ALIGN_HAPLO));
  if (PRESERVE) R_PreserveObject(Code);
  FREEint(M1);
  FREEint(M2);
  return Code;
}


SEXP codeHaplo(SEXP M) {// 2 x Haplo - Matrix, as vector gespeichert
  assert_shuffle();
  if (length(M) == 0) ERR("'M' has length 0.");
  Uint 
    nrow = 2,
    snps = (Uint) ncols(M);
  if (!isMatrix(M) || (Uint) nrows(M) != nrow)
    ERR("'M' not a matrix or of wrong size.");

  SEXP Code = create_codevector(HAPLO, snps, 1);
  ToInt(M);
  codeInner(Mint, snps, Align(Code, ALIGN_HAPLO));
  if (PRESERVE) R_PreserveObject(Code);
  FREEint(M);
  return Code;
}


SEXP decodeHaplo(SEXP CM) {
  if (Information == R_NilValue) Init();
  Uint 
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(CM, Information), INFO_INFO)),
    snps = info[SNPS], 
    bits = snps * BitsPerCode,
    fullBlocks = bits / BitsPerBlock,
    *cm = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != HAPLO) ERR("not a haplo coding");
  
  SEXP Ans;
  PROTECT(Ans = allocMatrix(INTSXP, BitsPerCode, snps));
  Uint *M = (Uint*) INTEGER(Ans);
  for (Uint i=0, endfor= BitsPerCode * snps; i<endfor; i++) M[i] = 0; // noetig?
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint j=0; j<fullBlocks; j++) {
    Uint
      *C = cm + j * UnitsPerBlock,
      *mm = M + j * BitsPerBlock;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = 1 << u * 2,
	shiftP1 = shift << 1;
      for (Uint k=0; k<UnitsPerBlock; k++, mm++) {
	*mm = (C[k] & shift) > 0;
	mm++;
	*mm = (C[k] & shiftP1) > 0;
      }
    }
  }
  if (fullBlocks * CodesPerBlock < snps) {
    Uint
      *C = cm + fullBlocks * UnitsPerBlock,
      *end = M + bits;
    M += BitsPerBlock * fullBlocks;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = 1 << u * 2,
	shiftP1 = shift << 1;
      for (Uint k=0; k<UnitsPerBlock && M < end; k++, M++) {
	*M = (C[k] & shift) > 0;
	M++; if (M == end) break;
	*M = (C[k] & shiftP1) > 0;
	}
      if (M == end) break;
    }
  }
 
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return Ans;
}


void haplo2geno(BlockType0 *X, Uint snps, Uint individuals) {
  assert_shuffle();
  Uint blocks = Blocks(snps);
  const BlockType
    EN = SET_EPI32(0xAAAAAAAA),
    NE = SET_EPI32(0x55555555); 

  for (Uint j=0; j<individuals; j++) {
    BlockType *x = X + j * blocks;
    for (Uint b=0; b<blocks; b++) {
      BlockType y = AND(x[b], EN),
	z = AND(x[b], NE);
      x[b] = OR(AND(y, SHL(z, 1L)), XOR(SHR(y, 1), z));
      // x[b] = OR(SHR(y, 1L), OR(z, AND(y, SHL(z, 1L))));
    }
  }
}
  
SEXP haplo2geno(SEXP H) {
 if (Information == R_NilValue) Init();
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(H, Information), INFO_INFO)),
    snps =  (Uint) info[SNPS];
  if (info[WHAT] != HAPLO) ERR("not a haplo coding");
  Uint
    *haplo =  Align(H, ALIGN_HAPLOGENO);
  if (info[INDIVIDUALS] != 1)
    ERR("currently only 1 individual allowed in haplo2geno");

  SEXP Code = create_codevector(GENO, snps, info[INDIVIDUALS]);
  Uint *code = algn(INTEGER(Code));
  assert((uintptr_t) code >= (uintptr_t) INTEGER(Code));
  MEMCOPY(code, haplo, Blocks(snps) * BytesPerBlock);
  haplo2geno((BlockType0 *) code, snps, 1);
  if (PRESERVE) R_PreserveObject(Code);
    return Code;  
}


void codeInnerGeno(Uint *M, Uint snps, Uint *code) {
  // Vektorisieren ?!
  Uint j,
    bits = snps * BitsPerCode,
    fullBlocks = (bits - 1)  / BitsPerBlock;

  for (j=0; j<fullBlocks; j++) {
    Uint *C = code + j * UnitsPerBlock,
      *mm = M + j * CodesPerBlock;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = u * 2;
      for (Uint k=0; k<UnitsPerBlock; k++, mm++) {
	int m = *mm;
	if (checking && (m < 0 || m > 2)) ERR("wrong values for SNP matrix");
	//if (m == 2) m = 3; // schneller? schoener?
	C[k] |= m << shift; 
      }
    }
  }
  if (fullBlocks * CodesPerBlock < snps) {
    Uint *end = M + bits / BitsPerCode,
      *C = code + UnitsPerBlock * fullBlocks;
    M += CodesPerBlock * fullBlocks;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = u * 2;
      for (Uint k=0; k<UnitsPerBlock && M < end; k++, M++) {
	int m = *M;
	if (checking && (m < 0 || m > 2)) ERR("wrong values for SNP matrix");
	//	*(p++) += (double) m;
	
	//if (m == 2) m = 3;
	C[k] |= m << shift;
      }
      if (M == end) break;
    }
  }
}


static Uint code_oldSNP = 0,
  code_oldIndiv = 0;

SEXP codeGeno(SEXP M) {// 2 x Haplo - Matrix, as vector gespeichert
  if (Information == R_NilValue) Init();
  if (length(M) == 0) ERR("'M' has length 0.");
  Uint
    what,
    snps,
    individuals;
  if (isMatrix(M)) {
    snps = nrows(M);
    individuals = ncols(M);
    if (snps < individuals && PL > 0) {
      if (code_oldSNP != snps || code_oldIndiv != individuals) {
	PRINTF("less SNPS than individuals: are you sure you have given an snp x individual matrix?\n");
      }
      code_oldSNP = snps;
      code_oldIndiv = individuals;
    }
    what = GENOMATRIX;
  } else {
    snps = length(M);
    individuals = 1;
    what = GENO;
  }

  SEXP Code = create_codevector(what, snps, individuals);
  ToInt(M);
  Uint 
    *C  = Align(Code, ALIGN_HAPLO),
    blocks = Blocks(snps),
    units = blocks * UnitsPerBlock;
  
  //double *p = REAL(VECTOR_ELT(getAttrib(Code, Information), INFO_P)) + INFO_P_P;
    
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
 for (Uint i=0; i<individuals; i++)
   codeInnerGeno(Mint + i * snps, snps, C + i * units);
  if (PRESERVE) R_PreserveObject(Code);
  FREEint(M);
 return Code;
}


SEXP decodeGenoXX(SEXP CM) { // 0.724 bei 25000 x 25000
 if (Information == R_NilValue) Init();
  Uint 
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(CM, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    bits = snps * BitsPerCode,
    fullBlocks = bits / BitsPerBlock,
    blocks = Blocks(snps),
    units = blocks * UnitsPerBlock,
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENO && info[WHAT] != GENOMATRIX) ERR("not a geno coding");
  if (info[WHAT] == GENO && info[INDIVIDUALS] != 1) BUG;
  
  SEXP Ans;
  if (info[WHAT] == GENO) PROTECT(Ans = allocVector(INTSXP, snps));
  else PROTECT(Ans = allocMatrix(INTSXP, snps, individuals));
  Uint *M = (Uint*) INTEGER(Ans);
  for (Uint i=0, endfor= snps; i<endfor; i++) M[i] = 0; // noetig?
 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint i=0; i<individuals; i++) {
    Uint *C = C0 + i * units,
      *mm = M + i * snps;
    for (Uint j=0; j<fullBlocks; j++, C += UnitsPerBlock) {
      for (Uint u=0; u<CodesPerUnit; u++) {
	Uint shift = 1 << u * 2,
	  shiftP1 = shift << 1;
	for (Uint k=0; k<UnitsPerBlock; k++, mm++) {
	  *mm = ((C[k] & shift) > 0) + ((C[k] & shiftP1) > 0);
	}
      }
    }
    
    if (fullBlocks < blocks) {
      Uint *end = M + (i+1) * snps;
      for (Uint u=0; u<CodesPerUnit; u++) {
	Uint shift = 1 << u * 2,
	  shiftP1 = shift << 1;
	for (Uint k=0; k<UnitsPerBlock && mm < end; k++, mm++) {
	  *mm = ((C[k] & shift) > 0) + ((C[k] & shiftP1) > 0);
	}
	if (mm == end) break;
      }
    }
  }
 
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);

  return Ans;
}

SEXP decodeGeno(SEXP CM) { // 0.45 bei 25000 x 25000
 if (Information == R_NilValue) Init();
  BlockType first2bits32 = SET_EPI32(0x00000003);
  Uint 
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(CM, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    bits = snps * BitsPerCode,
    fullBlocks = bits / BitsPerBlock,
    blocks = Blocks(snps),
    //    units = blocks * UnitsPerBlock,
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENO && info[WHAT] != GENOMATRIX) ERR("not a geno coding");
  if (info[WHAT] == GENO && info[INDIVIDUALS] != 1) BUG;
  
  SEXP Ans;
  if (info[WHAT] == GENO) PROTECT(Ans = allocVector(INTSXP, snps));
  else PROTECT(Ans = allocMatrix(INTSXP, snps, individuals));
  Uint *M = (Uint*) INTEGER(Ans);
  for (Uint i=0, endfor= snps; i<endfor; i++) M[i] = 0; // noetig?

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint i=0; i<individuals; i++) {
    Uint *mm = M + i * snps;
    BlockType *C = (BlockType *) C0 + i * blocks;
    for (Uint j=0; j<fullBlocks; j++, C++) {
      BlockType c = *C;
      for (Uint u=0; u<CodesPerUnit; u++, mm+=UnitsPerBlock) {	  
	BlockUnitType d;
	d.b = AND(c, first2bits32);
#define MD(i) mm[i] = d.u[i];
#if UnitsPerBlock == 4
	MD(0); MD(1); MD(2); MD(3);
#elif UnitsPerBlock == 8
	MD(0); MD(1); MD(2); MD(3); MD(4); MD(5); MD(6); MD(7);
#else	    
	for (Uint k=0; k<UnitsPerBlock; k++)  MD(k);
#endif	    
	c = SHR(c, BitsPerCode);
      }
    }
    if (fullBlocks < blocks) {
      Uint *end = M + (i+1) * snps;
       BlockType c = ((BlockType *) C0)[i * blocks + fullBlocks];
      for (Uint u=0; u<CodesPerUnit && mm < end; u++) {
	BlockUnitType d;
	d.b = AND(c, first2bits32);
	for (Uint k=0; k<UnitsPerBlock && mm < end; k++, mm++)  *mm = d.u[k];
	c = SHR(c, BitsPerCode);
      } 
    }
  }
 
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return Ans;
}


SEXP copyGeno(SEXP CM) {
  if (Information == R_NilValue) Init();
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(CM, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    what = info[WHAT];
  if (what != GENOMATRIX) ERR("not a geno matrix");

  //  printf("what = %d %d; %d %d\n", what, GENOMATRIX, snps, individuals);
  
  SEXP Code = create_codevector(what, snps, individuals);
  Uint 
    *from  = Align(CM, ALIGN_HAPLO),
    *to = Align(Code, ALIGN_HAPLO),
    blocks = Blocks(snps);
    
  MEMCOPY(to, from, blocks * BytesPerBlock * individuals);
  return Code;
}



SEXP decodeNthGeno(SEXP CM, SEXP NN) { // 0.45 bei 25000 x 25000
  if (Information == R_NilValue) Init();
  //  BlockType first2bits32 = SET_EPI32(0x00000003);
  Uint
    first2bits32 = 0x00000003,
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(CM, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    blocks = Blocks(snps),
    //    units = blocks * UnitsPerBlock,
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENOMATRIX) ERR("not a genoimatrix coding");
  
  SEXP Ans;
  if (lenN ==1) PROTECT(Ans = allocVector(INTSXP, individuals));
  else PROTECT(Ans = allocMatrix(INTSXP, lenN, individuals));
  Uint *M = (Uint*) INTEGER(Ans);

  for (Uint i=0; i<individuals; i++) {    
    Uint *mm = M + i * lenN;
    for (Uint ni=0; ni<lenN; ni++, mm++) {
      Uint j=N[ni]/ CodesPerBlock,
	jj = N[ni] % CodesPerBlock,
	u = jj / UnitsPerBlock,
	k = jj % UnitsPerBlock;     
      BlockType c = ((BlockType *) C0)[i * blocks + j];
      BlockUnitType d;
      d.b = c;
      // d.b = AND(SHR(c, BitsPerCode * u), first2bits32);  *mm = d.u[k];
      *mm = (d.u[k] >> (BitsPerCode * u)) & first2bits32;
    }
  }
 
  UNPROTECT(1);
   return Ans;
}


SEXP zeroNthGeno(SEXP CM, SEXP NN) { // 0.45 bei 25000 x 25000
  if (Information == R_NilValue) Init();
  Uint
    first2bits32 = 0x00000003,
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(CM, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    blocks = Blocks(snps),
    //    units = blocks * UnitsPerBlock,
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENOMATRIX) ERR("not a genoimatrix coding");
  
  SEXP Ans;
  if (lenN ==1) PROTECT(Ans = allocVector(INTSXP, individuals));
  else PROTECT(Ans = allocMatrix(INTSXP, lenN, individuals));
  Uint *M = (Uint*) INTEGER(Ans);

  for (Uint i=0; i<individuals; i++) {    
    Uint *mm = M + i * lenN;
    for (Uint ni=0; ni<lenN; ni++, mm++) {
      Uint j=N[ni]/ CodesPerBlock,
	jj = N[ni] % CodesPerBlock,
	u = jj / UnitsPerBlock,
	k = jj % UnitsPerBlock;     
      BlockType *C = ((BlockType *) C0) + i * blocks + j;
      BlockUnitType *d = (BlockUnitType*) C;
      d->u[k] &= ~(first2bits32 << (BitsPerCode * u));
    }
  } 
  UNPROTECT(1);
  return R_NilValue;
}



typedef union {
  BlockType b;
  HalfBlockType s[2];
} BlockTypeHalf ALIGNED;


void vectorGenoKahan(Real *V, BlockType0 *CM, Uint snps, Uint individuals,
			double *ans) {
  Uint blocks = Blocks(snps);
  BlockType    
    // NullEins = SET_EPI32(0x55555555),
    first2bits32 = SET_EPI32(0x00000003);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint i=0; i<individuals; i++) {
    BlockType
      *cm = CM + i * blocks;
    Real
      sum = ZEROREAL,
      corr = ZEROREAL,
      *v = V;
    for (Uint b=0; b<blocks; b++, cm++) {
      BlockTypeHalf d;
      BlockType c = *cm;
      //   d.b = SHR(c, 1);
      // d.b = AND(NullEins, d.b);      
      // c = XOR(c, d.b);      
      for (Uint j=0; j<CodesPer32; j++) {
	d.b = AND(c, first2bits32);
#ifndef DO_FLOAT	
	Real f = INT2REAL(d.s[0]);
#else
	Real f = INT2REAL(d.b);
#endif	
	f = MULTREAL(*v, f);

	Real fcorr = SUBREAL(f, corr),
	  roughsum = ADDREAL(sum, fcorr);
	corr = SUBREAL(SUBREAL(roughsum, sum), fcorr);	
	sum = roughsum;

#ifndef DO_FLOAT
 	v++;
	f = INT2REAL(d.s[1]);
	f = MULTREAL(*v, f);	
	fcorr = SUBREAL(f, corr);
	roughsum = ADDREAL(sum, fcorr);
	corr = SUBREAL(SUBREAL(roughsum, sum), fcorr);	
	sum = roughsum;
#endif	
	v++;
	c = SHR(c, 2L);
      }
    }
    real *s = (real*) &sum;
    ans[i] = (double) (s[0] + s[1]
#if defined AVX2 or defined DO_FLOAT
		     + s[2] + s[3]
#if defined AVX2 and defined DO_FLOAT    
		     + s[4] + s[5] + s[6] + s[7]
#endif		     
#endif		     
		       );
  }
}


void vectorGenoIntern(Real *V, BlockType0 *CM, Uint snps, Uint individuals,
			double *ans) {

  if (GLOBAL_UTILS->basic.kahanCorrection) {
    vectorGenoKahan(V, CM, snps, individuals, ans);
    return;
  }
  Uint blocks = Blocks(snps);
  BlockType    
    // NullEins = SET_EPI32(0x55555555),
    first2bits32 = SET_EPI32(0x00000003);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint i=0; i<individuals; i++) {
    BlockType
      *cm = CM + i * blocks;
    Real
      sum = ZEROREAL,
      *v = V;
    for (Uint b=0; b<blocks; b++, cm++) {
      BlockTypeHalf d;
      BlockType c = *cm;
      //      d.b = SHR(c, 1);
      //d.b = AND(NullEins, d.b);      
      // c = XOR(c, d.b);      
      for (Uint j=0; j<CodesPer32; j++) {
	d.b = AND(c, first2bits32);
#ifndef DO_FLOAT	
	Real f = INT2REAL(d.s[0]);
#else
	Real f = INT2REAL(d.b);
#endif	
	f = MULTREAL(*v, f);
	sum = ADDREAL(sum, f);
#ifndef DO_FLOAT
 	v++;
	f = INT2REAL(d.s[1]);
	f = MULTREAL(*v, f);
	sum = ADDREAL(sum, f);
#endif	
	v++;
	c = SHR(c, 2L);
      }
    }
    real *s = (real*) &sum;
    ans[i] = (double) (s[0] + s[1]
#if defined AVX2 or defined DO_FLOAT
		     + s[2] + s[3]
#if defined AVX2 and defined DO_FLOAT    
		     + s[4] + s[5] + s[6] + s[7]
#endif		     
#endif		     
		       );
  }
}

SEXP vectorGeno(SEXP V, SEXP Z) {
  if (Information == R_NilValue) Init();
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(Z, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    len = length(V);
  if (len != snps) ERR("vector 'V' not of correct length");

  real *v0 = (real *) CALLOC(RealAlign(len), sizeof(real)),
    *v = (real*) algnrl(v0);    
  switch(TYPEOF(V)) {
  case INTSXP : {
    ToInt(V);
    for (Uint i=0; i<len; i++) v[i] = (real) Vint[i];
    FREEint(V);
    break;
  }
  case REALSXP : {
    double *vp = REAL(V);
#ifdef DO_FLOAT
  for (Uint i=0; i<len; i++) v[i] = (real) vp[i];
#else
  MEMCOPY(v, vp, len * sizeof(double));
#endif
    break;
  }
  case LGLSXP : {
    int *vp = LOGICAL(V);
    for (Uint i=0; i<len; i++) v[i] = (real) vp[i];
    break;
  }
  default : ERR("'V' of unknown type.");
  }

  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, individuals));
  double *ans = REAL(Ans);
  for (Uint i=0, endfor= individuals; i<endfor; i++) ans[i] = 0.0; // noetig?
  vectorGenoIntern((Real *) v, (BlockType0 *) Align(Z, ALIGN_VECTOR),
		   snps, individuals, REAL(Ans));

  Free(v0);
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return(Ans);
}

 


/* 

void genoVectorIntern(BlockType0 *CM, double *V, Uint snps, Uint individuals,
			double *ans) {

#if defined AVX2
  BUG;
  //must be rewritten completely using

  //	const int mask1 = MOVEMASK((Float) d);
  //	d = SHL(d, 1);
  //	const int mask2 = MOVEMASK((Float) d);		
  //	Double f = ADDDOUBLE(BLEND(Null, factor, mask1),
  //	                     BLEND(Null, factor, mask2));
  //	 _mm256_movemask_ps _mm256_blend_pd
  //	 for AVX512 again different!_mm512_mask_inserti64x4_mm_insert_epi64
       
#endif
  
  if (GLOBAL_UTILS->basic.kahanCorrection) {
    //    genoVectorKahan(V, CM, snps, individuals, ans);
    // return;
  }
  Uint blocks = Blocks(snps);
  BlockType
    NullEins = SET_EPI32(0x55555555),
    first2bits32 = SET_EPI32(0x00000003),
    lastbit32 = SET_EPI32(0x80000000);
  double *Vmem = (double*) CALLOC((individuals + CodesPerBlock)/DoublesPerBlock,
				  BytesPerBlock),
    *ansmem = (double*) CALLOC((snps + CodesPerBlock) / DoublesPerBlock + 1,
			       BytesPerBlock),
    *ans_aligned = (double *) algn(ansmem);
   MEMCOPY(Vmem, V, individuals * sizeof(double));
   Double Null = ZERODOUBLE;

  for (Uint i=0; i<individuals; i++) {
    double *a = ans_aligned;
    BlockType
      *cm = CM + i * blocks;
    Double factor = LOAD1_DOUBLE(Vmem + i);
    for (Uint b=0; b<blocks; b++) {
     
      
      BlockType d,
	c = cm[b];
      for (Uint j=0; j<CodesPer32; j++) {	
	//	d.b = AND(c, first2bits32);
	BlockUnitType mask1,mask2;
	mask1.b = AND(c, lastbit32);
	c = SHL(c, 1L);	
	mask2.b = AND(c, lastbit32),
	c = SHL(c, 1L);
	
	
	const int mask1 = MOVEMASK((Float) d);
	d = SHL(d, 1);
	const int mask2 = MOVEMASK((Float) d);
	
	
	Double f = ADDDOUBLE(BLEND(Null, factor, mask1),
	                     BLEND(Null, factor, mask2));

	//	f = ADDDOUBLE(f,  _mm_blendv_pd((Double) SHR((BlockType)f, 3), (Double) SHR((BlockType)factor, 1), factor));
	
	STORE_DOUBLE(a, ADDDOUBLE(LOAD_DOUBLE(a), f));
	a += DoublesPerBlock;

       
	f = ADDDOUBLE(BLEND(Null, factor, mask1 >> 2),
	  BLEND(Null, factor, mask2 >> 2));

	//	f = ADDDOUBLE(f,  _mm_blendv_pd((Double) SHR((BlockType)f, 3), (Double) SHR((BlockType)factor, 1), factor));
	
	STORE_DOUBLE(a, ADDDOUBLE(LOAD_DOUBLE(a), f));
	a += DoublesPerBlock;
	
	
	//	c = SHR(c, 2L);	
 	c = SHL(c, 1L);	
      }
	
    }
  }

	if (false)   MEMCOPY(ans, ans_aligned, snps * sizeof(double));
   double *pans = ans_aligned;
#define dbytes (sizeof(double) * UnitsPerBlock)
   for (Uint b=0; b<blocks; b++, pans += CodesPerBlock) {
     for (Uint j=0; j<CodesPerUnit/2; j++) {
       double swap[UnitsPerBlock];
       MEMCOPY(swap, pans + j * UnitsPerBlock, dbytes);
       MEMCOPY(pans + j * UnitsPerBlock,
	       pans + CodesPerBlock - (j+1) * UnitsPerBlock, dbytes);
       MEMCOPY(pans + CodesPerBlock - (j+1) * UnitsPerBlock, swap, dbytes);
       //     printf("%d %d %ld %d\n",  j * UnitsPerBlock,CodesPerBlock  - (j+1) * UnitsPerBlock,  (uintptr_t) pans, CodesPerBlock * sizeof(double));
     }
     //    assert(b < 5);
   }

   if (false)
    for (Uint i=0; i<2*CodesPerBlock; i++) {
      //if (i==CodesPerBlock) printf("\n");	    
     //printf("%10g %10g \n", ans_aligned[i], ans[i]);
   }
   
 
   

	
  MEMCOPY(ans, ans_aligned, snps * sizeof(double));
  FREE(ansmem);
  FREE(Vmem);
}



void genoVectorIntern(BlockType0 *CM, double *V, Uint snps, Uint individuals,
			double *ans) {
  // nice code, but very slow !!
  if (GLOBAL_UTILS->basic.kahanCorrection) {
    //    genoVectorKahan(V, CM, snps, individuals, ans);
    // return;
  }
  Uint blocks = Blocks(snps);
  BlockType
    NullEins = SET_EPI32(0x55555555),
    first2bits32 = SET_EPI32(0x00000003);
  double *Vmem = (double*) CALLOC((individuals + CodesPerBlock)/DoublesPerBlock,
				  BytesPerBlock),
    *ansmem = (double*) CALLOC((snps + CodesPerBlock) / DoublesPerBlock,
			       BytesPerBlock),
    *ans_aligned = (double *) algn(ansmem);
   MEMCOPY(Vmem, V, individuals * sizeof(double));
   for (Uint i = 0; i<individuals; i++) Vmem[i] = V[snps - i];
 
  for (Uint i=0; i<individuals; i++) {
    double *a = ans_aligned;
    BlockType
      *cm = CM + i * blocks;
    Double factor = LOAD1_DOUBLE(V + i);
    for (Uint b=0; b<blocks; b++) {
      BlockTypeHalf d;
      BlockType c = cm[b];
   //   d.b = SHR(c, 1);
   //   d.b = AND(NullEins, d.b);      
   //   c = XOR(c, d.b);
      for (Uint j=0; j<CodesPer32; j++) {
	d.b = AND(c, first2bits32);
	Double f = MULTDOUBLE(factor, INT2DOUBLE(d.s[0]));
	STORE_DOUBLE(a, ADDDOUBLE(f, LOAD_DOUBLE(a)));
	a += DoublesPerBlock;

	f = MULTDOUBLE(factor, INT2DOUBLE(d.s[1]));	
	STORE_DOUBLE(a, ADDDOUBLE(f, LOAD_DOUBLE(a)));		
	a += DoublesPerBlock;
	
	c = SHR(c, BitsPerCode);	
      }
    }
  }
  MEMCOPY(ans, ans_aligned, snps * sizeof(double));
  FREE(ansmem);
  FREE(Vmem);
}

*/

void genoVectorIntern(BlockType0 *CM, double *V, Uint snps, Uint individuals,
			double *ans) {

  if (GLOBAL_UTILS->basic.kahanCorrection) {
    //    genoVectorKahan(V, CM, snps, individuals, ans);
    // return;
  }
  Uint blocks = Blocks(snps),
    blocksM1 = blocks - 1;
  double *end_ans = ans + snps;
  BlockType c,
    first2bits32 = SET_EPI32(0x00000003);
  for (Uint i=0; i<snps; ans[i++] = 0.0);

  for (Uint i=0; i<individuals; i++) {
    BlockType *cm = CM + i * blocks;
    double
      *a = ans,
      f0 = V[i],
      //  factor[4] = {0.0, f0, f0, 2.0 * f0}; // so also ok for haplo types!
      //                  mit der 00, 01, 11 Codierung
      factor[4] = {0.0, f0, 2.0 * f0, 0}; // nur fuer geno mit 00, 01, 10 Code
    for (Uint b=0; b<=blocksM1; b++) {
      // BlockType
	c = cm[b];
      Uint endfor = b < blocksM1 ? CodesPer32
			: (snps - blocksM1 * CodesPerBlock) / UnitsPerBlock;
      for (Uint j=0; j<endfor; j++) {
	BlockUnitType d;
	d.b = AND(c, first2bits32);
	a[0] += factor[d.u[0]]; 
	a[1] += factor[d.u[1]];
	a[2] += factor[d.u[2]];
	a[3] += factor[d.u[3]];
#if defined AVX2
	a[4] += factor[d.u[4]]; 
	a[5] += factor[d.u[5]];
	a[6] += factor[d.u[6]];
	a[7] += factor[d.u[7]];
	a += 8;
#else
	a += 4;
#endif 
	c = SHR(c, BitsPerCode);
      }
    }

    BlockUnitType d;
    d.b = AND(c, first2bits32);
    assert(end_ans - a < UnitsPerBlock && end_ans - a >= 0);
    for(Uint k=0; a < end_ans; a++) *a = factor[d.u[k++]];
  }  
}


SEXP genoVector(SEXP Z, SEXP V) {
  if (Information == R_NilValue) Init();
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(Z, Information), INFO_INFO)),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    len = length(V);
  if (len != individuals) ERR("vector 'V' not of correct length");
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, snps));
  genoVectorIntern((BlockType0 *) Align(Z, ALIGN_VECTOR), REAL(V), 
 		   snps, individuals, REAL(Ans));
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return(Ans);
}

 

BlockType INLINE VECTORADD_ALTERNATIVE(BlockType0 xx, BlockType0 yy) {
  // 1 Operation more because of highest bit
  // Idee: die werte direkt eintragen, statt ueber die Tabelle nachschauen:
  BlockType z, sum, y, x;
  // 2 * 1 oder 2 * 1
  x = OR(xx, yy);
  z = SHL(x, 1L);
  x = AND(x, z);
  x = AND(x, MASK2); // ergb
    
  z = AND(xx, yy);
  y = AND(z, MASK1); // 1 * 1
  x = OR(x, y);     // "addition" zu "2 * 1"
  
  z = AND(z, MASK2); // 2 * 2;

  // globaler shift left leider nicht moeglich da oberstes Bit runterfaellt
  // somit werden die beiden (*) Zeilen zusaetzlich benoetigt
  // ohne die beiden unteren Zeilen hat code nur 14 (statt VECTADD 14)
  // Operationen mit einem Zeitgewinn von somit 7% und einem damit notwendigern
  // Speichmehrbedarf von 1.5% bei SSE und 0.8% bei AVX
  y = SHL(z, 1L);
  y = OR(z, y);
  sum = AND(y, LH);

  y = SHR(x, 4L);
  z = SHR(z, 3L); // (*)
  y = OR(y, z);   // (*)
  y = AND(y, LH);
  
  return ADD8(sum, y);
}



BlockType INLINE VECTORADD(BlockType0 xx, BlockType0 yy) {
  BlockType z, sum, y;  
  z = AND(xx, yy);
  y = AND(z, LH); // 0.25 von 4
  sum = SHUFFLE8(and2scalar, y); // 0.5 von 4

  y = SHR(z, 4L);
  y = AND(y, LH); // 
  y = SHUFFLE8(and2scalar, y);// 
  sum = ADD8(sum, y); // 
  //   showblock(sum, CodesPerBlock);

#ifndef OLD_VECTORADD
  z = XOR(xx, yy);
  y = AND(z, LH);
  y = SHUFFLE8(xor2scalar, y);
  sum = ADD8(sum, y);

  y = SHR(z, 4L);
  y = AND(y, LH);
  y = SHUFFLE8(xor2scalar, y);
#else  // FOLGENDES IST NICHT SCHNELLER !!!
  z = XOR(xx, yy);
  y = SHR(z, 1L);
  z = AND(y, z);
  y = SHR(z, 3L);
  z = OR(z, y);
  y = AND(z, LH);
  y = SHUFFLE8(xor2scalar, y);
#endif
  return ADD8(sum, y);
}

Uint INLINE ADDUPVECTOR(BlockType0 x) {
#if defined AVX2
  BlockUnitType vsum;
  vsum.b = SAD8(x, ZERO);  
  return vsum.u[0]+ vsum.u[2] + vsum.u[4] + vsum.u[6];
#else 
  BlockType vsum = SAD8(x, ZERO); //
  return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
#endif
}
 
Uint scalar012(BlockType0 *x, BlockType0 *y, Uint blocks) {
#define x15  (255L / 16L) // 4 : max=2^2
  //                         4 : 2-bit-coding in 1 Byte                    
  Uint
    bundles = blocks / x15,
    ans = 0;

  for (Uint i=0; i<bundles; i++) { 
    BlockType
      *xx = x + i * x15,
      *yy = y + i * x15,
      sum = ZERO;
#if x15 == 15
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;

      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
#else
    for (Uint j=0; j<x15; j++) {
      //      sum = ADD8(sum, VECTORADD(LOADU(xx + j), LOADU(yy + j))) ;
      sum = ADD8(sum, VECTORADD(LOAD(xx++), LOAD(yy ++))) ;
      // sum = ADD8(sum, VECTORADD(xx[j], yy[j])); // T: gleichauf mit vorigem
      //      sum = ADD8(sum, VECTORADD(_mm_loadu_si128(xx + j), // langsamer
      //			_mm_loadu_si128(yy + j))); 
    }
 #endif
   ans += ADDUPVECTOR(sum); 
  }
  BlockType sum = ZERO,
    *endx = x + blocks,
    *xx = x + bundles * x15,
    *yy = y + bundles * x15; 
  for (; xx < endx; xx++, yy++) sum = ADD8(sum, VECTORADD(*xx, *yy));
  ans += ADDUPVECTOR(sum);
  return ans;
}

void allele_freq2(BlockType0* CGM, Uint snps, int individuals,
		double *sumP, double *sumP2) {
  // NOTE! twice allel frequency (indicated by the 2 at the end) 
  double n =  (double) individuals;
  if (snps > 536870912) ERR("too many SNPs");
  Uint blocks = Blocks(snps);
  
#define x31  (255L / 8L) // 2 : max{0,1,2}
  //                         4 : 2-bit-coding in 1 Byte                    

  Uint
    byteblocks = (snps - 1L) / BytesPerBlock + CodesPerByte,
    alignedByteBlocks = byteblocks + 1,
    units = byteblocks * BytesPerBlock,
    alignedUnits = units + UnitsPerBlock,
    // here: P= 2(p_i), not 2 p_i - 1 == expected values == means
    groupsM1 = (individuals - 1L) / x31;
  int
    *Sums0 = (int*) MALLOC(alignedUnits * BytesPerUnit),      
    *sum0 = (int*) MALLOC(alignedByteBlocks * BytesPerBlock);
  BlockType *sum = (BlockType0*) algn(sum0),
    *Sums  = (BlockType0*) algn(Sums0),
    // NullEins = SET_EPI32(0x55555555),
    first8bits = SET_EPI32(0x000000FF),
    first2bits8 = SET_EPI32(0x03030303);
  if (Sums == NULL || sum == NULL) ERR("allocation error");
  { BlockType zero = ZERO;
    Uint endfor = units / UnitsPerBlock;
    for (Uint k=0; k<endfor; Sums[k++] = zero); }
  // parallel funktioniert nicht gut!
  for (Uint i=0; i<=groupsM1; i++) {
    BlockType *cgm = CGM + i * blocks * x31;
    Uint endx = i < groupsM1 ? x31 : individuals - x31 * groupsM1;
    { BlockType zero = ZERO; for (Uint k=0; k < byteblocks; sum[k++] = zero);}
    for (Uint x=0; x<endx; x++) {
      BlockType *s = sum;
      for (Uint k=0; k<blocks; k++, cgm++) { // see permutation.tex (*)
	BlockType d,
	  c = *cgm;
       	//	d = SHR(c, 1);
	//	d = AND(NullEins, d);
	//	c = XOR(c, d);      
	for (Uint j=0; j<CodesPerByte; j++, s++) { // for-schleife aufloesen?
	  d = AND(c, first2bits8);
	  *s = ADD8(*s, d);
	  c = SHR(c, 2);
	}
      }
    }
    
    //    R < Simianer_sims_23_01.R --vanilla
	 
    BlockType *pS = (BlockType0*) Sums;      
    for (Uint k=0; k<byteblocks; k++) { // see permutation.tex (*)
      BlockType s0 = sum[k];
      for (Uint j=0; j<BytesPerUnit; j++, pS++) {
	*pS = ADD32(*pS, AND(s0, first8bits)); 
	s0 = SHR(s0, BitsPerByte);
      }
    }
  } // for groups
  FREE(sum0);
  
  Uint *S = (Uint*) Sums,
    alignedreal = byteblocks * BytesPerBlock * sizeof(real) + BytesPerBlock;
  Pd0 = (real*) MALLOC(alignedreal);
  real
    *Pd = (real *) algnrl(Pd0),
    *pPd = Pd;
  double
    sum_Psq = 0.0,
    sum_P = 0.0;
  if (GLOBAL_UTILS->basic.kahanCorrection) {
    double corr = 0.0,
      corrSq = 0.0;
    for (Uint b=0,i=0; b<blocks; b++, pPd += CodesPerBlock) {
      for (Uint c=0; c<CodesPerBlock; c++, i++) {
	double Pi = (double) S[i] / n;	
	pPd[permut_finv_g_h[c]] = (real) Pi;

	double fcorr = Pi - corr,
	  roughsum = sum_P + fcorr;
	corr = (roughsum - sum_P) - fcorr;	
	sum_P = roughsum;
	
	fcorr = Pi * Pi - corrSq;
	roughsum = sum_Psq + fcorr;
	corrSq = (roughsum - sum_Psq) - fcorr;	
	sum_Psq = roughsum;
      }
    }
  } else {
    for (Uint b=0,i=0; b<blocks; b++, pPd += CodesPerBlock) {
      for (Uint c=0; c<CodesPerBlock; c++, i++) {
	double Pi = (double) S[i] / n;     
	pPd[permut_finv_g_h[c]] = (real) Pi;
	sum_P += Pi;
	sum_Psq += Pi * Pi;
      }
    }
  }
    
  FREE(Sums0);
  if (sumP != NULL) *sumP = sum_P;
  if (sumP2 != NULL) *sumP2 = sum_Psq;    
}

      

 
void oldIrelationshipMatrix(BlockType0* CGM, Uint snps,
			Uint individuals,
			bool centred, bool normalized, double *ans) {
  Uint blocks = Blocks(snps);
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 16)  
#endif
  for (Uint i=0; i<individuals; i++) {
    for (Uint j=i; j<individuals; j++) {
      Uint Mij = scalar012(CGM + i * blocks, CGM + j * blocks, blocks);
      ans[i + j * individuals] = ans[i * individuals + j] = (double) Mij;
    }
  }
  
  if (centred || normalized) {
    double sum_P, sum_Psq;
    allele_freq2(CGM, snps, individuals, &sum_P,
		 &sum_Psq); // Pd0 set as side effect
    real *Pd = (real*) algnrl(Pd0);
    if (centred) {
      double halfppt = 0.5 * sum_Psq,
	*PGt =  (double*) MALLOC(individuals * sizeof(double));
      vectorGenoIntern((Real *) Pd, CGM, snps, individuals, PGt);
      for (Uint i=0; i<individuals; i++) {
	PGt[i] -= halfppt; // reduce GtG - PtG - GtP + PtP to GtG - PtG - GtP
      }
      for (Uint i=0; i<individuals; i++) {
	Uint k = i * (individuals + 1);// starting with diagonal element in ans
	double PGti = PGt[i];
	for (Uint j=i; j<individuals; j++, k++) {
	  ans[j * individuals + i] = ans[k] = ans[k] - PGt[j] - PGti;
	}
      } 
      FREE(PGt);
    }
    FREE(Pd0);
    if (normalized) {
      double sum_pi_qi = sum_P - 0.5 * sum_Psq;
      // sum_pi_qi == 2 sum pi(1-pi) with pi = Pi / 2
      Uint Isq = individuals * individuals;
      for (Uint i=0; i<Isq; i++) ans[i] /= sum_pi_qi;
    }
  }
}


uint64 sumGeno(BlockType0 *S, Uint blocks) {  
#define x31  (255L / 8L) // 2 : max{0,1,2}
  //                         4 : 2-bit-coding in 1 Byte
#define BytesPer64  (64 / BitsPerByte)
 //#define intgroup_size (2147483648L / (x31 * 8L * BytesPerUnit)) // 2^31 / x31
#define intgroup_size 3 // 2^31 / x31
#define intgroup_sizeM1 (intgroup_size - 1L)
  Uint
    intgroupsM1 = (blocks - 1L) / (x31 * intgroup_size),
    sizex31 = intgroup_size * x31;
  //    groupsM1 = (blocks - 1L) / x31;
  BlockType
    // NullEins = SET_EPI32(0x55555555),
    first2bits8 = SET_EPI32(0x03030303),
    first8bits = SET_EPI32(0x000000FF);
  BlockUnitType sum32;
   uint64 ans = 0;

  for (Uint X=0; X<=intgroupsM1; X++) { // this loop needed only for huge
    // data sizes of more than 50000 individuals and 10^5 SNPS each
    sum32.b = ZERO;
    BlockType *S_int = S + X * sizex31 ;
    Uint groupsM1 = X < intgroupsM1 ? intgroup_sizeM1
			: (blocks - 1L) / x31 - intgroupsM1 * intgroup_size;
    for (Uint i=0; i<=groupsM1; i++) {
      BlockType sum = ZERO,
	*cgm = S_int + i * x31;
      Uint endx = i < groupsM1 || X < intgroupsM1 ? x31 :
	blocks - intgroupsM1 * sizex31 - x31 * groupsM1;
      for (Uint x=0; x<endx; x++, cgm++) {
	BlockType d,
	  c = *cgm;
	//	d = SHR(c, 1);
	//	d = AND(NullEins, d);
	//	c = XOR(c, d);      
	for (Uint j=0; j<CodesPerByte; j++) { // for-schleife aufloesen ??
	  d = AND(c, first2bits8);
	  sum = ADD8(sum, d);
	  c = SHR(c, 2);
	}
      }
      for (Uint j=0; j<BytesPerUnit; j++) {
	sum32.b = ADD32(sum32.b, AND(sum, first8bits)); 
	sum = SHR(sum, BitsPerByte);
      }
    } // for groups

    for (Uint i=0; i<UnitsPerBlock; i++) {
      ans += (uint64) sum32.u[i];
    }
  }

  return ans;
}

void IIrelationshipMatrix(BlockType0* CGM, Uint snps,
			  Uint individuals, double * ans) {
  Uint blocks = Blocks(snps);
  double nSq = (double) individuals * (double) individuals;
       
  if ((double) snps * nSq * 4.0 > 4.5035996e+15) // 2^{manissen-bits = 52}
    ERR("matrix too large to calculate the relationship matrix");
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 16)  
#endif
  for (Uint i=0; i<individuals; i++) {
    for (Uint j=i; j<individuals; j++) {
      Uint Mij = scalar012(CGM + i * blocks, CGM + j * blocks, blocks);
      ans[i + j * individuals] = ans[i * individuals + j] = (double) Mij;
    }
  }
}
    
void CentredRelationshipMatrix(BlockType0* CGM, Uint snps,
			       Uint individuals, double *ans) {   
  Uint blocks = Blocks(snps);
  double sum = sumGeno(CGM, blocks * individuals);
  IIrelationshipMatrix(CGM, snps, individuals, ans);
  DoCentering(ans, individuals, true, true, sum);
}
    

void IrelationshipMatrix(BlockType0* CGM, Uint snps,
			 Uint individuals, SEXP Ans) {   
  Uint blocks = Blocks(snps);
  double sum = GLOBAL.relationship.normalized == True
    ? sumGeno(CGM, blocks * individuals)
    : RF_NAN;
  IIrelationshipMatrix(CGM, snps, individuals, REAL(Ans));
  DoCentering(Ans, individuals, sum);
}
    


SEXP relationshipMatrix(SEXP GM) {
  if (Information == R_NilValue) Init();
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(GM, Information), INFO_INFO)),
   individuals = info[INDIVIDUALS],
    *cgm = Align(GM, ALIGN_RELATION);
  if (info[WHAT] != GENOMATRIX) ERR("not a coded Z matrix");
 
  SEXP Ans; 
  PROTECT(Ans = allocMatrix(REALSXP, individuals, individuals));
  IrelationshipMatrix((BlockType0 *) cgm, info[SNPS], individuals, Ans);
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return Ans;
}
  



SEXP allele_freq(SEXP GM) {
  if (Information == R_NilValue) Init();
  Uint 
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(GM, Information), INFO_INFO)),
   individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    *cgm = Align(GM, ALIGN_ALLELE);
  if (info[WHAT] != GENOMATRIX) ERR("not a coded Z matrix");
  SEXP Ans; 
  PROTECT(Ans = allocVector(REALSXP, snps));
  allele_freq2((BlockType0 *) cgm, snps, individuals, NULL,
	       NULL); // Pd0 set as side effect
  real *Pd = (real*) algnrl(Pd0);
  double *ans = REAL(Ans);
  for (Uint i=0; i<snps; i++) ans[i] = 0.5 * (double) Pd[i];
  FREE(Pd0);
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
  return Ans;
}


SEXP createSNPmatrix(SEXP SNPs, SEXP Individuals) {  
  assert_shuffle();  
  Uint 
    n = ToInt0(Individuals),
    snps = ToInt0(SNPs);
  SEXP Code = create_codevector(GENOMATRIX, snps, n);
  if (PRESERVE) R_PreserveObject(Code);
  return Code;
}


SEXP fillSNPmatrix(SEXP Z, SEXP Idx, SEXP V) {
  assert_shuffle();
  ToInt(Idx);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(Z, Information), INFO_INFO)),
    *z = Align(Z, ALIGN_RELATION),
    *v = Align(V, ALIGN_HAPLO),
     len = length(Idx),
    blocks = Blocks(info[SNPS]),
    units = blocks * UnitsPerBlock,
    bytes = BytesPerBlock * blocks;
  if (info[WHAT] != GENOMATRIX) ERR("not a geno matrix that is filled");
  for (Uint i=0; i<len; i++, v += units) {
    Uint idx = Idxint[i] - 1;
    if (idx > info[INDIVIDUALS]) ERR("Idx out of bound");
    MEMCOPY(z + idx * units, v, bytes);
  }
  FREEint(Idx);
  return R_NilValue;
}




SEXP unlock(SEXP C) {
  BUG;
  const char *c[3] = WHAT_NAMES;
  if (Information == R_NilValue) Init();
  Uint 
    *info = (Uint *) INTEGER(VECTOR_ELT(getAttrib(C, Information), INFO_INFO));
  if ((SEXP) info == R_NilValue) ERR("not a coded object");
  if (PRESERVE) {
#ifdef SCHLATHERS_MACHINE  
   PRINTF("unlocking a %.50s at %ld\n", c[info[WHAT]], (uintptr_t) INTEGER(C));
#endif  
  } else warning("'unlock' is called although 'dolocking()' has not been called");
  R_ReleaseObject(C);
  return R_NilValue;
}


// population
#define INFO 0
#define BREEDING 1

// info
#define INFO_DUMMY 0
#define NSNPS 2
#define POSITION 3
#define SNP_POSITION 5
#define LENGTH_CHROM 6
#define LENGTH_TOTAL 7
#define SNPS_EQUIDIST 17
#define GENERATION_REFTABLE 18

// individual
#define RECOMBI1 0
#define RECOMBI2 1
#define MUTAT1 2
#define MUTAT2 3
#define ORIGIN1 4
#define ORIGIN2 5
//#define FATHER 6
//#define MOTHER 7
#define HAPLOTYP1 8 // in case of founder, otherwise empty
#define HAPLOTYP2 9
#define CODEDHAPLO HAPLOTYP1

#define MAX_CHROM 100

//SEXP getIndivid(SEXP Population, int generation, int sex, int nr) {}


Uint INLINE GetGenoCode(BlockType0 *code, Uint pos) {
  int block = pos / CodesPerBlock,
    unit = pos % UnitsPerBlock,
    elmt = (pos % CodesPerBlock) / UnitsPerBlock;
  BlockUnitType *c = (BlockUnitType*) code + block;
  return (c[0].u[unit] >> (2 * elmt)) & 0x00000003;
}


void INLINE PutGenoCode(Uint value, Uint pos, BlockType0 *code) {
  int 
    unit = pos % UnitsPerBlock,
    elmt = (pos % CodesPerBlock) / UnitsPerBlock;
  BlockUnitType *c = (BlockUnitType*) code + pos / CodesPerBlock;
  Uint
    shift = 2 * elmt;
  c[0].u[unit] = (c[0].u[unit] & (0xFFFFFFFF xor (0x00000003 << shift)))//clear
    | (value << shift);
}


Uint INLINE GetHaploCode(BlockType0 *code, Uint pos, Uint haplo) {
  int unit = pos % UnitsPerBlock,
    elmt = (pos % CodesPerBlock) / UnitsPerBlock;
  BlockUnitType *c = (BlockUnitType*) code + pos / CodesPerBlock;
  assert(haplo <= 1);
  assert(elmt < CodesPerUnit);
  return (c[0].u[unit] >> (2 * elmt + haplo)) & 0x00000001;
}


void INLINE PutHaploCode(Uint value, Uint pos, Uint haplo, BlockType0 *code) {
  int 
    unit = pos % UnitsPerBlock,
    elmt = (pos % CodesPerBlock) / UnitsPerBlock;
  BlockUnitType *c = (BlockUnitType*) code + pos / CodesPerBlock;
  Uint
    shift = 2 * elmt + haplo;
  if (value) c[0].u[unit] |= (0x00000001 << shift);
  else 
    c[0].u[unit] &= 0xFFFFFFFF xor (0x00000001 << shift);
}


void INLINE PutHaploOne(Uint pos, Uint haplo, BlockType0 *code) {
  int 
    unit = pos % UnitsPerBlock,
    elmt = (pos % CodesPerBlock) / UnitsPerBlock;
  BlockUnitType *c = (BlockUnitType*) code + pos / CodesPerBlock;
  Uint
    shift = 2 * elmt + haplo;
  c[0].u[unit] |= (0x00000001 << shift);
}


  
Uint getNrSNPposition(double cur_recomb, double *position, Uint Min, Uint Max) {
  Uint min = Min,
    max = Max;
  //  print("%d %d\n", min, max);
  while (min < max) {
    Uint inbetween = (min + max + 1) / 2;
    //printf("min=%d,max=%d %10g cr=%10g, pm=%10g\n", min, max,position[inbetween],
    //	   cur_recomb, position[max] );
    if (position[inbetween] <= cur_recomb) min = inbetween;
    if (position[inbetween] > cur_recomb) max = inbetween - 1;
  }
  if (min > max) BUG;
  //  printf("%10g %10g min=%d max=%d (start:%d %d) %10g %10g %10e\n", position[max], cur_recomb, min, max, Min, Max, position[Max-1], position[Max], cur_recomb -  position[Max]);
  if (position[max] > cur_recomb) {
    if (max == Min) return max - 1; // links des ersten Markers auf dem
    //                                 Chromosom; somit wird letzter marker 
    //                                 vom Vor-Chromoosom zugeordnet
    BUG;
  }
  if (max < Max && position[max+1] <= cur_recomb) BUG;
  return max;
}


#define ASSERT(WHERE,WHAT,NAME)						\
  assert(STRCMP(CHAR(STRING_ELT(getAttrib(WHERE, R_NamesSymbol), WHAT)), \
		NAME) == 0);
int IcomputeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		 int From_SNP, int To_SNP, SEXP Select, BlockType0 *ans){
  //  zaehler ++; printf("%d\n", zaehler);
  Uint
    individuals = length(Generation);
  if (individuals != (Uint) length(Sex) || individuals != (Uint) length(Nr))
    ERR("'gen', 'gender' and 'nr' do not have same length");
  if (From_SNP <= 0) ERR("'from_p' not valid.");

  Uint len_sel = length(Select);
  bool do_select = len_sel > 0;
  CondToInt(do_select, Select);
  Uint cur_select = NA_INTEGER;
    
  ASSERT(population, BREEDING, "breeding");
  SEXP breeding = VECTOR_ELT(population, BREEDING);

  ASSERT(population, INFO, "info")
  SEXP info = VECTOR_ELT(population, INFO);
  Uint len_info = length(info);

  static Uint gen_ref[MAX_GENE_INPUT] = { 0 };
  if (GENERATION_REFTABLE >= len_info ||
      length(VECTOR_ELT(info, GENERATION_REFTABLE)) == 0) {
    if (gen_ref[MAX_GENE_INPUT - 1] == 0) { // should never happen except uninit
      for (Uint i=0; i<MAX_GENE_INPUT; i++) gen_ref[i] = i;
      PRINTF("Generation reference table now initialized.\n"); 
    }
  } else {
    ASSERT(info, GENERATION_REFTABLE, "origin.gen");
    int *p = INTEGER(VECTOR_ELT(info, GENERATION_REFTABLE)); 
    for (Uint i=0; i<MAX_GENE_INPUT; i++) gen_ref[i] = (Uint) (p[i] - 1);
  }

  ASSERT(info, LENGTH_TOTAL, "length.total");
  Uint n_chrom = length(VECTOR_ELT(info, LENGTH_TOTAL)) - 1;
  if (n_chrom >= MAX_CHROM) ERR1("maximum chromosoms is %d", MAX_CHROM);
  double *length_total = REAL(VECTOR_ELT(info, LENGTH_TOTAL));
  
  ASSERT(info, LENGTH_CHROM, "length");
  double *LengthChrom = REAL(VECTOR_ELT(info, LENGTH_CHROM));
  
  ASSERT(info, NSNPS, "snp");
  Uint n_snps;
  SEXP Dummy = VECTOR_ELT(info, INFO_DUMMY);

  //   if (zaehler >= 11424) printf("AB\n");
 #ifdef DOCH_AUF_INFO_ABSPEICHERN  
   if (TYPEOF(Dummy) != INTSXP || length(Dummy) == 0){
    Uint mem = n_chromP1 + n_chromP1;
    SEXP dummy;
    PROTECT(dummy = allocVector(INTSXP, mem));    
    Uint *snps_perchrom = (Uint*) INTEGER(VECTOR_ELT(info, NSNPS)),
      *length_prec_chroms = (Uint*) INTEGER(dummy) + 0,
      *totalsnps_prec_chroms = (Uint*) INTEGER(dummy) + n_chromP1;
    length_prec_chroms[0] = 0;
    totalsnps_prec_chroms[0] = 0;
    for (Uint c=1; c<=n_chrom; c++) {
      length_prec_chroms[c] =  length_prec_chroms[c-1] + LengthChrom[c-1];
      totalsnps_prec_chroms[c] = totalsnps_prec_chroms[c-1] +
	snps_perchrom[c-1];
    }
    SET_VECTOR_ELT(info, INFO_DUMMY, dummy);
    UNPROTECT(1);
  }
  Uint
    *length_prec_chroms = (Uint*) INTEGER(VECTOR_ELT(info, INFO_DUMMY)) + 0,
    *totalsnps_prec_chroms = length_prec_chroms + n_chromP1;
#else
  if (TYPEOF(Dummy) != LGLSXP) {
    SEXP dummy;
    PROTECT(dummy = allocVector(LGLSXP, 1));
    LOGICAL(dummy)[0] = true;
    SET_VECTOR_ELT(info, INFO_DUMMY, dummy);
    UNPROTECT(1);  
   }
   
  double length_prec_chroms[MAX_CHROM + 1];
  SEXP snps_perchrom = VECTOR_ELT(info, NSNPS);
  ToInt(snps_perchrom);
  Uint
    totalsnps_prec_chroms[MAX_CHROM + 1];
  length_prec_chroms[0] = 0;
  totalsnps_prec_chroms[0] = 0;
  for (Uint c=1; c<=n_chrom; c++) {
    length_prec_chroms[c] =  length_prec_chroms[c-1] + LengthChrom[c-1];
    totalsnps_prec_chroms[c] = totalsnps_prec_chroms[c-1] +
      snps_perchromint[c-1];
  }
#endif

  // if (zaehler >= 11424) printf("ABC\n");
   n_snps = totalsnps_prec_chroms[n_chrom];
   assert(n_snps > 0);
  Uint fromSNP = From_SNP - 1,
    toSNP = To_SNP==NA_INTEGER || To_SNP >= (int) n_snps ? (Uint) MAXINT
          : (Uint) To_SNP -1;
  if (toSNP < fromSNP) ERR("'to_p' smaller than 'from_p'");
  if (fromSNP % CodesPerBlock != 0)
    ERR1("(from_p-1) must be devisible by %d", CodesPerBlock);
 if (toSNP < MAXINT && (toSNP + 1)  % CodesPerBlock != 0)
    ERR1("'to_p' must be devisible by %d or 'NA'", CodesPerBlock);


  Uint cur_nsnps = 0;
  if (do_select) {
    cur_select = Selectint[0] - 1;
    Uint toSNP_P1 = toSNP + 1;
    for (Uint k=0; k<len_sel; k++)
      cur_nsnps += Selectint[k] > fromSNP && Selectint[k] <= toSNP_P1;
  } else cur_nsnps = 1L + MIN(n_snps - 1L, toSNP) - fromSNP;
  if (ans == NULL) return cur_nsnps;
  Uint blocks = Blocks(cur_nsnps);

  ASSERT(info, SNPS_EQUIDIST, "snps.equidistant");
  bool equidist = (bool) LOGICAL(VECTOR_ELT(info, SNPS_EQUIDIST))[0];
  
  ASSERT(info, POSITION, "position");
  SEXP position = VECTOR_ELT(info, POSITION);

  ASSERT(info, SNP_POSITION, "snp.position");
  double *snp_positions = REAL(VECTOR_ELT(info, SNP_POSITION)),
    last_snp_position = snp_positions[n_snps - 1],
    first_snp_position = snp_positions[0];

  ToInt(Generation);
  ToInt(Sex);
  ToInt(Nr);

  //  if (zaehler >= 11424) printf("ABD\n");
  for (Uint i=0; i<individuals; i++) {
    //    print("i=%d %d\n", i, individuals);
    //    printf("%d %d\n", i, individuals);
    if (length(breeding) < (int) Generationint[i]) ERR("value of 'generation' too large.");
    SEXP
      g = VECTOR_ELT(breeding, Generationint[i]-1);
    if (length(g) < (int) Sexint[i]) ERR("'sex' out of range.");
    SEXP 
      s = VECTOR_ELT(g, (int) Sexint[i]-1);
    if (length(s) < (int) Nrint[i]) ERR("value of 'nr' too large..");  
    SEXP
      Indiv = VECTOR_ELT(s, Nrint[i]-1);
    BlockType *haplo = ans + i * blocks,
      NullEins[2] = { SET_EPI32(0x55555555),
		      SET_EPI32(0xAAAAAAAA)};
    for (Uint hapnr=0; hapnr<=1; hapnr++) {
      //printf("h=%d \n",  hapnr);
      //      if (zaehler >= 11424) printf("haplo=%u\n", hapnr);      
      if (TYPEOF(VECTOR_ELT(Indiv, RECOMBI1 + hapnr)) != REALSXP)
	ERR("real valued vector for recombination information expected");
      double *recombi = REAL(VECTOR_ELT(Indiv, RECOMBI1 + hapnr));
      Uint
	n_recomb = length(VECTOR_ELT(Indiv, RECOMBI1 + hapnr)) - 1;
      SEXP origin = VECTOR_ELT(Indiv, ORIGIN1 + hapnr);
      ToInt(origin);
      Uint
	n_mutat = length(VECTOR_ELT(Indiv, MUTAT1 + hapnr)),
	// separation of the following 3 allows for selection but also
	// for mutation error such as insertion and deletion
	p = 0, // position what is written if there is no selection
	q = 0, // position what is written actually
	psel = 0, pmutat = 0,
	till = 0,
	from = 0; // current recombination starting point
      bool mutations = n_mutat > 0;
      
      SEXP mutat = mutations ? VECTOR_ELT(Indiv, MUTAT1 + hapnr) : NULL;
      CondToInt(mutations, mutat);
      
      for (Uint index=0; index <= n_recomb; index++) {
	// printf("n_rec=%d %d\n", n_recomb, index);
	
	double cur_recomb = recombi[index];
	Uint cur_chrom; // starting with number 0

	//	if (zaehler >= 11424) printf("index=%d %d\n", index, n_recomb);

	if (cur_recomb <= first_snp_position) {
	  from = 0;
	  continue;
	}
	
	for (cur_chrom=0; cur_chrom<n_chrom; cur_chrom++)
	  if (length_total[cur_chrom + 1] >= cur_recomb) break;
	if (equidist) {
	  if (cur_recomb > last_snp_position) till = n_snps - 1;
	  else {
	    double halfbetween = REAL(VECTOR_ELT(position, cur_chrom))[0];
	    till = totalsnps_prec_chroms[cur_chrom] +
	      CEIL((cur_recomb - length_prec_chroms[cur_chrom] - halfbetween) /
		   (2.0 * halfbetween)) - 1;
	  }
	} else {

	  //	  printf("cur_chrom =%d\n", cur_chrom);
	  Uint min = cur_chrom == 0 ? 0 : totalsnps_prec_chroms[cur_chrom],
	    max = cur_chrom == n_chrom ? n_chrom - 1
	    : totalsnps_prec_chroms[cur_chrom + 1] - 1;
	    till = getNrSNPposition(cur_recomb, snp_positions, min, max);
	    // printf("%d %d %d last=%d till=%d\n", cur_recomb, min, max, last_snp_position, till);
	}

	if (index == 0) {
	  from = till + 1;
	  assert(from == 0);
	  continue;
	}
	
	assert((int) toSNP != NA_INTEGER && toSNP <= MAXINT);
	//printf("%d %d %d %d\n", till, from, toSNP, fromSNP);
	if (till >= from && from <= toSNP && till >= fromSNP) {
	  Uint von = MAX(fromSNP, from),
	    bis = MIN(toSNP, till),
	    o_info[4];
 	  decodeOrigin(originint[index-1], o_info);
	  //
	  assert(bis < n_snps);
	  assert(von <= bis);

	  //  SEXP Orig_gen = VECTOR_ELT(breeding, o_info[ORIGIN_GENERATION]), 
	  SEXP Orig_gen = VECTOR_ELT(breeding,
				     gen_ref[o_info[ORIGIN_GENERATION]]), 
	    O2 = VECTOR_ELT(Orig_gen, o_info[ORIGIN_SEX]),
	    Orig = VECTOR_ELT(O2, o_info[ORIGIN_NR]);
	  Uint o_haplotype = o_info[ORIGIN_HAPLO];
	  SEXP coded_haplo = VECTOR_ELT(Orig, CODEDHAPLO);
	  if (n_snps <= CodesPerBlock*3)//Coding selbst + hinten & vorne PLatz
	    ERR2("number of of snps is %d, but should be more than %d so that the programme can be used.", n_snps, CodesPerBlock * 3);
	  if ((Uint) length(coded_haplo) == n_snps) {
	    coded_haplo = codeHaplo2(VECTOR_ELT(Orig, HAPLOTYP1),
				     VECTOR_ELT(Orig, HAPLOTYP2));
	    SET_VECTOR_ELT(Orig, CODEDHAPLO, coded_haplo);
	  }
 
	  //printf("von=%d bis=%d snps=%d index=%d recomb=%d g=%d s=%d n=%d h=%d\n", von, bis, n_snps, index, n_recomb,o_info[ORIGIN_GENERATION],o_info[ORIGIN_SEX],o_info[ORIGIN_NR], o_info[ORIGIN_HAPLO] );
	  //printf("here\n");
	  //BUG;

	  //if (zaehler >= 11424) printf("ok\n");
	  BlockType *o_haplo = (BlockType0*) Align(coded_haplo, ALIGN_CODED);
	  //	  if (zaehler >= 11424 || zaehler >= 7220) printf("ok %d!\n", zaehler);
	  //	  printf("her rre\n");


	  for (Uint op=von; op<=bis; ) {
	    //
	    
	    //  printf("haplo=%d index=%d von=%d %d %d\n", hapnr, index, von, bis, op);
	    // printf("%d %d %.50s\n", (op % CodesPerBlock), (q % CodesPerBlock),  do_select);
	    if ( (op % CodesPerBlock) == (q % CodesPerBlock) && !do_select) {
	      // if (zaehler >= 11000) printf(".");
	      // blockwise copying	      
	      Uint 
		bnr = q / CodesPerBlock,
		start = q % CodesPerBlock,
		o_rest = 1 + bis - op,
		rest = CodesPerBlock - start;
	      if (rest > o_rest) rest = o_rest;
	      BlockType o_code = o_haplo[op / CodesPerBlock];
	      ////	      printf("rest = %d\n", rest);

	      //  if (zaehler >= 11424)   printf("a\n");
	      
	      if (mutations) {
		//if (zaehler >= 11424) printf("mutat = %d\n", n_mutat);
		Uint end = op + rest; // excluded, R-Coding !!
		if (mutatint[pmutat] <= end) { // mutat hat R codierung !!
		  // insbesondere bei select koennen Mutation uebersprungen
		  // werden. Somit muss zunaechst pmutat gegebenfalls erhoeht
		  // werden:
		  while (pmutat < n_mutat && mutatint[pmutat] <= op) pmutat++;
		  while (pmutat < n_mutat && mutatint[pmutat] <= end) {
		    //		    printf("%d %d\n", pmutat, n_mutat);
		    Uint idx = (mutatint[pmutat] - 1) % CodesPerBlock;
		    o_code ^= BitMaskStart[idx] & BitMaskEnd[idx];
		    pmutat++;
		  }
		  mutations = pmutat < n_mutat;
		}
	      }
	      //  if (zaehler >= 11424) printf("xxxx\n");
	      o_code &= BitMaskStart[start] & BitMaskEnd[start + rest - 1] &
		NullEins[o_haplotype];
	      if (o_haplotype != hapnr) {
		if (o_haplotype > hapnr) o_code = SHR(o_code, 1);
		else o_code = SHL(o_code, 1);
	      }

	      if (false) 
		if (bis==von+rest-1 && (bis == n_snps-1 || bis == 2 * n_snps - 1)) {
		//showblock(o_code, 1);
		// printf("haplo=%d index=%d von=%d %d %d rest=%d\n", hapnr, index, von, bis, op, rest);
		  // printf("%d %d %.50s\n", (op % CodesPerBlock), (q % CodesPerBlock),  do_select);
	      }
	      //   if (zaehler >= 11424) printf("fast ende %d\n", zaehler);
	      haplo[bnr] |= o_code;
	      // printf("%d %ld\n", bnr, haplo[bnr]);
	      //	      showblock((Uint*) haplo + bnr, 1);
	      //  if (zaehler >= 11424) printf("ende\n");
	      op += rest;
	      p += rest; // p is now endpoint + 1 !
	      q += rest;
	    } else {
	      BUG;
	      // if (op % 1000 == 0) printf(";");
	      if (false)	      
		printf("von %u %u %ld %ld z=%d %d %d %d %d\n",//
		       von, bis, 
		       op % CodesPerBlock, q % CodesPerBlock, zaehler,
		       from, till, index, n_recomb);
	   // BUG;
	      // bitwise copying
	      if (do_select) {
		if (psel >= len_sel || cur_select != p) {
		  op++;
		  p++;
		  continue;
		}
		cur_select = Selectint[psel++];
	      }
	      Uint H = GetHaploCode((BlockType0 *) o_haplo, op, o_haplotype);
	      if (mutations) {
		// note mutat:R coding; op:C-Coding
		while (pmutat < n_mutat && mutatint[pmutat] <= op) pmutat++;
		if (mutatint[pmutat] == op + 1) {
		  H = 1- H;
		  pmutat++;
		}
		mutations = pmutat < n_mutat;
	      }
	      if (H) PutHaploOne(q, hapnr, haplo);
	      //printf("H = %d\n", H);
	      op++;
	      p++;
	      q++;
	    }
	  }// for loop through all snps between two recombination points
	} // if space between two recombination points is not empty
	from = till + 1;
      } // for n_recomb
      FREEint(mutat);
      FREEint(origin);
    } // hapnr 
  }

  FREEint(snps_perchrom);
  FREEint(Select);
  FREEint(Generation);
  FREEint(Sex); 
  FREEint(Nr); 

  
  return 0;
}

    
 SEXP computeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		 SEXP From_SNP, SEXP To_SNP, SEXP Select, SEXP Geno){  
  if (Information == R_NilValue) Init();
  bool geno = LOGICAL(Geno)[0];

  Uint
    individuals = length(Generation);
  if (!geno && individuals != 1)
    ERR("haplo type is computed only for a single individuum.");
  Uint snps = IcomputeSNPS(population, Generation, Sex, Nr,
			   INTEGER(From_SNP)[0],
			   INTEGER(To_SNP)[0], Select, NULL);
   SEXP Ans = create_codevector(geno ? (individuals == 1 ? GENO : GENOMATRIX)
			       : HAPLO, snps, individuals);    
  BlockType *ans = (BlockType0*) Align(Ans, ALIGN_COMPUTE);
  IcomputeSNPS(population, Generation, Sex, Nr, INTEGER(From_SNP)[0],
	       INTEGER(To_SNP)[0], Select, ans);
  if (geno) haplo2geno(ans, snps, individuals);
  
  if (PRESERVE) R_PreserveObject(Ans);
 
  return Ans;
 }



SEXP IsolveRelMat(Uint individuals, double *Atau, double tau,
		  double *Vec, double beta, bool matrix) {
  const char *info[3] = {"rest", "yhat", "rel.matrix"};
  double
     *A = NULL,
    *pA = NULL;

  if (tau <= 0) ERR("'tau' must be positive");

  Uint elmts = 2 + matrix,
    iP1 = individuals + 1,
    i2 = individuals * individuals;
  SEXP Ans, namevec;
  PROTECT(Ans = allocVector(VECSXP, elmts));
  PROTECT(namevec = allocVector(STRSXP, elmts));
  setAttrib(Ans, R_NamesSymbol, namevec);
  for (Uint k=0; k<elmts; k++) SET_STRING_ELT(namevec, k, mkChar(info[k]));
 
  if (matrix) {
    SEXP RA;
    PROTECT(RA = allocMatrix(REALSXP, individuals, individuals));
    SET_VECTOR_ELT(Ans, 2, RA);
    MEMCOPY(REAL(RA), Atau, sizeof(double) * i2);
    pA = REAL(RA);
  } else {
    A =(double *) MALLOC(sizeof(double) * i2),
    MEMCOPY(A, Atau, sizeof(double) * i2);
    pA = A;
  }
  for (Uint i=0; i<i2; i+=iP1) Atau[i] += tau;

  SEXP rest, yhat;
  PROTECT(rest = allocVector(REALSXP, individuals));
  SET_VECTOR_ELT(Ans, 0, rest);
  double *r = REAL(rest);
  MEMCOPY(r, Vec, sizeof(double) * individuals);
  int err;
  err = Ext_solvePosDef(Atau, individuals, true, r, 1, NULL, NULL);
  if (err != NOERROR) {
    errorstring_type errorstring;
    Ext_getErrorString(errorstring);
    ERR1("error occurred when solving the system (%.50s)", errorstring);
  }
  PROTECT(yhat = allocVector(REALSXP, individuals));
  SET_VECTOR_ELT(Ans, 1, yhat);
  double *y = REAL(yhat);
  xA(r, pA, individuals, individuals, y);
  FREE(A);
  for (Uint i=0; i<individuals; y[i++] += beta);
  UNPROTECT(4 + matrix);
  if (PRESERVE) R_PreserveObject(Ans);
   
  return(Ans);
}


SEXP compute(SEXP popul, SEXP Generation, SEXP Sex, SEXP Nr,
	     SEXP Tau, SEXP Vec, SEXP Betahat,
	     SEXP Select, SEXP Matrix_return) {

  if (Information == R_NilValue) Init();
  Uint snps = IcomputeSNPS(popul, Generation, Sex, Nr, 1, NA_INTEGER,
			   Select, NULL),
    individuals = length(Generation),
    mem = UnitsAlign(snps, individuals);
  int *G = (int *) CALLOC(mem, BytesPerUnit);
  if (G == NULL) ERR1("mem allocation (%d bytes)", mem * BytesPerUnit);
  BlockType *g = (BlockType*) algn(G);
  if ((Uint) length(Vec) != individuals) ERR("'vec' not of correct length");

  IcomputeSNPS(popul, Generation, Sex, Nr,
	       1, NA_INTEGER, // R coded indices !!
	       Select, g);
   
  Uint i2 = individuals * individuals;
  double
    *A = (double *) MALLOC(sizeof(double) * i2);
  if (A == NULL) ERR("mem allocation");
  CentredRelationshipMatrix(g, snps, individuals, A);
  assert(A[0] >= 0.0);
  FREE(G);

  SEXP Ans = IsolveRelMat(individuals, A, REAL(Tau)[0], REAL(Vec),
			  REAL(Betahat)[0], LOGICAL(Matrix_return)[0]);
  FREE(A);

  return Ans;
}


SEXP solveRelMat(SEXP A, SEXP Tau, SEXP Vec, SEXP Beta, SEXP Destroy) {
  Uint individuals = nrows(A),
      bytes = individuals * individuals * sizeof(double);
  bool destroy = LOGICAL(Destroy)[0];
  double *a = NULL;

  if (!destroy) {
    a = (double *) MALLOC(bytes);
    MEMCOPY(a, REAL(A), bytes);
  }
  SEXP Ans = IsolveRelMat(individuals, destroy ? REAL(A) : a,
			  REAL(Tau)[0], REAL(Vec), REAL(Beta)[0], false);
  FREE(a);

  return Ans;
}
