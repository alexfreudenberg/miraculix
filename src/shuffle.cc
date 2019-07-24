
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

#define BitsPerCode 2L
#define PseudoSSE 1     // in case IntrinsicBase does not return anything better

#include "IntrinsicsBase.h"
#include "intrinsics.h"
#include <General_utils.h>
#include "xport_import.h"
#include "options.h"
#include "MX.h"
#include "haplogeno.h"
#include "Haplo.h"

#define checking true
#define INLINE inline


#if defined DO_FLOAT
#define real float
#define Real Float
#define ZEROREAL ZEROFLOAT 
#define INT2REAL INT2FLOAT 
#define MULTREAL MULTFLOAT 
#define ADDREAL ADDFLOAT  
#define SUBREAL SUBFLOAT
#if defined ANY_SIMD
#define REALVALUE .f
#define IR(X, Y) INT2FLOAT(f.u128[X], d.u128[Y]);
#else
#define REALVALUE 
#define IR(X, Y) INT2FLOAT(f, d);
#endif // not ANY_SIMD
#else
#define real double
#define Real Double
#define ZEROREAL ZERODOUBLE 
#define INT2REAL INT2DOUBLE 
#define MULTREAL MULTDOUBLE 
#define ADDREAL ADDDOUBLE  
#define SUBREAL SUBDOUBLE  
#if defined ANY_SIMD
#define REALVALUE .d
#define IR(X, Y) INT2DOUBLE(f.d128[X], d.m64[Y]);
#else
#define REALVALUE 
#define IR(X, Y) INT2DOUBLE(f, d.m64[Y]);
#endif // not ANY_SIMD

#endif // not DO_FLOAT


int xxx = CodesPerBlock;

//
ALL_INLINER

/*
Uint xxxx() { return Blocks(5); }


Uint inline *algn(int *X) {assert(algn_general(X, BytesPerBlock)>=(uintptr_t)X); return (Uint *) algn_general(X, BytesPerBlock); } 
									

Uint inline Blocks(Uint X) { printf("X=%d zzzCPB=%ld %ld %ld\n", X, CodesPerBlock, CodesPerUnit, UnitsPerBlock); return 1L + (X - 1L) / CodesPerBlock; } 

Uint inline XBlocks(Uint X) { printf("X=%d aaaCPB=%ld %ld %ld\n", X, CodesPerBlock, CodesPerUnit, UnitsPerBlock); return 1L + (X - 1L) / CodesPerBlock; };

//#define Units(snps) (Blocks(snps) * UnitsPerBlock)

real inline *algnrl(real *X) {assert(algn_general(X, BytesPerBlock)>=(uintptr_t)X); return (real *)  algn_general(X, BytesPerBlock); } 
 
Uint inline RealAlign(Uint X) { return BytesPerBlock / sizeof(real) + (1L +  (X - 1L) / CodesPerBlock) * CodesPerBlock; }


//g++ -fno-omit-frame-pointer  -fsanitize=undefined,address -fno-omit-frame-pointer  -std=gnu++11 -I"/usr/local/lib/R/include" -DNDEBUG  -I"/home/schlather/R/x86_64-pc-linux-gnu-library/3.6/RandomFieldsUtils/include" -I/usr/local/include  -mavx  -fpic  -I/usr/lib/jvm/java-8-openjdk-amd64/include -I/usr/lib/jvm/java-8-openjdk-amd64/include/linux/ -fno-omit-frame-pointer  -fsanitize=undefined,address -fno-omit-frame-pointer    -c xport_import.cc -o xport_import.o // i, 32, falsch!

//g++ -I/usr/local/lib/R/include -I/usr/local/include -I/usr/lib/gcc/x86_64-linux-gnu/5/include/ -O2 -Wformat -Wformat-security -Wformat-y2k -Wuninitialized -Wall -Wextra -Wshadow -Wpointer-arith -pedantic -g -Wswitch-default  -fopenmp -fPIC -pthread -pipe -Wformat-overflow=2  -march=native -Wno-unused-variable -Wno-unused-function  -fno-omit-frame-pointer  -fsanitize=undefined,address -fno-omit-frame-pointer    -march=native -c miraculix/src/options.cc -o miraculix.CRAN/src/options.o // rr, 64 OK!

// clang++  -std=gnu++11 -I"/usr/local/lib/R/include" -DNDEBUG  -I"/home/schlather/R/x86_64-pc-linux-gnu-library/3.6/RandomFieldsUtils/include" -I/usr/local/include -fopenmp -mavx  -fpic  -I/usr/lib/jvm/java-8-openjdk-amd64/include -I/usr/lib/jvm/java-8-openjdk-amd64/include/linux/     -c options.cc -o options.o





*/


int yyy = CodesPerBlock;



// [1 + (X - 1) / UnitsPerBlock] * UnitsPerBlock + (UnitsPerBlock - 1)
// last part: to make finally sure that first used element of R-vector
// is aligned (maximum distance to 0-elemtn: UnitsPerBlock - 1)
// first part: then memory is read in blocks; so next number of Bytes that is
// is divisible by Unitsperblock * Bytesperunit is taken:
// 1 + (X - 1) / UnitsPerBlock]



BlockType and2scalar, xor2scalar, LH
#ifndef ANY_SIMD
  , SHUFFLEMASK1, SHUFFLEMASK2
#endif
  ;

#if defined SSSE3 or defined AVX2
#define ShuffleOK 1
#elif defined ShuffleOK
#undefine ShuffleOK
#endif


Uint CodesPerBlockShuffle() { return CodesPerBlock; }
Uint UnitsPerIndivShuffle(Uint snps) { return Blocks(snps) * UnitsPerBlock; }


void showblock(Uint *X, Uint snps) { 
  Uint *x = X,
    bits = snps * BitsPerCode,
    endfor = 1L + (Uint) ((bits - 1L) / BitsPerBlock);
  
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
void showblock(BlockType0 X, Uint blocks) {
  showblock((Uint *) &X, blocks * CodesPerBlock);
}

void showblock(Uint *X, Uint snps, Uint individuals) { 
  for (Uint i=0; i<individuals; i++) {
#ifdef SCHLATHERS_MACHINE    
    // PRINTF("%ld\n", (uint ptr_t) (X + i * Blocks(snps) * UnitsPerBlock));
#endif 
    showblock(X + i * Blocks(snps) * UnitsPerBlock, snps);
  }
}

Uint permut_finv_g_h[CodesPerBlock];
bool shuffleNotInit = true;
void InitShuffle() {
#if not defined ShuffleOK
  // die PseudoSSE scheint so furchtbar langsam zu sein, dass
  // es sich bis auf die Haplo-Sachen nicht lohnt die Fehler rauszuholen
  // to do ?! Fehler rausholen?
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert(BytesPerBlock == sizeof(BlockType0));
  if (!shuffleNotInit ) BUG;
 
//                   0  1  2  3  4  5  6  7   8  9 10 11 12 13 14 15
//#define andshuffle 0, 1, 4, 0, 1, 2, 5, 0,  4, 5, 8, 0, 0, 0, 0, 0
//#define xorshuffle 0, 0, 0, 2, 0, 0, 0, 2,  0, 0, 0, 2, 2, 2, 2, 4
//    SETREV8(as[i], andshuffle);
//    SETREV8(xs[i], xorshuffle);
    SETREV8(and2scalar, 0, 1, 4, 0, 1, 2, 5, 0,  4, 5, 8, 0, 0, 0, 0, 0);
    SETREV8(xor2scalar, 0, 0, 0, 2, 0, 0, 0, 2,  0, 0, 0, 2, 2, 2, 2, 4);

  SET32(LH, 0x0F0F0F0F);
#ifndef ANY_SIMD
  SET32(SHUFFLEMASK1, 0x01010101);
  SET32(SHUFFLEMASK2, 0x02020202);
#endif  

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
#endif    
}



void assert_shuffle() {
  Uint method = GLOBAL.genetics.method;
#ifdef PseudoSSEconfirmed
    ERR("snpcoding=Shuffle needs at least SSSE3");
#endif    
  // die PseudoSSE scheint so furchtbar langsam zu sein, dass
  // es sich bis auf die Haplo-Sachen nicht lohnt die Fehler rauszuholen
  // to do ?! Fehler rausholen?
  if (method != Shuffle) { // will work if aligned: ( && method != TwoBit)
    PRINTF("NOTE: 'RFoptions(snpcoding=Shuffle)' is necessary for this function.\nSo, 'RFoptions' is set accordingly by the package.\n");
    GLOBAL.genetics.method = Shuffle;
    // ERR("Set 'RFoptions(snpcoding=Shuffle)' first.");
  }
  if (shuffleNotInit) InitShuffle();
}

void assert_shuffle(SEXP M) {
  Uint method = MethodOf(M);
  if (method != Shuffle && method != Haplo)
    ERR("Function only available for method 'Shuffle'\n");
  if (shuffleNotInit) InitShuffle();
}



SEXP create_codevector(Uint snps, Uint individuals) {
  assert_shuffle();  
  Uint  memPerIndiv = Blocks(snps) * UnitsPerBlock;
  SEXP Code = CreateEmptyCodeVector(snps, individuals, memPerIndiv);
  if (PRESERVE) {BUG; R_PreserveObject(Code);}
  return(Code);
}


Ulong sumGenoShuffle(Uint *S, Uint snps, Uint individuals) {
 Ulong ans = 0;
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  Ulong blocks = (Ulong) Blocks(snps) * individuals;

#define x31  (255L / 8L) // 2 : max{0,1,2}
  //                         4 : 2-bit-coding in 1 Byte
#define BytesPer64  (64 / BitsPerByte)
 //#define intgroup_size (2147483648L / (x31 * 8L * BytesPerUnit)) // 2^31 / x31
#define intgroup_size 3 // 2^31 / x31
#define intgroup_sizeM1 (intgroup_size - 1L)
  Ulong
    intgroupsM1 = (blocks - 1L) / (x31 * intgroup_size),
    sizex31 = intgroup_size * x31;
  //    groupsM1 = (blocks - 1L) / x31;
  BlockType first2bits8, first8bits;
    // NullEins = SET32(0x55555555),
  SET32(first2bits8, 0x03030303);
  SET32(first8bits, 0x000000FF);
  BlockUnitType sum32;
 
  for (Ulong X=0; X<=intgroupsM1; X++) { // this loop needed only for huge
    // data sizes of more than 50000 individuals and 10^5 SNPS each
    ZERO(sum32 VI);
    BlockType *S_int = (BlockType0 *) S + X * sizex31 ;
    Ulong groupsM1 = X < intgroupsM1 ? intgroup_sizeM1
			: (blocks - 1L) / x31 - intgroupsM1 * intgroup_size;
    for (Ulong i=0; i<=groupsM1; i++) {
      BlockType sum,
	*cgm = S_int + i * x31;
      ZERO(sum);
      Uint endx = i < groupsM1 || X < intgroupsM1 ? x31 :
	blocks - intgroupsM1 * sizex31 - x31 * groupsM1;
      for (Uint x=0; x<endx; x++, cgm++) {
	BlockType d,
	  c = *cgm;
	//	d = SHR32(c, 1);
	//	d = AND(NullEins, d);
	//	c = XOR(c, d);      
	for (Uint j=0; j<CodesPerByte; j++) { // for-schleife aufloesen ?
	  AND(d, c, first2bits8);
	  ADD8(sum, sum, d);
	  SHR32(c, c, 2);
	}
      }
      for (Uint j=0; j<BytesPerUnit; j++) {
	BlockType L1;
	AND(L1, sum, first8bits);
	ADD32(sum32 VI, sum32 VI, L1); 
	SHR32(sum, sum, BitsPerByte);
      }
    } // for groups

    for (Uint i=0; i<UnitsPerBlock; i++) {
      ans += (Ulong) sum32.u32[i];
    }
  }

  //printf("sumgeno = %lu\n", ans);
#endif
  return ans;
}


void ReUseAsShuffle(SEXP Code) {
  Uint *info = GetInfo(Code),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS];
  Ulong memInUnitsPerIndiv = Blocks(snps) * UnitsPerBlock;  
  InternReUseAs(Code, Shuffle, snps, individuals, memInUnitsPerIndiv, false);  
}
	
void haplo2genoShuffle(Uint *X, Uint snps, Uint individuals, Uint *Ans) {
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  Uint blocks = Blocks(snps);
  BlockType EN, NE;
  SET32(EN, 0xAAAAAAAA);
  SET32(NE, 0x55555555); 

  for (Uint j=0; j<individuals; j++) {
    Uint jb = j * blocks;
    BlockType *x = ((BlockType0 *) X)  + jb,
      *a = ((BlockType0 *) Ans)  + jb;
    for (Uint b=0; b<blocks; b++) {
      BlockType y,z,dummy;
      AND(y, x[b], EN);
      AND(z, x[b], NE);
      SHL32(dummy, z, 1L);
      AND(dummy, y, dummy);
      SHR32(y, y, 1);
      XOR(y, y, z);
      
      OR(a[b], dummy, y); 
      // a[b] = OR(SHR32(y, 1L), OR(z, AND(y, SHL32(z, 1L))));
    }
  }
#endif
}





void zeroNthGenoShuffle(SEXP CM, SEXP NN) { // 0.45 bei 25000 x 25000
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert_shuffle(CM);
  Uint
    first2bits32 = 0x00000003,
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    blocks = Blocks(snps),
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENOMATRIX) ERR("not a genomicmatrix coding");
  
  for (Ulong i=0; i<individuals; i++) {    
    BlockType *C1 = ((BlockType0 *) C0) + i * blocks;
    for (Uint ni=0; ni<lenN; ni++) {
      Uint j=N[ni] / CodesPerBlock,
	jj = N[ni] % CodesPerBlock,
	u = jj / UnitsPerBlock,
	k = jj % UnitsPerBlock;     
      BlockType *C = C1 + j;
      ((BlockUnitType*) C)->u32[k] &= ~(first2bits32 << (BitsPerCode * u));
    }
  }
#endif
}




SEXP get_matrixN_shuffle(SEXP CM, SEXP NN) { // 0.45 bei 25000 x 25000
  SEXP Ans = R_NilValue;
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert_shuffle(CM);
  //  BlockType first2bits32 = SET32(0x00000003);
  Uint
    first2bits32 = 0x00000003,
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    blocks = Blocks(snps),
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENOMATRIX) ERR("not a genoimatrix coding");
  
  if (lenN ==1) PROTECT(Ans = allocVector(INTSXP, individuals));
  else PROTECT(Ans = allocMatrix(INTSXP, lenN, individuals));
  Uint *M = (Uint*) INTEGER(Ans);

  for (Ulong i=0; i<individuals; i++) {    
    Uint *mm = M + i * lenN;
    for (Uint ni=0; ni<lenN; ni++, mm++) {
      Uint nni = N[ni];
      Uint j=nni / CodesPerBlock,
	jj = nni % CodesPerBlock,
	u = jj / UnitsPerBlock,
	k = jj % UnitsPerBlock;     
      BlockType c = ((BlockType0 *) C0)[i * blocks + j];
      BlockUnitType d;
      d VI = c;
      // d VI = AND(SHR32(c, BitsPerCode * u), first2bits32);  *mm = d.u32[k];
      *mm = (d.u32[k] >> (BitsPerCode * u)) & first2bits32;
    }
  }
 
  UNPROTECT(1);
#endif
  
  return Ans;
}


void vectorGenoKahan(Real *V, BlockType0 *CM, Uint snps, Uint individuals,
		     double *ans) {
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  Uint blocks = Blocks(snps);
  BlockType first2bits32;
  // NullEins = SET32(0x55555555),
  SET32(first2bits32, 0x00000003);


#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint i=0; i<individuals; i++) {
    BlockType
      *cm = CM + i * blocks;
    Real sum, corr,
      *v = V;
    ZEROREAL(sum);
    ZEROREAL(corr);
    for (Uint b=0; b<blocks; b++, cm++) {
      BlockUnitType d,f;
      BlockType c = *cm;
      //   d VI = SHR32(c, 1);
      // d VI = AND(NullEins, d VI);      
      // c = XOR(c, d VI);      
      for (Uint j=0; j<CodesPerUnit; j++) {
	AND(d VI, c, first2bits32);
	Real  dummy, fcorr, roughsum;

#define VGK					\
	MULTREAL(f REALVALUE, *v, f REALVALUE);			\
	SUBREAL(fcorr, f REALVALUE, corr);		\
	ADDREAL(roughsum, sum, fcorr);		\
	SUBREAL(dummy, roughsum, sum);		\
	SUBREAL(corr, dummy, fcorr);		\
	sum = roughsum

	
	IR(0,0);
#ifdef AVX2
	IR(1,1);
#endif	
	VGK;

#if not defined DO_FLOAT
#ifdef AVX2
	IR(0,2); IR(1,3);
#else
	IR(0,1);
#endif
	v++;
	VGK;
#endif
	
	v++;
	SHR32(c, c, 2L);
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

  #endif
}


void vectorGenoIntern(Real *V, BlockType0 *CM, Uint snps, Uint individuals,
		      double *ans) {
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else

  if (GLOBAL_UTILS->basic.kahanCorrection) {
    vectorGenoKahan(V, CM, snps, individuals, ans);
    return;
  }
  Uint blocks = Blocks(snps);
  BlockType    first2bits32;
  // NullEins = SET32(0x55555555),
  SET32(first2bits32, 0x00000003);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Uint i=0; i<individuals; i++) {
    BlockType
      *cm = CM + i * blocks;
    Real sum,
      *v = V;
    ZEROREAL(sum);
    for (Uint b=0; b<blocks; b++, cm++) {
      BlockUnitType d,f;
      BlockType c = *cm;
      //      d VI = SHR32(c, 1);
      //d VI = AND(NullEins, d VI);      
      // c = XOR(c, d VI);      
      for (Uint j=0; j<CodesPerUnit; j++) {
	AND(d VI, c, first2bits32);

#define VGI							\
	MULTREAL(f REALVALUE, *v, f REALVALUE);			\
	ADDREAL(sum, sum, f REALVALUE);				\
	
	IR(0, 0);
#ifdef AVX2
	IR(1, 1);
#endif	
	VGI;

#if not defined DO_FLOAT
#ifdef AVX2
	IR(0, 2); IR(1, 3);
#else
	IR(0, 1);
#endif
	v++;
	VGI;
#endif

	
	v++;
	SHR32(c, c, 2L);
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

  #endif
}

SEXP vectorGeno(SEXP V, SEXP Z) {
  SEXP Ans = R_NilValue;
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert_shuffle(Z);
  Uint
    *info = GetInfo(Z),
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
    Rint *vp = LOGICAL(V);
    for (Uint i=0; i<len; i++) v[i] = (real) vp[i];
    break;
  }
  default : ERR("'V' of unknown type.");
  }

  PROTECT(Ans = allocVector(REALSXP, individuals));
  double *ans = REAL(Ans);
  for (Uint i=0, endfor= individuals; i<endfor; i++) ans[i] = 0.0; // noetig?
  vectorGenoIntern((Real *) v, (BlockType0 *) Align(Z, ALIGN_VECTOR),
		   snps, individuals, REAL(Ans));

  Free(v0);
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
#endif  
  return(Ans);
}

//entering genoVI CpB=64
//X=39900 aaaCPB=64 16 4
//snps=39900 CpB=64  blocks = 624 xxx=64 64

//entering genoVI CpB=64
//X=39900 zzzCPB=64 16 4
//snps=39900 CpB=64  blocks = 624 xxx=64 64

//entering genoVI CpB=64
//X=39900 iii CPB=32 =  16 * 2
//snps=39900 CpB=64  blocks = 1247 xxx=64 64


void genoVectorIntern(BlockType0 *CM, double *V, Uint snps, Uint individuals,
		      double *ans) {
 #if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else

  if (GLOBAL_UTILS->basic.kahanCorrection) {
    //    genoVectorKahan(V, CM, snps, individuals, ans);
    // return;
  }
  printf("entering genoVI CpB=%d\n", CodesPerBlock);
  Uint blocks = Blocks(snps),
    blocksM1 = blocks - 1;

  printf("snps=%d CpB=%d  blocks = %d xxx=%d %d\n", snps,  CodesPerBlock,blocks,xxx, yyy); BUG;
  
  double *end_ans = ans + snps;
  BlockType c,first2bits32;
  SET32(first2bits32, 0x00000003);
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

      
      Uint endfor = b < blocksM1 ? CodesPerUnit
		: (snps - blocksM1 * CodesPerBlock) / UnitsPerBlock;
      //      printf("endfor = %d blocksM1=%d upb=%d cp=%d\n", endfor, blocksM1, UnitsPerBlock, CodesPerUnit);


      printf("snps=%d, blocksM1=%d endf=%d CpB=%d BpC=%d\n", snps, blocksM1, endfor, CodesPerBlock, BitsPerCode);
#if defined AVX2
      printf("AVX2\n");
#else
        printf("KEIN AVX2\n");    
#endif 
    
	// a+=4 * endfor; continue; printf("Rd corrigier\n");

      for (Uint j=0; j<endfor; j++) {
	BlockUnitType d;
	AND(d VI, c, first2bits32);
	a[0] += factor[d.u32[0]]; 
	a[1] += factor[d.u32[1]];
	a[2] += factor[d.u32[2]];
	a[3] += factor[d.u32[3]];
#if defined AVX2
	a[4] += factor[d.u32[4]]; 
	a[5] += factor[d.u32[5]];
	a[6] += factor[d.u32[6]];
	a[7] += factor[d.u32[7]];
	a += 8;

	/* todo speed might be improved using
mm256_bitshuffle_epi64_mask
_mm_blendv_epi8 (__m128i a, __m128i b, __m128i mask)
Rint _mm_cmpestra (__m128i a, Rint la, __m128i b, Rint lb, const Rint imm8)
__m128i _mm_hadd_epi32 (__m128i a, __m128i b)
__m64 _mm_maddubs_pi16 (__m64 a, __m64 b)
void _mm_maskmove_si64 (__m64 a, __m64 mask, char* mem_addr)
Rint _mm_movemask_ps (__m128 a)
Rint _m_pmovmskb (__m64 a)
Rint _mm256_movemask_pd (__m256d a)
_mm256_mask_bitshuffle_epi64_mask (__mmask32 k2, __m256i b, __m256i c)
__m128d _mm_mask_permute_pd (__m128d src, __mmask8 k, __m128d a, const Rint imm8)
__m256d _mm256_permute_pd (__m256d a, Rint imm8)
__m256d _mm256_permute2f128_pd (__m256d a, __m256d b, Rint imm8)
__m256d _mm256_permute4x64_pd (__m256d a, const Rint imm8)
	*/
	
#else
	a += 4;
#endif 
	SHR32(c, c, BitsPerCode);
      }
    }

    BlockUnitType d;
    AND(d VI, c, first2bits32);

    printf("end_ans - a = %ld %d \n", end_ans-a, UnitsPerBlock);
    
    assert(end_ans >= a && end_ans - a < UnitsPerBlock);
    for(Uint k=0; a < end_ans; a++) {
      //   printf("%f %f\n", factor[d.u32[k]], factor[d.u32[k+4]]);
      *a += factor[d.u32[k++]];
    }
  }
#endif
}


SEXP genoVector(SEXP Z, SEXP V) {
  SEXP Ans = R_NilValue;
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert_shuffle(Z);
  Uint
    *info = GetInfo(Z),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    len = length(V);
  if (len != individuals) ERR("vector 'V' not of correct length");
  PROTECT(Ans = allocVector(REALSXP, snps));
  genoVectorIntern((BlockType0 *) Align(Z, ALIGN_VECTOR), REAL(V), 
 		   snps, individuals, REAL(Ans));
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
#endif  
  return(Ans);
}


#if defined ShuffleOK
BlockType INLINE VECTORADD(BlockType0 xx, BlockType0 yy) {
  BlockType z, sum, y;  
  AND(z, xx, yy);
  AND(y, z, LH); // 0.25 von 4
  SHUFFLE8(sum, and2scalar, y); // 0.5 von 4

  SHR32(y, z, 4L);
  AND(y, y, LH); // 
  SHUFFLE8(y, and2scalar, y);// 
  ADD8(sum, sum, y); // 
  //   showblock(sum, CodesPerBlock);

#ifndef OLD_VECTORADD
  XOR(z, xx, yy);
  AND(y, z, LH);
  SHUFFLE8(y, xor2scalar, y);
  ADD8(sum, sum, y);

  SHR32(y, z, 4L);
  AND(y, y, LH);
  SHUFFLE8(y, xor2scalar, y);
#else  // FOLGENDES IST NICHT SCHNELLER !!!
  XOR(z, xx, yy);
  SHR32(y, z, 1L);
  AND(z, y, z);
  SHR32(y, z, 3L);
  OR(z, z, y);
  AND(y, z, LH);
  SHUFFLE8(y, xor2scalar, y);
#endif

  ADD8(sum, sum, y);
  return sum;
}

/*
#else  // no SIMD available
BlockType INLINE VECTORADD(BlockType0 xx, BlockType0 yy) {
  // 1 Operation more because of highest bit
  // Idee: die werte direkt eintragen, statt ueber die Tabelle nachschauen:
  BlockType z, sum, y, x;
  // 2 * 1 oder 2 * 1
  OR(x, xx, yy);
  SHL32(z, x, 1L);
  AND(x, x, z);
  AND(x, x, SHUFFLEMASK2); // ergb
    
  AND(z, xx, yy);
  AND(y, z, SHUFFLEMASK1); // 1 * 1
  OR(x, x, y);     // "addition" zu "2 * 1"
  
  AND(z, z, SHUFFLEMASK2); // 2 * 2;

  // globaler shift left leider nicht moeglich da oberstes Bit runterfaellt
  // somit werden die beiden (*) Zeilen zusaetzlich benoetigt
  // ohne die beiden unteren Zeilen hat code nur 14 (statt VECTADD 14)
  // Operationen mit einem Zeitgewinn von somit 7% und einem damit notwendigern
  // Speichmehrbedarf von 1.5% bei SSSE3 und 0.8% bei AVX
  SHL32(y, z, 1L);
  OR(y, z, y);
  AND(sum, y, LH);

  SHR32(y, x, 4L);
  SHR32(z, z, 3L); // (*)
  OR(y, y, z);   // (*)
  AND(y, y, LH);

  ADD8(sum, sum, y);
  
  return sum;
}
*/

#endif //  INLINE VECTORADD


#if defined ShuffleOK
Uint INLINE ADDUPVECTOR(BlockType0 x) {
  BlockType zero;
  ZERO(zero);
#if defined AVX2
  BlockUnitType vsum;
  SAD8(vsum VI, x, zero);  
  return vsum.u32[0]+ vsum.u32[2] + vsum.u32[4] + vsum.u32[6];
#else 
  BlockType vsum;
  SAD8(vsum, x, zero); //
  Uint L1, L2;
  EXTRACT16(L1, vsum, 0);
  EXTRACT16(L2, vsum, 4);
  return L1 + L2;
#endif
}
#endif
 
Uint scalar012(BlockType0 *x, BlockType0 *y, Uint blocks) {
  Uint ans = 0;
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
 #define x15  (255L / 16L) // 4 : max=2^2
  //                         4 : 2-bit-coding in 1 Byte                    
  Uint
    bundles = blocks / x15;

  for (Uint i=0; i<bundles; i++) { 
    BlockType L1, L2, sum,
      *xx = x + i * x15,
      *yy = y + i * x15;
    ZERO(sum);
#if x15 == 15L
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));

    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));

    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));

#else
    for (Uint j=0; j<x15; j++) {
      LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    }
#endif
    ans += ADDUPVECTOR(sum); 
  }
  BlockType sum ,
    *endx = x + blocks,
    *xx = x + bundles * x15,
    *yy = y + bundles * x15;
  ZERO(sum);
  for (; xx < endx; xx++, yy++) ADD8(sum, sum, VECTORADD(*xx, *yy));
  ans += ADDUPVECTOR(sum);
#endif
  return ans;
}




SEXP get_matrixshuffle(SEXP CM) { // 0.45 bei 25000 x 25000
  SEXP Ans = R_NilValue;
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert_shuffle(CM);
  BlockType first2bits32;
  SET32(first2bits32, 0x00000003);
  Uint 
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    bits = snps * BitsPerCode,
    fullBlocks = bits / BitsPerBlock,
    blocks = Blocks(snps),    
    *C0 = Align(CM, ALIGN_HAPLO);
  if (info[WHAT] != GENOMATRIX) ERR("not a geno coding");
   
  PROTECT(Ans = allocMatrix(INTSXP, snps, individuals));
  Uint *M = (Uint*) INTEGER(Ans);
  for (Uint i=0, endfor= snps; i<endfor; i++) M[i] = 0; // noetig?

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Ulong i=0; i<individuals; i++) {
    Uint *mm = M + i * snps;
    BlockType *C = (BlockType0 *) C0 + i * blocks;
    for (Uint j=0; j<fullBlocks; j++, C++) {
      BlockType c = *C;
      for (Uint u=0; u<CodesPerUnit; u++, mm+=UnitsPerBlock) {	  
	BlockUnitType d;
	AND(d VI, c, first2bits32);
#define MD(i) mm[i] = d.u32[i];
#if UnitsPerBlock == 4
	MD(0); MD(1); MD(2); MD(3);
#elif UnitsPerBlock == 8
	MD(0); MD(1); MD(2); MD(3); MD(4); MD(5); MD(6); MD(7);
#else	    
	for (Uint k=0; k<UnitsPerBlock; k++) MD(k);
#endif	    
	SHR32(c, c, BitsPerCode);
      }
    }
    if (fullBlocks < blocks) {
      Uint *end = M + (i + 1L) * snps;
      BlockType c = ((BlockType0 *) C0)[i * blocks + fullBlocks];
      for (Uint u=0; u<CodesPerUnit && mm < end; u++) {
	BlockUnitType d;
	AND(d VI, c, first2bits32);
	for (Uint k=0; k<UnitsPerBlock && mm < end; k++, mm++)  *mm = d.u32[k];
	SHR32(c, c, BitsPerCode);
      } 
    }
  }
 
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
#endif  
  return Ans;
}



SEXP matrix_start_shuffle(Uint individuals, Uint snps, SEXP file) {  
  SEXP Code = create_codevector(snps, individuals);
  start_info(Code, file, UnitsPerBlock, CodesPerBlock);
  return(Code);
}
 



void codeInnerGeno(Uint *M, Uint snps, Uint *code) {
 #if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
 // Vektorisieren ?!
  Uint j,
    bits = snps * BitsPerCode,
    //fullBlocks = (bits - 1)  / BitsPerBlock;
    fullBlocks =  bits / BitsPerBlock; // 30.5.2019

 
  if (false && checking) {
    for (j=0; j<fullBlocks; j++) {
      Uint *C = code + j * UnitsPerBlock,
	*mm = M + j * CodesPerBlock;
      for (Uint u=0; u<CodesPerUnit; u++) {
	Uint shift = u * 2;
	for (Uint k=0; k<UnitsPerBlock; k++, mm++) {
	  Uint m = *mm;
	  C[k] |= m << shift; 
	}
      }
    }
  } else {
    BlockType0 *Code = (BlockType0 *) code;
    for (j=0; j<fullBlocks; j++) {
      BlockType0 C;
      ZERO(C);
      Uint *mm =  (M + j * CodesPerBlock);
      for (Uint u=0; u<CodesPerUnit; u++, mm+=UnitsPerBlock) {
	Uint shift = u * 2;
	BlockType0 m;
	LOADU(m, (BlockType0 *) mm);
	SHL32(m, m, shift);
	OR(C, C, m);
      }
      Code[j] = C;
    }
  }
  
  Uint snps_done = fullBlocks * CodesPerBlock;
  if (snps_done < snps) {
    Uint *end = M + snps,
      *C = code + UnitsPerBlock * fullBlocks;
    M += snps_done;
    for (Uint u=0; u<CodesPerUnit; u++) {
      Uint shift = u * 2;
      for (Uint k=0; k<UnitsPerBlock && M < end; k++, M++) {
	Uint m = *M;
	//	if (checking && (m < 0 || m > 2)) ERR("wrong values for SNP matrix");
	C[k] |= m << shift;
      }
      if (M == end) break;
    }
  }
 
#endif
  
}



void matrix_shuffle(Uint *M, Uint start_individual, Uint end_individual, 
		    Uint start_snp, Uint end_snp, Uint Mnrow,
		    SEXP Ans,
		    double VARIABLE_IS_NOT_USED *G) {
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  if (start_snp % CodesPerBlock != 0) BUG;
  Uint
    *info = GetInfo(Ans),
    //    individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    cur_snps = end_snp - start_snp,
    *pM = M,
    *endM = M + Mnrow * (end_individual - start_individual),
    blocks = Blocks(snps);
  Ulong units = (Ulong) blocks * UnitsPerBlock;
  Uint *ans = ((Uint*) Align(Ans, ALIGN_HAPLOGENO)
	       ) + start_snp / CodesPerUnit + start_individual * units;    

  //printf("Mnrow = %d %d\n", Mnrow, end_individual - start_individual);
  for (; pM<endM; pM += Mnrow, ans += units) {
    codeInnerGeno(pM, cur_snps,ans);
  }
  #endif
}


SEXP matrix_coding_shuffle(Uint *M, Uint snps, Uint individuals) {  
  SEXP Code = matrix_start_shuffle(individuals, snps, R_NilValue);
  matrix_shuffle(M, 0, individuals, 0, snps, snps, Code, NULL);    
  return Code;
}



void matrixshuffle_mult(Uint* CGM, Uint snps, Uint individuals, double *ans) {
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  Uint blocks = Blocks(snps);
    
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 16)  
#endif
  for (Uint i=0; i<individuals; i++) {
    for (Uint j=i; j<individuals; j++) {
      Uint Mij = scalar012(((BlockType0 *) CGM) + i * blocks,
			    ((BlockType0 *) CGM) + j * blocks, blocks);
      ans[i + j * individuals] = ans[i * individuals + j] = (double) Mij;
    }
  }
  #endif
}
 



Uint *AlignShuffle(SEXP Code, Uint nr, bool test) {
  return AlignTest(Code, nr, test); }


real *Pd0 = NULL; 					   
void allele_freq2(BlockType0* CGM, Uint snps, Uint individuals,
		  double *sumP, double *sumP2) {
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  // NOTE! twice allel frequency (indicated by the 2 at the end) 
  double n =  (double) individuals;
  if (snps > 536870912) ERR("too many SNPs");
  Uint blocks = Blocks(snps);

  
#define x31  (255L / 8L) // 2 : max{0,1,2}
  //                         4 : 2-bit-coding in 1 Byte                    

  Uint
    byteblocks = (snps - 1L) / BytesPerBlock + CodesPerByte,
    alignedByteBlocks = byteblocks + 1L,
    units = byteblocks * BytesPerBlock,
    alignedUnits = units + UnitsPerBlock,
    // here: P= 2(p_i), not 2 p_i - 1 == expected values == means
    groupsM1 = (individuals - 1L) / x31; 
  Rint
    *Sums0 = (Rint*) MALLOC(alignedUnits * BytesPerUnit),      
    *sum0 = (Rint*) MALLOC(alignedByteBlocks * BytesPerBlock);// snps + zuschlag
  BlockType first8bits, first2bits8,
    *sum = (BlockType0*) algn(sum0),
    *Sums  = (BlockType0*) algn(Sums0);
  // NullEins = SET32(0x55555555),
  SET32(first8bits, 0x000000FF);
  SET32(first2bits8, 0x03030303);
  if (Sums == NULL || sum == NULL) ERR("allocation error");
  { BlockType zero;
    ZERO(zero);
    Uint endfor = units / UnitsPerBlock;
    for (Uint k=0; k<endfor; Sums[k++] = zero); }
  // parallel funktioniert nicht gut!
  for (Uint i=0; i<=groupsM1; i++) {
    BlockType *cgm = CGM + i * blocks * x31;
    Uint endx = i < groupsM1 ? x31 : individuals - x31 * groupsM1;
    {
      BlockType zero;
      ZERO(zero);
      for (Uint k=0; k < byteblocks; sum[k++] = zero);
    }
    for (Uint x=0; x<endx; x++) {
      BlockType *s = sum;
      for (Uint k=0; k<blocks; k++, cgm++) { // see permutation.tex (*)
	BlockType d,
	  c = *cgm;
       	//	d = SHR32(c, 1);
	//	d = AND(NullEins, d);
	//	c = XOR(c, d);      
	for (Uint j=0; j<CodesPerByte; j++, s++) { // for-schleife aufloesen?
	  AND(d, c, first2bits8);
	  ADD8(*s, *s, d);
	  SHR32(c, c, 2);
	}
      }
    }
    
    //    R < Simianer_sims_23_01.R --vanilla
	 
    BlockType *pS = (BlockType0*) Sums;      
    for (Uint k=0; k<byteblocks; k++) { // see permutation.tex (*)
      BlockType s0 = sum[k];
      for (Uint j=0; j<BytesPerUnit; j++, pS++) {
	BlockType L1;
	AND(L1, s0, first8bits);
	ADD32(*pS, *pS, L1); 
	SHR32(s0, s0, BitsPerByte);
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

#endif
}

      


SEXP allele_freqShuffle(SEXP GM) {
  SEXP Ans = R_NilValue; 
#if not defined ShuffleOK
  ERR("snpcoding=Shuffle needs at least SSSE3.");
#else
  assert_shuffle(GM);
  Uint 
    *info = GetInfo(GM),
   individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    *cgm = Align(GM, ALIGN_ALLELE);
  if (info[WHAT] != GENOMATRIX) ERR("not a coded Z matrix");
  PROTECT(Ans = allocVector(REALSXP, snps));
  allele_freq2((BlockType0 *) cgm, snps, individuals, NULL,
	       NULL); // Pd0 set as side effect
  real *Pd = (real*) algnrl(Pd0);
  double *ans = REAL(Ans);
  for (Uint i=0; i<snps; i++) ans[i] = 0.5 * (double) Pd[i];
  FREE(Pd0);
  UNPROTECT(1);
  if (PRESERVE) R_PreserveObject(Ans);
#endif
  return Ans;
}


