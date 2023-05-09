
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

#define BitsPerCode 2
#define MY_VARIANT 256

#define CORE_FACTOR 2

#include "Basic_miraculix.h"
#include"intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "Files.h"
#include "MX.h"
#include "Template.h"
#include "bitBase.h"
#include "2bit.h"

#if defined AVX2

ASSERT_SIMD(2bit256, avx2);

#include "haplogeno.h"
#include "intrinsicsIntern.h"
#include "2bitIntern.h"


typedef void (*scalarVectorGeno)(double *V, unit_t *Code, Long curLen, Long ld,
				 double *Ans);


sumGeno_start(2v256) 
  // printf("sumgeno2v256!\n");  
  return sumGenoIntern(code, snps, individuals, lda, cores, sums);
}   



void TwoBithaplo2geno256(unit_t *X, Long snps, Long individuals,
			 Long lda,  int cores, unit_t *Ans, Long ldAns) {
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  TwoBithaplo2genoIntern(X, snps, individuals, lda, cores, Ans, ldAns);
}


FileBinary(256)
SCALAR012(256)




void coding_2v256_4Byte(unit_t *M, Long ldM,
		       Long start_individual, Long end_individual, 
		       Long start_snp, Long end_snp,
		       int VARIABLE_IS_NOT_USED cores,
		       double VARIABLE_IS_NOT_USED *G,
		       unit_t * Ans,
		       Long VARIABLE_IS_NOT_USED individuals, Long ldAns) {
  BUG; // hat einen Fehler und ist auch nicth schneller
  
#define SizeFactor (BitsPerBlock / 128)
  if ((2 * start_snp) % CodesPerUnit != 0) // OK guaranted by calling function
    BUG;
     
  Long
     cur_snps = end_snp - start_snp,
    codes_atonce = CodesPerByte * UnitsPerBlock, // OK
    atonce128 = codes_atonce / SizeFactor,
    blocks = cur_snps / codes_atonce,
    rest128 = (cur_snps % codes_atonce) / atonce128,
    rest16 = cur_snps % atonce128;
  bool reste = rest128 + rest16 > 0;

  //  printf("cur_snps=%d %d %d %d 128once=%d\n", cur_snps, blocks, rest128, rest16,codes_atonce / SizeFactor);
  
  unit_t
    *Code= Ans + start_snp / CodesPerUnit + start_individual * (Long) ldAns;// OK
  Long deltaIndiv = end_individual - start_individual;

  BlockType two = SET32(2),
    wronga=ZERO(), 
    wrongb=ZERO(), 
    wrongc=ZERO(), 
    wrongd=ZERO();
  
#define Block4th uint64_t 
  const BlockType idx =  SET64(0L + (4L << 32L));
  for (Long ii=0; ii<deltaIndiv; ii++) {
    Block4th *code = (uint64_t *) (Code + ii * ldAns); // 256 / 4
    BlockType0 *pM = (BlockType0*) (M + ii * ldM);
    for (Long jj=0; jj<blocks; jj++) {
      const BlockType
	a0 = LOAD(pM++),
	b0 = LOAD(pM++),
	c0 = LOAD(pM++),
	d0 = LOAD(pM++);
      wronga = OR(wronga,GREATER_THAN(a0, two));
      wrongb = OR(wrongb,GREATER_THAN(b0, two));
      wrongc = OR(wrongc,GREATER_THAN(c0, two));
      wrongd = OR(wrongd,GREATER_THAN(d0, two));
      const BlockType
	a =  _mm256_permute2x128_si256(a0, b0, 0 + (2 << 4)),
	c =  _mm256_permute2x128_si256(c0, d0, 0 + (2 << 4)),
	b =  _mm256_permute2x128_si256(a0, b0, 1 + (3 << 4)),
	d =  _mm256_permute2x128_si256(c0, d0, 1 + (3 << 4));
      const BlockType
	b1 = SHL32(b, BitsPerByte),
	c1 = SHL32(c, 2 * BitsPerByte),
	d1 = SHL32(d, 3 * BitsPerByte);
      const BlockType
	ab = OR(a, b1),
	cd = OR(c1, d1);
      const BlockType
	e = OR(ab, cd),
	f = SHR64(e, 30),
	g = OR(e, f),	
	h = SHUFFLE32(g, 2),	
	i = SHL32(h, 4),
	j = OR(i, g);
      const BlockType k = _mm256_permutevar8x32_epi32(j, idx);
      code[jj] = (Block4th) _mm256_extract_epi64(k, 0);
    }

    wronga = OR(wronga, wrongb);
    wrongc = OR(wrongc, wrongd);
    UnionType wrongx; // OK
    wrongx.vi = OR(wronga, wrongc);
    bool wrong = wrongx.u64[0] | wrongx.u64[1] | wrongx.u64[2] | wrongx.u64[3];

    if (reste)
      wrong |= coding_2v128((unit_t*) pM, rest128, rest16,
			    (uint32_t*) (code + blocks));
    
    if (wrong) ERR0("2bit#256: values outside 0,1,2.");
    
  } // for delteIndiv
}


/* Plink: FOR AVX512 masked-blend?!
  00 -> 00
  01 : missing
  10 -> 01
  11 -> 10

#define x0000 0
#define x0001 1
#define x0010 2
#define x0011 3
#define x0100 4
#define x0101 5
#define x0110 6
#define x0111 7
#define x1000 8
#define x1001 9
#define x1010 A
#define x1011 B
#define x1100 C
#define x1101 D
#define x1110 E
#define b1111 F
#define BB1(A,B,C) A##B##C
#define BB0(A,B,C)  BB1(A,B,C) 
#define B(X,Y) BB0(0x,x##X,x##Y)



  LOWER 4 BIT
  01 -> 00, dcba -> dbca, 
  ------------
  0000 -> 00000000
  0001 -> 00000000
  0010 -> 00010000
  0011 -> 00010001
  0100 -> 00000000
  0101 -> 00000000
  0110 -> 00010000
  0111 -> 00010001
  1000 -> 00100000
  1001 -> 00100000
  1010 -> 00110000
  1011 -> 00110001
  1100 -> 00100010
  1101 -> 00100010
  1110 -> 00110010
  1111 -> 00110011

  UPPER 4 BIT
  LOWER 4 BIT SHL 2

##########################################

  2 bit

  LOWER 4 BIT
  10 -> 11; dcba -> dbca,
  ------------
  0000 -> 00000000
  0001 -> 00000001
  0010 -> 00010001
  0011 -> 00010001 // does not exist
  0100 -> 00000010
  0101 -> 00000011
  0110 -> 00010011
  0111 -> 00010011 // does not exist
  1000 -> 00100010
  1001 -> 00100011
  1010 -> 00110011
  1011 -> 00110011 // does not exist
  1100 -> 00100010 // does not exist
  1101 -> 00100011 // does not exist
  1110 -> 00110011 // does not exist
  1111 -> 00110011 // does not exist

  LOWER 4 BIT SHL 2

##########################################


const BlockType plinkTable =
  SETREV8(
	  B(0000,0000),
	  B(0000,0000),
	  B(0001,0000),
	  B(0001,0001),
	  B(0000,0000),
	  B(0000,0000),
	  B(0001,0000),
	  B(0001,0001),
	  B(0010,0000),
	  B(0010,0000),
	  B(0011,0000),
	  B(0011,0001),
	  B(0010,0010),
	  B(0010,0010),
	  B(0011,0010),
	  B(0011,0011)),
  TwoBitTable =
  SETREV8(	    
	  B(0000,0000),
	  B(0000,0001),
	  B(0001,0001),
	  B(0001,0001), //
	  B(0000,0010),
	  B(0000,0011),
	  B(0001,0011),
	  B(0001,0011), // does not exist
	  B(0010,0010),
	  B(0010,0011),
	  B(0011,0011),
	  B(0011,0011), // does not exist
	  B(0010,0010), // does not exist
	  B(0010,0011), // does not exist
	  B(0011,0011), // does not exist
	  B(0011,0011)); // does not exist
*/



#define code2blockN(N,code) BlockType *code##N = (BlockType0 *) (Code + N * ld)
#define ZeroN(N,sum) Doubles sum##N = ZERODOUBLE()
#define LoadN(N, a, code) const BlockType a##N = LOADU(code##N++)
#define andConstN(N, b, a, Const) const BlockType b##N = AND(a##N, Const)
#define shrConstN(N, c, b, Const) const BlockType c##N = SHR64(b##N, Const)
#define AndN(N, d, a, c) const BlockType d##N = AND(a##N, c##N)
#define OrN(N, d, a, c)	BlockType d##N = OR(a##N, c##N)
#define blend0N(N,g,f,e) const Doubles g##N = BLENDvDOUBLE(zero, f, (Doubles)e##N)
#define AddUpN(N,sum, g) sum##N = ADDDOUBLE(g##N, sum##N)
#define shl1N(N,c) c##N = SHL64(c##N, 1)
#define addAnsN(N, Ans, s, sum) \
  s = (double*) &sum##N;  Ans[N] += s[0] + s[1] + s[2] + s[3]


#define code2block(code) code2blockN(0,code)
#define Zero(sum) ZeroN(0, sum)
#define Load(a, code) LoadN(0, a, code)
#define andConst(b,a,Const) andConstN(0, b, a, Const)
#define shrConst(c, b, Const) shrConstN(0, c, b, Const)
#define And(d, a, c) AndN(0, d, a, c)
#define Or(d, a, c) OrN(0, d, a, c)
#define blend0(g, f, e)	blend0N(0,g,f,e)
#define AddUp(sum, g) AddUpN(0, sum, g)
#define shl1(c)	shl1N(0,c)
#define addAns(Ans, s, sum) addAnsN(0, Ans, s, sum)


#define TRAFO_PLINK\
  And(d, a, c);	   \
  Or(e, d, b)

#define TRAFO_TWOBIT\
  Or(e, a, c)


#define scalar1xX(NAME,TRAFO)						\
  void scalar1x##NAME(double *V, unit_t *Code, Long curLen, Long ld,	\
		      double *Ans) {					\
  Doubles *v = (Doubles*) V;						\
  Doubles zero = ZERODOUBLE();						\
  code2block(code);							\
  Zero(sum);								\
  for (Long i=0; i<curLen; i+=CodesPerBlock) {				\
    Load(a, code);							\
    andConst(b, a, E1N1);						\
    shrConst(c, b, 1);							\
    TRAFO;								\
    for (Long j=0; j<(CodesPerBlock / doubles); j++) {			\
      const Doubles f = LOADuDOUBLE((double*) v++);			\
      blend0(g, f, e);							\
      AddUp(sum, g);							\
      shl1(e);								\
      blend0(h, f, e);							\
      AddUp(sum, h);							\
      shl1(e);								\
    }									\
  }									\
  double *s;								\
  addAns(Ans, s, sum);							\
  }


#define SCALAR1xX(NAME)							\
  scalar1xX(NAME,TRAFO_TWOBIT)						\
  scalar1xX(NAME##plink,TRAFO_PLINK)

SCALAR1xX(1)




#define BroadcastN(N, code, V) const Doubles code##N = BROADCAST(V++)
#define blendMultiN(N, g, f, e) \
  const Doubles g##N = BLENDvDOUBLE(zero, f##N, (Doubles) e##N)
#define AddN(N, k, g, h) const Doubles k##N = ADDDOUBLE(g##N, h##N)
#define LoadDoubleN(N,m,ans) const Doubles m##N = LOADuDOUBLE((double*)(ans + N))
#define miniAddN(N, n, m, l) const Doubles n##N = ADDDOUBLE(m##N, l##N)
#define StoreDoubleN(N, ans, n) STOREuDOUBLE((double*) (ans++), n##N)


#define Broadcast(code, V) BroadcastN(0, code, V)
#define blendMulti(g, f, e) blendMultiN(0, g, f, e)
#define Add(k, g, h) AddN(0, k, g, h)
#define LoadDouble(m, ans) LoadDoubleN(0, m, ans); LoadDoubleN(1, m, ans); \
  LoadDoubleN(2, m, ans); LoadDoubleN(3, m, ans);
#define miniAdd(n, m, l) miniAddN(0, n, m, l); miniAddN(1, n, m, l);	\
  miniAddN(2, n, m, l); miniAddN(3, n, m, l)
#define StoreDouble(ans, n) StoreDoubleN(0, ans, n); StoreDoubleN(1, ans, n); \
  StoreDoubleN(2, ans, n); StoreDoubleN(3, ans, n)

#define addquer_shift(l, k, e)				\
  const Doubles l = k##0;				\
  shl1(e)					       


#define GVblock(N)					\
  blendMulti(g##N, multiV, e);				\
  shl1(e);						\
  blendMulti(h##N, multiV, e);				\
  Add(k##N, g##N, h##N);				\
  addquer_shift(l##N, k##N, e)


#define ScalarGV(NAME,TRAFO, BLOCKS, ATONCE)				\
  void scalarGV##NAME(double *V, unit_t *Code,Long curLen,Long ld, double*Ans){ \
    code2block(code);							\
    Doubles *ans = (Doubles*) Ans;					\
    Doubles zero = ZERODOUBLE();						\
    Broadcast(multiV, V);						\
    /*double *x = NULL;	*/						\
    for (Long i=0; i<curLen; i+=CodesPerBlock) {				\
      Load(a, code);							\
      andConst(b, a, E1N1);						\
      shrConst(c, b, 1);						\
      TRAFO;								\
      for (Long j=0; j<(CodesPerBlock / doubles / ATONCE); j++) {	\
	BLOCKS;								\
	LoadDouble(m, ans);						\
	miniAdd(n, m, l);						\
	StoreDouble(ans, n);						\
      }									\
    }									\
  }

#define SCALAR_GV(NAME, BLOCKS, ATONCE) \
  ScalarGV(NAME,TRAFO_TWOBIT, BLOCKS, ATONCE)	\
  ScalarGV(NAME##plink, TRAFO_PLINK, BLOCKS, ATONCE)	

SCALAR_GV(1, GVblock(0); GVblock(1); GVblock(2); GVblock(3), 4)

#if defined code2block
#undef code2block
#undef Zero
#undef Load
#undef andConst
#undef shrConst
#undef And
#undef Or
#undef blend0
#undef AddUp
#undef shl1
#undef addAns
#undef Broadcast
#undef addquer_shift
#undef Add
#undef LoadDouble
#undef miniAdd
#undef StoreDouble
#undef blendMulti
#endif


#define code2block(code) code2blockN(0,code); code2blockN(1,code)
#define Zero(sum) ZeroN(0, sum); ZeroN(1, sum)
#define Load(a, code) LoadN(0, a, code); LoadN(1, a, code)
#define andConst(b,a,Const) andConstN(0, b, a, Const); andConstN(1, b, a, Const)
#define shrConst(c, b, Const) shrConstN(0,c, b, Const); shrConstN(1,c, b, Const)
#define And(d, a, c) AndN(0, d, a, c); AndN(1, d, a, c)
#define Or(d, a, c) OrN(0, d, a, c); OrN(1, d, a, c)
#define blend0(g, f, e)	blend0N(0,g,f,e); blend0N(1,g,f,e)
#define AddUp(sum, g) AddUpN(0, sum, g); AddUpN(1, sum, g) 
#define shl1(c)	shl1N(0,c); shl1N(1,c)	
#define addAns(Ans, s, sum) addAnsN(0, Ans, s, sum); addAnsN(1, Ans, s, sum)


SCALAR1xX(2)


#define Broadcast(code, V) BroadcastN(0, code, V); BroadcastN(1, code, V) 
#define blendMulti(g, f, e) blendMultiN(0, g, f, e); blendMultiN(1, g, f, e)
#define Add(k, g, h) AddN(0, k, g, h); AddN(1, k, g, h)
#define LoadDouble(m, ans) LoadDoubleN(0, m, ans); LoadDoubleN(1, m, ans); \
  LoadDoubleN(2, m, ans); LoadDoubleN(3, m, ans);
#define miniAdd(n, m, l) miniAddN(0, n, m, l); miniAddN(1, n, m, l);	\
  miniAddN(2, n, m, l); miniAddN(3, n, m, l)
#define StoreDouble(ans, n) StoreDoubleN(0, ans, n); StoreDoubleN(1, ans, n); \
  StoreDoubleN(2, ans, n); StoreDoubleN(3, ans, n)

#define addquer_shift(l, k, e)				\
  const Doubles l = ADDDOUBLE(k##0, k##1);		\
  shl1(e)

SCALAR_GV(2, GVblock(0); GVblock(1); GVblock(2); GVblock(3), 4)


#if defined code2block
#undef code2block
#undef Zero
#undef Load
#undef andConst
#undef shrConst
#undef And
#undef Or
#undef blend0
#undef AddUp
#undef shl1
#undef addAns
#undef Broadcast
#undef addquer_shift
#undef Add
#undef LoadDouble
#undef miniAdd
#undef StoreDouble
#undef blendMulti
#endif


#define code2block(code)			\
  code2blockN(0,code); code2blockN(1,code);	\
  code2blockN(2,code); code2blockN(3,code)
#define Zero(sum)				\
  ZeroN(0, sum); ZeroN(1, sum);			\
  ZeroN(2, sum); ZeroN(3, sum)
#define Load(a, code)			\
  LoadN(0, a, code); LoadN(1, a, code); \
  LoadN(2, a, code); LoadN(3, a, code)
#define andConst(b,a,Const)				\
  andConstN(0, b, a, Const); andConstN(1, b, a, Const);	\
  andConstN(2, b, a, Const); andConstN(3, b, a, Const)
#define shrConst(c, b, Const)				\
  shrConstN(0,c, b, Const); shrConstN(1,c, b, Const);	\
  shrConstN(2,c, b, Const); shrConstN(3,c, b, Const)
#define And(d, a, c)				\
  AndN(0, d, a, c); AndN(1, d, a, c);		\
  AndN(2, d, a, c); AndN(3, d, a, c)
#define Or(d, a, c)				\
  OrN(0, d, a, c); OrN(1, d, a, c);		\
  OrN(2, d, a, c); OrN(3, d, a, c)
#define blend0(g, f, e)				\
  blend0N(0,g,f,e); blend0N(1,g,f,e);		\
  blend0N(2,g,f,e); blend0N(3,g,f,e)
#define AddUp(sum, g)				\
  AddUpN(0, sum, g); AddUpN(1, sum, g);		\
  AddUpN(2, sum, g); AddUpN(3, sum, g) 
#define shl1(c)					\
  shl1N(0,c); shl1N(1,c);			\
  shl1N(2,c); shl1N(3,c)	
#define addAns(Ans, s, sum)				\
  addAnsN(0, Ans, s, sum); addAnsN(1, Ans, s, sum);	\
  addAnsN(2, Ans, s, sum); addAnsN(3, Ans, s, sum)

SCALAR1xX(4)







vectorGeno_header(double, double, double, 2v256) {
  int cores = GreaterZero(opt->cores);
  Long 
    m = tuning->logSnpGroupSize,
    n = tuning->indivGroupSize;
  Long miniRows = tuning->miniRows <= 0 ? 1 : tuning->miniRows;  
  Long miniCols = tuning->miniCols <= 0 ? 4 : tuning->miniCols;
  Long tilerepetV = 1;


  miniCols = 1; // 4.812 (  4.838) sec for calculating V * RelMatrix (V256)
  miniCols = 2; // 3.543 (  3.569) sec for calculating V * RelMatrix (V256)
  miniCols = 4; // 3.044 (  3.070) sec for calculating V * RelMatrix (V256)
  
  
  Long
    //    ldAns = cols,
    bytes = sizeof(*Ans) * cols,
    tileColsCode = n;                                 // of result
  Long tileRows = m <= 0 || (1L << m) > rows ? rows : (1L << m);// for rowPos
    Long codeIncr = miniCols * lda;
  if (tileColsCode <= 0) {
    tileColsCode = ( ( ((repetV - 1L) / tilerepetV + 1L) * cols - 1L )
		     / (cores * CORE_FACTOR) + 1L );
  }
  tileColsCode = tileColsCode > cols ? cols : tileColsCode;
  // n <= 0 || n > cols ? cols : n; // of result  
  tileRows = (tileRows / CodesPerBlock) * CodesPerBlock;
  if (tileRows == 0) {
    PRINTF("rows=%ld CpB=%u\n", rows, CodesPerBlock);
    ERR0("Choice of 'logSnpGroupSize' too small.\n");
  }


  //  printf("cC %d %d   %d %d\n", cols , tileColsCode, rows,tileRows);
  if (CodesPerBlock % miniCols != 0) BUG;
  
  if  (CodesPerBlock == 0) BUG;
  //  printf("%d %d %d\n", rows, CodesPerBlock, tileRows);
  if (rows % tileRows != 0) {
    PRINTF("vG: rows=%ld not a multiple of tileRows=%ld\n", rows, tileRows);
    BUG;
    }
  if (tileRows % CodesPerBlock != 0) {
    PRINTF("vG: tileRows=%ld not a multiple of CpB=%d\n", tileRows, CodesPerBlock);
    BUG;
  }
  if (cols % tileColsCode != 0) {
    PRINTF("vG: cols=%ld not a multiple of tileCols=%ld\n", cols, tileColsCode);
    BUG;
  }
  if (tileColsCode % miniCols != 0) {
    PRINTF("vG: tileCols=%ld not a multiple of miniCols=%ld\n", tileColsCode, miniCols);
    BUG;
  }
  
  scalarVectorGeno scalarVG;
  if (coding == TwoBitGeno || coding == TwoBitGenoTransposed) {
    switch(miniCols) {
    case 1 : scalarVG = scalar1x1; break;
    case 2 : scalarVG = scalar1x2; break;
    case 4 : scalarVG = scalar1x4; break;
    default : BUG;
    }
  } else if (coding == Plink || coding == PlinkTransposed) {
   switch(miniCols) {
    case 1 : scalarVG = scalar1x1plink; break;
    case 2 : scalarVG = scalar1x2plink; break;
    case 4 : scalarVG = scalar1x4plink; break;
    default : BUG;
    }
  } else BUG;


  
  MEMSET(Ans, 0,  bytes);
#define  CodesPerDouble (SizeOfDouble * BitsPerByte / BitsPerCode)
  assert(repetV == 1); // check v, ob repetV fehlt!
  double *v = (double*) MALLOC(sizeof(double) * rows);
  Long endj = rows / CodesPerBlock;
  if (rows != endj * CodesPerBlock) BUG;
  for (Long j=0; j<endj; j++) {
    double
      *out = v + j * CodesPerBlock,
      *in = V + j * CodesPerBlock;
    for (Long i=0; i<CodesPerDouble; i++) {
      out[ (CodesPerBlock - doubles + 0) - i * doubles] = in[i + 0];
      out[ (CodesPerBlock - doubles + 1) - i * doubles] = in[i + 32];
      out[ (CodesPerBlock - doubles + 2) - i * doubles] = in[i + 64];
      out[ (CodesPerBlock - doubles + 3) - i * doubles] = in[i + 96];
    }
  } 
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) collapse(2)
#endif
  for (Long tR=0; tR<repetV; tR += tilerepetV) {
    for (Long tC=0; tC<cols; tC += tileColsCode) {
      Long Rend = tR + tilerepetV < repetV ? tR + tilerepetV : repetV,
	i_end = Rend - miniRows + 1;      
      Long Cend = tC + tileColsCode < cols ? tC + tileColsCode : cols,
	j_end = Cend - miniCols + 1;	  
      for (Long rowPos = 0; rowPos < rows; rowPos += tileRows) {
	Long curLen = tileRows < rows - rowPos? tileRows : rows - rowPos;
	Long i = tR;
	for ( ; i<i_end; i+= miniRows) {
	  Long j = tC;
	  //	  printf("i=%ld++%ld<%ld %d %ld %f\n", i, miniRows, i_end, ldV, rowPos, v[0]);
	  double *VV0 = v + i * ldV + rowPos,
	    *A = Ans + i * ldAns;
	  unit_t *Code = code + j * lda + rowPos / CodesPerUnit; // OK
	  for ( ; j<j_end; j+=miniCols, Code+=codeIncr) {
	    scalarVG(VV0, Code, curLen, lda, A + j);
	  }
	}	
      }
    }
  }
  FREE(v);
}



#define Broadcast(code, V)				\
  const Doubles code##0 = BROADCAST(V++);	\
  const Doubles code##1 = BROADCAST(V++);	\
  const Doubles code##2 = BROADCAST(V++);	\
  const Doubles code##3 = BROADCAST(V++)

#define blendMulti(g, f, e)					\
  const Doubles g##0 = BLENDvDOUBLE(zero, f##0, (Doubles) e##0);	\
  const Doubles g##1 = BLENDvDOUBLE(zero, f##1, (Doubles) e##1);	\
  const Doubles g##2 = BLENDvDOUBLE(zero, f##2, (Doubles) e##2);	\
  const Doubles g##3 = BLENDvDOUBLE(zero, f##3, (Doubles) e##3) 


#define Add(k, g, h)						\
   const Doubles k##0 = ADDDOUBLE(g##0, h##0);		\
   const Doubles k##1 = ADDDOUBLE(g##1, h##1);		\
   const Doubles k##2 = ADDDOUBLE(g##2, h##2);		\
   const Doubles k##3 = ADDDOUBLE(g##3, h##3)

#define addquer_shift(l, k, e)				\
  const Doubles l##0 = ADDDOUBLE(k##0, k##1);		\
  const Doubles l##1 = ADDDOUBLE(k##2, k##3);		\
  shl1(e);						\
  const Doubles l  = ADDDOUBLE(l##0, l##1)		\

#define LoadDouble(m, ans)				\
  const Doubles m##0 = LOADuDOUBLE((double*) (ans + 0));	\
  const Doubles m##1 = LOADuDOUBLE((double*) (ans + 1))
 
#define miniAdd(n, m, l)			\
  const Doubles n##0 = ADDDOUBLE(m##0, l##0);	\
  const Doubles n##1 = ADDDOUBLE(m##1, l##1)

#define StoreDouble(ans, n)		\
  STOREuDOUBLE((double*) (ans++), n##0);	\
  STOREuDOUBLE((double*) (ans++), n##1)


SCALAR_GV(4, GVblock(0); GVblock(1), 2)

#if defined code2block
#undef code2block
#undef Zero
#undef Load
#undef andConst
#undef shrConst
#undef And
#undef Or
#undef blend0
#undef AddUp
#undef shl1
#undef addAns
#undef Broadcast
#undef addquer_shift
#undef Add
#undef LoadDouble
#undef miniAdd
#undef StoreDouble
#undef blendMulti
#undef atOnce
#endif


#define code2block(code) \
  code2blockN(0,code); code2blockN(1,code); code2blockN(2,code)
#define Zero(sum) ZeroN(0, sum); ZeroN(1, sum); ZeroN(2, sum)
#define Load(a, code) LoadN(0, a, code); LoadN(1, a, code); LoadN(2, a, code); 
#define andConst(b,a,Const) \
  andConstN(0, b, a, Const); andConstN(1, b, a, Const); andConstN(2, b, a,Const)
#define shrConst(c, b, Const) \
  shrConstN(0,c, b, Const); shrConstN(1,c, b, Const); shrConstN(2,c, b, Const)
#define And(d, a, c) AndN(0, d, a, c); AndN(1, d, a, c); AndN(2, d, a, c)
#define Or(d, a, c) OrN(0, d, a, c); OrN(1, d, a, c); OrN(2, d, a, c)
#define blend0(g, f, e) blend0N(0,g,f,e); blend0N(1,g,f,e); blend0N(2,g,f,e)
#define AddUp(sum, g) AddUpN(0, sum, g); AddUpN(1, sum, g); AddUpN(2, sum, g)
#define shl1(c) shl1N(0,c); shl1N(1,c); shl1N(2,c)
#define addAns(Ans, s, sum)\
  addAnsN(0, Ans, s, sum); addAnsN(1, Ans, s, sum); addAnsN(2, Ans, s, sum)


#define Broadcast(code, V) \
  BroadcastN(0, code, V); BroadcastN(1, code, V); BroadcastN(2, code, V) 
#define blendMulti(g, f, e) \
  blendMultiN(0, g, f, e); blendMultiN(1, g, f, e); blendMultiN(2, g, f, e)
#define Add(k, g, h) AddN(0, k, g, h); AddN(1, k, g, h); AddN(2, k, g, h)
#define LoadDouble(m, ans) LoadDoubleN(0, m, ans); LoadDoubleN(1, m, ans)
#define miniAdd(n, m, l) miniAddN(0, n, m, l); miniAddN(1, n, m, l);
#define StoreDouble(ans, n) StoreDoubleN(0, ans, n); StoreDoubleN(1, ans, n)

#define addquer_shift(l, k, e)				\
  const Doubles l##0 = ADDDOUBLE(k##0, k##1);		\
  shl1(e);						\
  const Doubles l  = ADDDOUBLE(l##0, k##2)		\


SCALAR_GV(3, GVblock(0); GVblock(1), 3)


//vectorGeno_header(LongDouble, LongDouble, LongDouble, 2v256) {
//  BUG;
  //  printf("xxx= %f %d\n", (double) (V[0] + Ans[0]), code[rows + cols + lda + repetV + ldV + ldAns]);
//}


#define  CodesPerDouble (SizeOfDouble * BitsPerByte / BitsPerCode)
genoVector_header(double, double, double, 2v256) {
  int cores = GreaterZero(opt->cores);
  Long
    tilerepetV = 1,
    tileColsCodebundle = cols; // should be improved
				    //    ldV = cols,
  Long 
    m = tuning->logSnpGroupSize,
    n = tuning->indivGroupSize;
  Long miniRows = tuning->miniRows <= 0 ? 1 : tuning->miniRows;  
  Long miniCols = tuning->miniCols <= 0 ? 4 : tuning->miniCols;

  // Wageningen/run rows=51200 cols=102400 cores=2 coding=2 mode=1 repetV=1 variant=32 SNPmode=5 cmp=0
  // miniCols = 1; // 3.057 (  3.083) sec for calculating V * RelMatrix (V256)
  //  miniCols = 2; // 2.816 (  2.842) sec for calculating V * RelMatrix (V256)
  //  miniCols = 4; // 2.75 (  2.820) sec for calculating V * RelMatrix (V256)
  
  Long bytes = sizeof(*Ans) * rows,
     tileColsCode = n <= 0 || n > cols ? cols : n, // of result
    tileRows = m > 0 ? 1<<m : rows / (cores * CORE_FACTOR);//for rowPos
  Long codeIncr = miniCols * lda;
  tileRows = ((tileRows > rows ? rows : tileRows) / CodesPerBlock) * CodesPerBlock;
  if (tileRows == 0) {
    PRINTF("gV_double256: rows=%ld CpB=%u\n", rows, CodesPerBlock);
    ERR0("Choice of 'logSnpGroupSize' too small.\n");
  }
  scalarVectorGeno scalarGV = NULL;
  //  printf("miniCols = %d\n", miniCols);
  // BUG;
 if (coding == TwoBitGeno) {
   switch(miniCols) {
   case 1 : scalarGV = scalarGV1; break;
   case 2 : scalarGV = scalarGV2; break;
   case 3 : scalarGV = scalarGV3; break;
   case 4 : scalarGV = scalarGV4; break;
   default : BUG;
   }
 } else if (coding == Plink) {
   switch(miniCols) {
   case 1 : scalarGV = scalarGV1plink; break;
   case 2 : scalarGV = scalarGV2plink; break;
   case 3 : scalarGV = scalarGV3plink; break;
   case 4 : scalarGV = scalarGV4plink; break;
   default : BUG;
   }
 }
 assert(scalarGV != NULL);
 
   if (rows % CodesPerBlock != 0) {
    PRINTF("rows=%ld not a multiple of CpB=%u\n", rows, CodesPerBlock);
    BUG;
  }
  if (tileRows % CodesPerBlock != 0) {
    PRINTF("tileRows=%ld not a multiple of CpB=%d\n", tileRows, CodesPerBlock);
    BUG;
  }
  if (cols % tileColsCode != 0) {
    PRINTF("cols=%ld not a multiple of tileCols=%ld\n", cols, tileColsCode);
    BUG;
  }
  if (tileColsCode % miniCols != 0) {
    PRINTF("tileCols=%ld not a multiple of miniCols=%ld\n", tileColsCode, miniCols);
    BUG;
  }
  
  MEMSET(Ans, 0,  bytes);

#ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(cores) schedule(static) collapse(3) 
#endif
  for (Long tR=0; tR<repetV; tR += tilerepetV) {
    for (Long rowPos = 0; rowPos < rows; rowPos += tileRows) {
      for (Long tCbundle=0; tCbundle<cols; tCbundle += tileColsCodebundle) {
	// todo: speicher aufmachen!
	Long tCbundleEnd = tCbundle + tileColsCodebundle < cols ?
	  tCbundle + tileColsCodebundle : cols;	
    
	for (Long tC=tCbundle; tC<tCbundleEnd; tC += tileColsCode) {
	  Long Rend = tR + tilerepetV < repetV ? tR + tilerepetV : repetV,
	    i_end = Rend - miniRows + 1, 
	    curLen =tileRows < rows - rowPos ? tileRows : rows - rowPos;
	  Long Cend = tC + tileColsCode < cols ? tC + tileColsCode : cols,
	    j_end = Cend - miniCols + 1;
	  Long i = tR;
	  for ( ; i<i_end; i+= miniRows) { //auch unnoetig, falls immer 1 Vector
	    //	  printf("i=%d < %d\n", i), i_end;
	    Long j = tC;
	    //	  printf("i=%ld %d %ld++%ld<%d %f\n", i, ldV, rowPos, tileRows, rows, v[0]);
	    double *vv = V + i * ldV,
	      *a = Ans + i * ldAns + rowPos;
	    unit_t *Code = code + j * lda + rowPos / CodesPerUnit; // OK
	    //  printf("j=%d %d miniC=%d, rowPos=%d %d %d xurLen=%d\n", j, ld, miniCols, rowPos, rows, rows, curLen);	  
	    for ( ; j<j_end; j+=miniCols, Code+=codeIncr, vv+=miniCols) {
	      // if (curLen < 0 || curLen > 10000) {printf("curlen=%d\n", curLen); BUG;}
	      scalarGV(vv, Code, curLen, lda, a);
	    }
	    //	  printf("j done\n");
	  }	
	  //	printf("i done\n");
	}
	//      printf("tC done\n");
      } 
      //      printf("tC bundle done\n");
    }
    //	printf("rowPos done\n");
  }

  //  printf("loop done\n");

  Long endj = rows / CodesPerBlock;
  if (rows != endj * CodesPerBlock) BUG;
  double in[CodesPerBlock];

  for (Long j=0; j<endj; j++) {
    for (Long r=0; r<repetV; r++) {
      double *out = Ans + j * CodesPerBlock + r * ldAns;//??? r * cols? 
      MEMCOPY(in, out, CodesPerBlock * SizeOfDouble);
      for (Long i=0; i<CodesPerDouble; i++) {
 	for (Long k=0; k<doubles; k++) {
	  //     	printf("%d rows=%d %d < -  %d CpD=%d\n", i, rows, i + k * 32, (CodesPerBlock - doubles + k) - i * doubles,  CodesPerDouble);
	  //	printf("%d rows=%d %d < - %d\n", i, rows, i + k * 32, (CodesPerBlock - doubles + k) - i * doubles);
	  out[i + k * 32] = in[ (CodesPerBlock - doubles + k) - i * doubles];
	}
      }
    }
  }  
}
 

void allele_sum_2v256(unit_t *Code, Long snps, Long individuals,
		      Long lda, Long LDAbitalign,
		      basic_options VARIABLE_IS_NOT_USED *opt,
		      Ulong *Ans) {
  allele_sum_I_2bit(Code, snps, individuals, lda, LDAbitalign, NULL, NULL, Ans);
}


 
#else // !defined AVX
#include "avx_miss.h"


genoVector_header(double, double, double, 2v256) Sv
genoVector_header(LongDouble, LongDouble, LongDouble, 2v256) Sv
vectorGeno_header(double, double, double, 2v256) Sv
vectorGeno_header(LongDouble, LongDouble, LongDouble, 2v256) Sv


DeclareVersion2(256,Su,Sv,Sm)
SIMD_MISS(2bit256, avx2);

#endif  // defined AVX2
// trafo2Geno1Geno128
