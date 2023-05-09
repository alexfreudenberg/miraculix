
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


#ifndef miraculix_2bitIntern_H
#define miraculix_2bitIntern_H 1


#if defined SHUFFLE8

const BlockType
  to_scalar =  SETREV8( 0, 1, 4, 2,  1, 2, 5, 3,  4, 5, 8, 6,  2, 3, 6, 4);

static Ulong scalar(BlockType0 *x, BlockType0 *y, Long blocks,
		   Long VARIABLE_IS_NOT_USED Delta) {
  const BlockType zero = ZERO();
  Ulong ans = 0;

  for (char k=0; k<=1; k++) {
    Long bundles;
    int minibundle;
    if (k) {
      minibundle = (int) (blocks % loop_2bit);
      x += bundles * loop_2bit;
      y += bundles * loop_2bit;
      bundles = minibundle > 0;
    } else {
      bundles = blocks  / loop_2bit;
      minibundle = loop_2bit;
    }
    
    for (Long i=0; i<bundles; i++) { 
      //printf("k=%d i=%d bundles=%d\n", k, i, bundles);
     BlockType
	sum = ZERO(),
	*xx = x + i * minibundle,
	*yy = y + i * minibundle;     
      for (int j=0; j<minibundle; j++) {
	//	printf("k=%d i=%d j=%d %d\n", k, i, j, minibundle);
	const BlockType a0 = LOAD(xx++);
	const BlockType b0 = LOAD(yy++);
	const BlockType z0 = AND(a0, b0);
	const BlockType z1 = XOR(a0, b0);
	
	const BlockType z2 = AND(z1, E1N1);
	const BlockType z3 = SHR32(z2, 1L);
	const BlockType z4 = AND(z3, z1);
	
	const BlockType z5 = SHL32(z4, 1); 
	const BlockType z6 = OR(z4, z5);
	const BlockType z7 = OR(z6, z0); 
	
	const BlockType z8 = AND(z7, O1F1); 
	const BlockType z9 = SHR32(z7, 4);
	const BlockType z10 = AND(z9, O1F1); 
	
	const BlockType z11 = SHUFFLE8(to_scalar, z8); 
	const BlockType z12 = SHUFFLE8(to_scalar, z10); 
	sum = ADD8(sum, z11);
	sum = ADD8(sum, z12); 
      }
      const BlockType vsum = SAD8(sum, zero);
      ans += HORIZ_ADD64(vsum);
    }
  }
  return ans;
}


static void scalar_2x2(BlockType0 *x, BlockType0 *y,
		       BlockType0 *u, BlockType0 *v,
		       Long blocks,
		       Ulong *ansA, Ulong *ansB, Ulong *ansC, Ulong *ansD
		       ) {
  const BlockType zero = ZERO();
  Ulong a=0, b=0, c=0, d=0;
  //#define y15 (255L / 16L)  // 16L = 4 * 2^2 (a byte full of 1's)
  //#define y15 (255L / 16L)  // warum ist dies schneller ?????

  for (char k=0; k<=1; k++) {
    Long  bundles = NA_LONG;
    int minibundle;
    if (k) {
      minibundle = (int) (blocks % loop_2bit);
      x += bundles * loop_2bit;
      y += bundles * loop_2bit;
      u += bundles * loop_2bit;
      v += bundles * loop_2bit;
      bundles = minibundle > 0;
    } else {
      bundles = blocks  / loop_2bit;
      minibundle = loop_2bit;
    }
    for (Long i=0; i<bundles; i++) { 
      BlockType
	asum = ZERO(),
	bsum = ZERO(),
	csum = ZERO(),
	dsum = ZERO(),    
	*xx = x + i * minibundle,
	*yy = y + i * minibundle,
	*uu = u + i * minibundle,
	*vv = v + i * minibundle;
      for (int j=0; j<minibundle; j++) {
	const BlockType aa = LOAD(xx++);
	const BlockType bb = LOAD(yy++);
	const BlockType cc = LOAD(uu++);
	const BlockType dd = LOAD(vv++);
	
	const BlockType a0 = AND(aa, bb); // x,y
	const BlockType b0 = AND(cc, dd); // u,v
	const BlockType a1 = XOR(aa, bb);
	const BlockType b1 = XOR(cc, dd);
	
	const BlockType c0 = AND(aa, dd); // x, v
	const BlockType d0 = AND(cc, bb); // u, y
	const BlockType c1 = XOR(aa, dd);
	const BlockType d1 = XOR(cc, bb);
	
	AND4_C(2, 1, E1N1); //
	SHR32_4(3, 2, 1L); //
	AND4(4, 3, 1);
	
	SHL32_4(5, 4, 1L); //
	OR4(6, 4, 5); //
	OR4(7, 6, 0); //
	
	AND4_C(8, 7, O1F1); 
	SHR32_4(9, 7, 4);//
	AND4_C(10, 9, O1F1); //
	
	SHUFFLE8_4(11, to_scalar, 8); //
	SHUFFLE8_4(12, to_scalar, 10);//
	
	ADDUP8_4(sum, 11);      
	ADDUP8_4(sum, 12); 
      }
      SAD8_4(vsum, sum, zero);
      HORIZ_ADD64_4(vsum);
     }
  }
  *ansA = a;
  *ansB = b;
  *ansC = c;
  *ansD = d;
}

#define SCALAR012(NR)						\
  Ulong scalar_2v##NR(unit_t *x, unit_t *y, Long blocks, Long Delta) {	\
    return scalar((BlockType0*) x, (BlockType0*) y, blocks, Delta);	\
  }									\
  void multi_2v##NR(unit_t *x, unit_t *y,				\
		    Long VARIABLE_IS_NOT_USED individuals,		\
		    Long blocks, Long VARIABLE_IS_NOT_USED lda,	\
		   double *ans) {			\
    *ans += (double) scalar((BlockType0*) x, (BlockType0*) y, blocks, 0); \
  }									\
  void multi2x2_2v##NR(unit_t *x, unit_t *y, Long individuals, \
		       Long blocks, Long lda,		   \
		       double *ans) {					\
    Ulong A=0, B=0, C=0, D=0;						\
    scalar_2x2((BlockType0 *) x, (BlockType0 *) y,			\
	       (BlockType0 *) (x + lda), (BlockType0 *) (y + lda),	\
	       blocks, &A, &B, &C, &D);					\
    *ans += (double) A;							\
    ans[individuals + 1] += (double) B;					\
    ans[individuals] += (double) C;					\
    ans[1] += (double) D;						\
}



static sumGeno_start(Intern) 
  Ulong *cS = sums;
  if (cS == NULL) cS = (Ulong*) MALLOC(individuals * sizeof(Ulong));

#define maxint 2147483647  
#define x31 (255L / (2 * CodesPerByte)) // 2 : max{0,1,2}

  const Long
    blocks = Blocks(snps),
    rest = blocks % x31,
    kend = rest > 0,
    blocks31 = blocks - rest;

   const BlockType
    zero = ZERO(),
    N6E2 = SET32(0x03030303);

  // to do: faster implementation?!: 4bit then 8bit then 64 bit
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long ll=0; ll<individuals; ll++) {
    BlockType sum64 = ZERO(),
      *pC = ((BlockType0 *)(code + ll * lda));
    for (Long k=0; k<=kend; k++) {
      Long bend, qend;
      if (k) {
	bend = 1;
	qend = rest;
      } else {
	bend = blocks31;
	qend = x31;
      }
      for (Long b=0; b<bend; b+=x31) {
	BlockType sum = ZERO();      
	for (Long q=0; q<qend; q++, pC++) {
	  BlockType c = *pC;
	  for (Uint j=0; j<CodesPerByte; j++) { // for-schleife aufloesen ?
	    sum = ADD8(sum, AND(c, N6E2));
	    c = SHR32(c, 2);
	  }
	}
	sum64 = ADD64(sum64, SAD8(sum, zero));	
      }
    }
    cS[ll] = HORIZ_ADD64(sum64);
  }
  
  for (Long i=0; i<individuals; total += cS[i++]);
  
  if (cS != sums) { FREE(cS); }
  
  return total;
}
	
#endif  // SHUFFLE8



#define x127 (255 / 2)
static void allele_sum_I_2bit(unit_t *code, Long snps, Long individuals,
			       Long lda, Long LDAbitalign,
			      double *sumP, double *sumP2, Ulong *Pd) {
  if (4.0 * (double) snps * (double) individuals * (double) individuals
      > (double) LONG_MAX)
    ERR0("allele_freq exceeded maximum size of SNP matrix");
  
  Long
     lda127 = lda * x127; 
  Long
    blocks = Blocks(snps),
    byteblocks = (snps - 1L) / BytesPerBlock + CodesPerByte,
    // here: P= 2(p_i), not 2 p_i - 1 == expected values == means
    groupsM1 = (individuals - 1L) / x127,
    rest =  individuals - x127 * groupsM1,
    uniblocks = BytesPerUnit * byteblocks,// OK; close to "snps / UnitsPerBlock"
    bytesSums = uniblocks * BytesPerBlock + LDAbitalign / BitsPerByte,
    bytessum = byteblocks * BytesPerBlock + LDAbitalign / BitsPerByte;
  int 
    *Sums0 = (int*) MALLOC(bytesSums),
    *sum0 = (int*) MALLOC(bytessum),
     *lastSums0 = (int*) ((char*) Sums0 + bytesSums - BytesPerBlock),
    *lastsum0 = (int*) ((char*) sum0 + bytessum - BytesPerBlock);
  if (Sums0 == NULL || sum0 == NULL) ERR0("allocation error");
  const BlockType
    O6F2 = SET32(0x000000FF),
    N6E2= SET32(0x03030303);
  BlockType0
    zero = ZERO(),
    *sum = (BlockType0*) algn_generalL(sum0, LDAbitalign),    
    *Sums  = (BlockType0*) algn_generalL(Sums0, LDAbitalign);
  for (Long k=0; k<uniblocks; Sums[k++] = zero); 

  for (Long i=0; i<=groupsM1; i++) { // parallel funktioniert nicht gut!
    unit_t *pC0 = code + i * lda127;
    Long endq = i < groupsM1 ? x127 : rest;
    for (Long k=0; k < byteblocks; sum[k++] = zero);
   
    for (Long q=0; q<endq; q++) {
      BlockType0 *pC = (BlockType0 *) (pC0 + q * lda);
      BlockType *s0 = sum;
      for (Long k=0; k<blocks; k++, pC++) { // see permutation.tex (*)
	assert(pC < (BlockType0*)(code + lda * individuals-UnitsPerBlock));// OK
	BlockType c = LOAD(pC);
 	for (Uint j=0; j<CodesPerByte; j++, s0++) {
	  assert(s0 < (BlockType0*) lastsum0);
	  // achtung: in sum ist die SNP-Reihenfolge 0, 4, 8 ,...1, 5, ...
	  BlockType s =  LOAD(s0);
	  s = ADD8(s, AND(c, N6E2));
	  STORE(s0, s);
	  c =SHR32(c, 2);
	}
      }
    }

    // add Byte-Results up into 33 Bit-integer
    BlockType0 *pS0 = (BlockType0*) Sums;      
    for (Long k=0; k<byteblocks; k++) { // see permutation.tex (*)
      BlockType s0 = LOAD(sum + k);
      for (int j=0; j<BytesPerUnit; j++, pS0++) { // OK
	// achtung nochmals gemischt:
	assert(pS0 < (BlockType0*) lastSums0);
	BlockType pS = LOAD(pS0);
	pS = ADD32(pS, AND(s0, O6F2));
	STORE(pS0, pS);
	s0 = SHR32(s0, BitsPerByte);
      }
    }
  } // for groups
  FREE(sum0);

   
  Uint *S = (Uint*) Sums;			 
  Ulong  sum_Psq = 0,
    sum_P = 0;
  for (Long j=0; j<snps; j++) {
    Long block = j / CodesPerBlock,  
      k = j % CodesPerBlock, // 128er Bloecke
      // 32 Bytes jeweils erste 2 Bit; dann zweite 2 Bit; etc. --> l	
      l = k / CodesPerByte + BytesPerBlock * (k % CodesPerByte),
      // 1 32 Byte-Vektor: 0.,4.,8.,12.,16.20.,24.28, werden zuerst
      // gelesen; dann  1.,5.,9.,13.,17.21.,25.29 etc
      // o gibt diese Nummerierung an
      // (l / BytesPerBlock) gibt die Nummer des 32 Bytes-Vektors an       
      m = (l / BytesPerBlock) * BytesPerBlock, 
      o =(l % BytesPerBlock)/BytesPerUnit+(l % BytesPerUnit)*UnitsPerBlock;// OK
    Ulong s = S[ block * CodesPerBlock + m  +  o];
    
    Pd[j] = s;
    sum_P += s;
    sum_Psq += s * s;
  }
  
  FREE(Sums0);
  double n = (double) individuals;
  if (sumP != NULL) *sumP = (double) sum_P / n;
  if (sumP2 != NULL) *sumP2 = (double) sum_Psq / (n * n);
}




static void TwoBithaplo2genoIntern(unit_t *code, Long snps, Long individuals,
				   Long lda,
				   int VARIABLE_IS_NOT_USED cores,
				   unit_t *ans, Long ldAns) {
  // genotypes 0 coded as 00
  //           1          01
  //           2          10
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  Long blocks = Blocks(snps);
 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif	  
  for (Long j=0; j<individuals; j++) {
    BlockType0 *xx = (BlockType0*) (code  + j * lda),
      *aa = (BlockType0 *) (ans  + j * ldAns);
    for (Long b=0; b<blocks; b++) {
      const BlockType a = LOAD(xx + b),
	y = AND(a, E1N1), // linkes Bit 
	z = AND(a, N1E1); // rechtes Bit
      const BlockType
	yR1 = SHR32(y, 1),
	zL1 = SHL32(z, 1);
      const BlockType
	leftBit = AND(y, zL1),
	rightBit = XOR(yR1, z);
      STORE(aa + b, OR(leftBit, rightBit));
    }
  }
}


#if defined CODE_ONE
#undef CODE_LASTZERO
#undef CODE_ONE
#undef CODE_TWO
#endif
#define CODE_LASTZERO 0
#define CODE_ONE 1
#define CODE_TWO 2


#endif // 2bitIntern.h
