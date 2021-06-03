
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2019 -- 2019  Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, writne to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#ifndef miraculix_packedIntern_H
#define miraculix_packedIntern_H 1

/*
static void showblock(Uint *X, Uint snps) { 
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



static void showblock(BlockType0 X) { showblock((Uint *) &X, CodesPerBlock); }
static void showblock(BlockType0 X, Uint blocks) {
  showblock((Uint *) &X, blocks * CodesPerBlock);
}

static void showblock(Uint *X, Uint snps, Uint individuals) { 
  for (Ulong i=0; i<individuals; i++) {
#ifdef SCHLATHERS_MACHINE    
    // PRINTF("%ld\n", (uint ptr_t) (X + i * Blocks(snps) * UnitsPerBlock));
#endif 
    showblock(X + i * Blocks(snps) * UnitsPerBlock, snps);
  }
}
*/
AVAILABLE


#define PackedInitIntern(Version)				\
  void PackedInit##Version() {			\
    assert(BytesPerBlock == sizeof(BlockType0));	\
    SET32(EinsMask4Bit, 0x11111111);			\
    SET32(CodeMask4Bit, 0x33333333);			\
    SET32(Null4Eins4,   0x0F0F0F0F);			\
  }



static Ulong sumGenoIntern(Uint *S, Uint snps, Uint individuals) { 
  Ulong ans = 0,
    unitsPerIndiv = UPI(snps),
    units = unitsPerIndiv  * individuals, // wird ueber alle
  // units hinweg summiert, ob in Gebrauch oder nicht, da
  // Null gesetzt falls nicht in Gebrauch. Einfacher als ueber
  // die echten Bloecke mit Fallunterscheidung.
    blocks = units / UnitsPerBlock;

#define maxint 2147483647  
#define x31 (255L / (2 * CodesPerByte)) // 2 : max{0,1,2}
#define GroupLoops (maxint / (x31 * BytesPerUnit * 2 * CodesPerBlock))
#define GroupLoopsM1 (GroupLoops - 1L)
#define BlocksPerGroup (GroupLoops * x31)

  Ulong
    intgroupsM1 = (blocks - 1L) / BlocksPerGroup,
    groupRestM1 = (blocks - 1L) / x31 - intgroupsM1 * GroupLoops,//BlocksPerGroup,
    rest = blocks - intgroupsM1 * BlocksPerGroup - x31 * groupRestM1;
  
  //    groupsM1 = (units - 1L) / x31;
  BlockType first2bits8, first8bits;
    // NullEins = SET32(0x55555555),
  SET32(first2bits8, 0x03030303);
  SET32(first8bits, 0x000000FF);
  BlockUnitType sum32;

  // to do: faster implementation?!: 4bit then 8bit then 64 bit
#ifdef DO_PARALLEL
  #pragma omp parallel for reduction(+:ans) num_threads(CORES)
#endif	  
  for (Ulong X=0; X<=intgroupsM1; X++) { // this loop needed only for huge
    // data sizes of more than 50000 individuals and 10^5 SNPS each
    ZERO(sum32 VI);
    BlockType *S_int = ((BlockType0 *) S) + X * BlocksPerGroup;
    Ulong groupsM1 = X < intgroupsM1 ? GroupLoopsM1 : groupRestM1;
    for (Ulong i=0; i<=groupsM1; i++) {
      // printf("X=%lu %lu i=%lu %lu\n", X, intgroupsM1, i, groupsM1);
      //if(X!=0)BUG;
      BlockType sum,
	*cgm = S_int + i * x31;
      ZERO(sum);
      // printf(" %d %d %d %d snps=%d %d\n", i, groupsM1, X , intgroupsM1, snps, individuals);
      Uint endx = i < groupsM1 || X < intgroupsM1 ? x31 : rest;
      for (Uint x=0; x<endx; x++, cgm++) {
	BlockType d,
	  c = *cgm;
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

  return ans;
}
	
static void haplo2genoIntern(Uint *X, Uint snps, Uint individuals,
			     Uint unitsPerIndiv,
			     Uint *Ans) {
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  Uint blocks = Blocks(snps);
  BlockType EN, NE;
  SET32(EN, 0xAAAAAAAA);
  SET32(NE, 0x55555555); 

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif	  
  for (Ulong j=0; j<individuals; j++) {
    Ulong jb = j * unitsPerIndiv;
    BlockType *x = (BlockType0*) (X  + jb),
      *a = (BlockType0 *) (Ans  + jb);
    for (Uint b=0; b<blocks; b++) {
      BlockType y,z,dummy;
      AND(y, x[b], EN);
      AND(z, x[b], NE);
      // to do: 2 further commands suffice: using 2bit arithmetics (&,&,>>1,add)
      //  SHL32(y, y, 1L); ADD32(a[b], y, z) // untested
      SHL32(dummy, z, 1L);
      AND(dummy, y, dummy);
      SHR32(y, y, 1);
      XOR(y, y, z);    
      OR(a[b], dummy, y); 
      // a[b] = OR(SHR32(y, 1L), OR(z, AND(y, SHL32(z, 1L))));
    }
  }
}


#if defined DO_FLOAT
#define real float

#else // ! DO_FLOAT
#define real double

#endif

real static inline *algnrl(real *X) {			       
  assert(algn_general(X, BytesPerBlock)>=(uintptr_t) X);
  return (real *)  algn_general(X, BytesPerBlock);
}	


#define VECTORADD					\
  LOAD(L1, xx); xx++;					\
  LOAD(L2, yy); yy++;					\
  SHR32(s, L2, 1);					\
  AND(s, s, L1); /* ( 1 * 2) >> 1*/				\
  SHR32(z, L1, 1);					\
  AND(z, z, L2); /* (2 * 1) >> 1 */					\
  OR(s, s, z); /* ( 1 * 2 + 2  * 1) >> 1 (uncleaned)*/			\
  AND(L1, L1, L2); /* 1 * 1 ; 2 * 2 (half) (uncleaned)*/		\
  SHR32(z, L1, 1); /* 1 * 1 ; 2 * 2 (half) (uncleaned)*/		\
  OR(s, s, z); /*  (1 * 2 + 2 * 1 + 2 * 2 (half)) >> 1 (uncleaned) */	\
  									\
  AND(z, s, EinsMask4Bit); ADD8(sum1A, sum1A, z);/*AND(sum1,s,CodeMask4Bit);*/ \
  SHR32(s, s, 2); AND(s, s, EinsMask4Bit); ADD8(sum1B, sum1B, s);	\
  AND(z, L1, CodeMask4Bit); ADD8(sum2A, sum2A, z);			\
  SHR32(L1, L1, 2); AND(L1, L1, CodeMask4Bit); ADD8(sum2B, sum2B, L1);

#define ZEROsum_1 ZERO(sum1A); ZERO(sum1B); 
#define ZEROsum_2 ZERO(sum2A); ZERO(sum2B);

#define BLOCKADD				\
  {						\
    BlockType L1, L2, z;					\
    VECTORADD; VECTORADD; VECTORADD; VECTORADD; VECTORADD;	\
    VECTORADD; VECTORADD;					\
  }


// 20 % of time; g
#define ADDUP4BIT_1						\
  SHR32(s, sum1A, 4); AND(s, s, Null4Eins4); ADD8(SB1, SB1, s);	\
  AND(sum1A, sum1A, Null4Eins4); ADD8(SB1, SB1, sum1A);		\
  SHR32(s, sum1B, 4); AND(s, s, Null4Eins4); ADD8(SB1, SB1, s);	\
  AND(sum1B, sum1B, Null4Eins4); ADD8(SB1, SB1, sum1B);	       

// c
#define ADDUP4BIT_2						\
  SHR32(s, sum2A, 4); AND(s, s, Null4Eins4); ADD8(SB2, SB2, s);	\
  AND(sum2A, sum2A, Null4Eins4); ADD8(SB2, SB2, sum2A);		\
  SHR32(s, sum2B, 4); AND(s, s, Null4Eins4); ADD8(SB2, SB2, s);	\
  AND(sum2B, sum2B, Null4Eins4); ADD8(SB2, SB2, sum2B)


static Uint scalar012(BlockType0 *x, BlockType0 *y, Uint blocks) {
  Uint ans = 0;
#define x28  (7 * 4) // (max=2) * 7 < 16; 4 = 255 / (2 * 7)
  
  Uint
    bundles = blocks / x28;

  BlockType
    *xx = x,
    *yy = y;
  for (Uint i=0; i<bundles; i++) { 
    BlockType SB1, SB2, sum1A, sum1B, sum2A, sum2B, s;
    ZERO(SB1);
    ZERO(SB2);
    ZEROsum_1; ZEROsum_2; BLOCKADD; ADDUP4BIT_2;    // 7-fold summation of SB2
    ZEROsum_2; BLOCKADD; ADDUP4BIT_1; ADDUP4BIT_2; // 14-fold summation of SB1
    ZEROsum_1; ZEROsum_2; BLOCKADD; ADDUP4BIT_2;
    ZEROsum_2; BLOCKADD; ADDUP4BIT_1; ADDUP4BIT_2;
    ans += 2 * ADDUPVECTOR(SB1) + ADDUPVECTOR(SB2); 
  }

  BlockType SB1, SB2, sum1A, sum1B, sum2A, sum2B, s;
  int end = blocks - bundles * x28 - 1;//auf keinen Fall Uint, da -1 moeglich !!
  ZERO(SB1);
  ZERO(SB2);
  ZEROsum_1; ZEROsum_2;
  for (int i=0; i <= end; i++) {//auf keinen Fall Uint, da -1 moeglich !!
    BlockType L1, L2, z;
    VECTORADD;
    if (i % 7 == 6 || i == end) {
      ADDUP4BIT_1; ADDUP4BIT_2;
      ZEROsum_1; ZEROsum_2;
    }
  }

  ans += 2 * ADDUPVECTOR(SB1) + ADDUPVECTOR(SB2);
  return ans;
}



// int core2_func (void) __attribute__ ((__target__ ("arch=core2")));
// https://dmalcolm.fedorapeople.org/gcc/2015-08-31/rst-experiment/declaring-attributes-of-functions.html#function-attributes
// # pragma GCC target ("string"...)
// https://dmalcolm.fedorapeople.org/gcc/2015-08-31/rst-experiment/pragmas-accepted-by-gcc.html

static void crossprodIntern(Uint *CGM, Uint snps, Uint individuals,
			    double *ans) {
  Uint
    blocks = Blocks(snps),
    unitsPerIndiv = UPI(snps);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(dynamic, 16)  
#endif
  for (Ulong i=0; i<individuals; i++) {
    BlockType0 *x = (BlockType0 *) (CGM + i * unitsPerIndiv);
    double *ans0 = ans + i,
      *ans1 = ans + i * individuals;
    for (Ulong j=i; j<individuals; j++) {
      Uint Mij = scalar012(x, (BlockType0 *) (CGM + j * unitsPerIndiv), blocks);
      ans0[j * individuals] = ans1[j] = (double) Mij;
    }
  }
}

   
static SEXP matrix_start_Intern(Uint snps, Uint individuals, SEXP VARIABLE_IS_NOT_USED file){
  SEXP Code;
  PROTECT(Code = CreateEmptyCodeVector(snps, individuals, MY_METHOD));
  // start_info(Code, file);
  UNPROTECT(1);
  return Code;
}




static void allele_freq_I(Uint *CGM, Uint snps, Uint individuals,
			 double *sumP, double *sumP2, real **Pd0) {
  // NOTE! twice allel frequency (indicated by the 2 at the end) 
  double n =  (double) individuals;
  if (individuals > 536870912)
    ERR("too many individualss to calculate allele frequency");
  Uint
    blocks = Blocks(snps),
    byteblocks = (snps - 1L) / BytesPerBlock + CodesPerByte,
    //snpBlocks = byteblocks * BytesPerBlock,
    //alignedSnpBlocks = snpBlocks + UnitsPerBlock,
    // here: P= 2(p_i), not 2 p_i - 1 == expected values == means
    groupsM1 = (individuals - 1L) / x31,
    rest =  individuals - x31 * groupsM1,
    unitsPerIndiv = UPI(snps),
    unitsPerIndiv31 = unitsPerIndiv * x31; 
  int
    //    *Sums0 = (int*) MALLOC(alignedSnpBlocks * BytesPerUnit),      
    *Sums0 = (int*) MALLOC(BytesPerUnit * byteblocks * BytesPerBlock
			   + BytesPerBlock),      
    *sum0 = (int*) MALLOC(byteblocks * BytesPerBlock + BytesPerBlock);// snps + zuschlag
  BlockType first8bits, first2bits8,
    *sum = (BlockType0*) algn_general(sum0, BytesPerBlock256),
    *Sums  = (BlockType0*) algn_general(Sums0, BytesPerBlock256);
  // NullEins = SET32(0x55555555),
  SET32(first8bits, 0x000000FF);
  SET32(first2bits8, 0x03030303);
  if (Sums == NULL || sum == NULL) ERR("allocation error");

  /// TODO
  /// 31 schleife hat ncihts mit "for (Uint x=0; x<endx; x++) " zu tun!!
  /// desahlt 31er schleife zu 127er!!!!
 
  { BlockType zero;
    ZERO(zero);
    Uint endfor = BytesPerUnit * byteblocks;
    for (Uint k=0; k<endfor; Sums[k++] = zero); //
  }
  
  for (Ulong i=0; i<=groupsM1; i++) { // parallel funktioniert nicht gut!
    BlockType *cgm = (BlockType0 *) (CGM + i * unitsPerIndiv31);
    Uint endx = i < groupsM1 ? x31 : rest;
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
	for (Uint j=0; j<CodesPerByte; j++, s++) { // for-schleife aufloesen?
	  AND(d, c, first2bits8);
	  ADD8(*s, *s, d);
	  SHR32(c, c, 2);
	}
      }

      //      printf("i=%d %d %d endx=%d\n", i, blocks, CodesPerByte, endx);

      //      BlockUnitType *z = (BlockUnitType0*) s;
      //      for (int k=0; k<32; k++) printf("%d ", z->u8[k]);
      
    }
  	 
    BlockType *pS = (BlockType0*) Sums;      
    for (Uint k=0; k<byteblocks; k++) { // see permutation.tex (*)
      BlockType s0 = sum[k];
      for (Uint j=0; j<BytesPerUnit; j++, pS++) {
	BlockType L1;
	AND(L1, s0, first8bits);
	ADD32(*pS, *pS, L1);

	//BlockUnitType *z = (BlockUnitType0*) pS;
	//	printf("%d %d %d %d\n", z->u32[0], z->u32[1], z->u32[2], z->u32[3]);
	
	SHR32(s0, s0, BitsPerByte);
      }
    }
  } // for groups
  FREE(sum0);
  
  Uint *S = (Uint*) Sums,
    alignedreal = byteblocks * BytesPerBlock * sizeof(real) + BytesPerBlock;
  (*Pd0) = (real*) MALLOC(alignedreal);
  real *Pd = (real *) algnrl(*Pd0);
  double
    sum_Psq = 0.0,
    sum_P = 0.0;
  if (OPTIONS_UTILS->basic.kahanCorrection) {
    ERR("Kahan not possible");
    double corr = 0.0,
      corrSq = 0.0;
    for (Uint j=0; j<snps; j++) {
      double Pi = (double) S[j] / n;
   	Pd[j] = (real) Pi;

	double fcorr = Pi - corr,
	  roughsum = sum_P + fcorr;
	corr = (roughsum - sum_P) - fcorr;	
	sum_P = roughsum;
	
	fcorr = Pi * Pi - corrSq;
	roughsum = sum_Psq + fcorr;
	corrSq = (roughsum - sum_Psq) - fcorr;	
	sum_Psq = roughsum;
    }
  } else {
    //    printf("\n");
    for (Uint j=0; j<snps; j++) {
      Uint block = j / CodesPerBlock,
	k = j % CodesPerBlock, // 128er Bloecke
	// 32 Bytes jeweils erste 2 Bit; dann zweite 2 Bit; etc. --> l	
	l = k / CodesPerByte + BytesPerBlock * (k % CodesPerByte),
	// 1 32 Byte-Vektor: 0.,4.,8.,12.,16.20.,24.28, werden zuerst
	// gelesen; dann  1.,5.,9.,13.,17.21.,25.29 etc
	// o gibt diese Nummerierung an
	// (l / BytesPerBlock) gibt die Nummer des 32 Bytes-Vektors an       
	m = (l / BytesPerBlock) * BytesPerBlock,
	o = (l % BytesPerBlock) / BytesPerUnit + (l % BytesPerUnit) * UnitsPerBlock;
	double Pi = (double) S[ block * CodesPerBlock + m  +  o] / n;

	Pd[j] = (real) Pi;
      sum_P += Pi;
      sum_Psq += Pi * Pi;
    }
  }
    
  FREE(Sums0);
  if (sumP != NULL) *sumP = sum_P;
  if (sumP2 != NULL) *sumP2 = sum_Psq;
}


static SEXP allele_freqIntern(SEXP GM) {
  SEXP Ans = R_NilValue; 
  Uint 
    *info = GetInfo(GM),
    individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    *cgm = Align256(GM, ALIGN_ALLELE);
  real *Pd0 = NULL;
  if (!isGeno(info[METHOD])) ERR("not a coded Z matrix");
  PROTECT(Ans = allocVector(REALSXP, snps));
  allele_freq_I(cgm, snps, individuals, NULL, NULL, &Pd0); 
  real *Pd = (real*) algnrl(Pd0);
  double *ans = REAL(Ans);
  for (Uint i=0; i<snps; i++) ans[i] = 0.5 * (double) Pd[i];
  FREE(Pd0);
  UNPROTECT(1);
  return Ans;
}



#endif // packedIntern.h
