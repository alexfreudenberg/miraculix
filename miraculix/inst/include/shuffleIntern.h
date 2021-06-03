
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


#ifndef miraculix_shuffleIntern_H
#define miraculix_shuffleIntern_H 1

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

#define ShuffleInitIntern(Version)					\
  void ShuffleInit##Version() {					\
    assert(BytesPerBlock == sizeof(BlockType0));			\
     /*                   0  1  2  3  4  5  6  7   8  9 10 11 12 13 14 15 */ \
    /* #define andshuffle 0, 1, 4, 0, 1, 2, 5, 0,  4, 5, 8, 0, 0, 0, 0, 0 */ \
    /* #define orshuffle 0, 0, 0, 2, 0, 0, 0, 2,  0, 0, 0, 2, 2, 2, 2, 4 */ \
    /*    SETREV8(as[i], andshuffle); */				\
    /*    SETREV8(xs[i], orshuffle); */					\
    SETREV8(and2scalar, 0, 1, 4, 0, 1, 2, 5, 0,  4, 5, 8, 0, 0, 0, 0, 0); \
    SETREV8(or2scalar, 0, 0, 0, 2, 0, 0, 0, 2,  0, 0, 0, 2, 2, 2, 2, 4); \
    SET32(LH, 0x0F0F0F0F);						\
  }

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
  //  XOR(z, xx, yy);
  OR(z, xx, yy);
  AND(y, z, LH);
  SHUFFLE8(y, or2scalar, y);
  ADD8(sum, sum, y);

  SHR32(y, z, 4L);
  AND(y, y, LH);
  SHUFFLE8(y, or2scalar, y);
#else  // FOLGENDES IST NICHT VIEL SCHNELLER & HAT FEHLER!!!
  XOR(z, xx, yy);
  SHR32(y, z, 1L);
  AND(z, y, z);
  SHR32(y, z, 3L);
  OR(z, z, y);
  AND(y, z, LH);
  SHUFFLE8(y, or2scalar, y);
#endif

  ADD8(sum, sum, y);
  return sum;
}


static Uint scalar012(BlockType0 *x, BlockType0 *y, Uint blocks) {
  Uint ans = 0;
  //#define x15 (255L / 16L)  // 16L = 4 * 2^2 (a byte full of 1's)
#define x15 (127L / 16L)  // warum ist dies schneller ?????
  Uint
    bundles = blocks  / x15;

  for (Uint i=0; i<bundles; i++) { 
    BlockType L1, L2, sum,
      *xx = x + i * x15,
      *yy = y + i * x15;
    ZERO(sum);
#if x15 == 7L
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));

    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
    LOAD(L1, xx); xx++; LOAD(L2, yy); yy++; ADD8(sum, sum, VECTORADD(L1, L2));
#elif x15 == 15L
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
  return ans;
}


static void crossprodIntern(Uint *CGM, Uint snps, Uint individuals,
			    double *ans) {
  Uint
    blocks = Blocks(snps),
    unitsPerIndiv = UPI(snps);

  // printf("blocks=%d\n", blocks);
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



#endif // shuffleIntern.h
