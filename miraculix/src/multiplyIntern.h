
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


#ifndef miraculix_multiplyIntern_H
#define miraculix_multiplyIntern_H 1


static BlockType CodeMask16;
#define MultiplyInitIntern(Version)			\
  void MultiplyInit##Version() {			\
    assert(BytesPerBlock == sizeof(BlockType0));	\
    SET16(CodeMask16, CodeMask);			\
  }


#define MULTIPLY							\
  AND(a, s1, CodeMask16); AND(b, s2, CodeMask16);  MADD16(dummy, a, b); \
  ADD32(sum.vi, sum.vi, dummy);

static Uint scalar012(BlockType0 *x, BlockType0 *y, Uint blocks) {
  BlockType0 *endx = x + blocks;
  BlockUnitType sum;
  ZERO(sum.vi);

  Uint zaehler = 0;

  for(; x < endx; x ++, y ++) {
    zaehler++;
    BlockType dummy, a, b, s1, s2;
    LOAD(s1, x);
    LOAD(s2, y);
       
    MULTIPLY;
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;   
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;
    SHR16(s1, s1, 2); SHR16(s2, s2, 2); MULTIPLY;
  }

 
  
  return sum.u32[0] + sum.u32[1] + sum.u32[2] + sum.u32[3]
#if defined AVX2
    + sum.u32[4] + sum.u32[5] + sum.u32[6] + sum.u32[7]
#endif    
    ;
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


#endif
