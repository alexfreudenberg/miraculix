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


// gcc -mavx -o hello_avx hello_avx.c

// _MM_ALIGN16 Uint ccc;
#define BitsPerCode 3L
#define MY_METHOD ThreeBit

#include "Bit23.intern.h"
#include "haplogeno.h"
#include "Haplo.h"
#include "align.h"
#include "xport_import.h"

#define TWO_GENUINEBITSPERMINIBLOCK 32768  // testing only
#define TABLE_SIZE  two_genuineBitsPerMiniblock
#define nr_results 4

static double rev_geno_code[7] = {0, RF_NA, RF_NA, 1, RF_NA, RF_NA, 2}; 
static BlockType geno_code[nr_genotypes] = {0, 3, 6}, // 0, 1, 2
  result_code[nr_results] = {0, 3, 2, 6}; //
static Uint result_value[nr_results] = {0, 1, 2, 4};

// Codierung: 0->0; 1->3; 2->6
static table_type *TABLE3 = NULL;

/*
void static printbits(BlockType x) { //
  printbits(x, sizeof(BlockType) * BitsPerByte, BitsPerCode);//
}
*/

Ulong static inline genuineBlocks(Ulong snps) {				
  return 1L + (snps - 1L) / genuineCodesPerBlock; }				

Uint BytesPerBlock3() { return BytesPerBlock; }
Uint CodesPerBlock3() { return genuineCodesPerBlock; }
Uint UnitsPerIndiv3(Uint snps) { return genuineBlocks(snps) * UnitsPerBlock; }
#define UPI UnitsPerIndiv3
Uint BitsPerCode3() { return BitsPerCode; }

void Init3() {
  if (MiniblocksPerBlock != 4) BUG;
  assert(BytesPerBlock == sizeof(BlockType0));
  initiate_tableI(&TABLE3, TABLE_SIZE, CodesPerMiniblock, BitsPerCode,
		  result_code, result_value, nr_results);
}


// these defintions must match extactly as.genomicmatrix
SEXP matrix_start3 START23(Init3)

  
Ulong sumGeno3(Uint *S, Uint snps, Uint individuals) {	
  Ulong	unitsPerIndiv = UPI(snps),					
    sum = 0L;							
  Uint counter = 0;						
  for (Uint i=0; i<individuals; i++, S += unitsPerIndiv) { /* ok */	
    for (Ulong j=0; j<unitsPerIndiv; j++) {				
      Uint s = S[j];						
      for (Uint u=0; u<CodesPerUnit; u++) {			
	sum += rev_geno_code[s & CodeMask];
	s >>= BitsPerCode;					
	if (++counter >= CodesPerMiniblock) {			
	  s >>= deltaMiniblockBits;				
	  counter = 0;						
	}							
      }								
    }							  
  }							  
  return sum;						  
  }

  

void haplo2geno3(Uint *X, Uint snps, Uint individuals, Uint unitsPerIndiv,
		 Uint *ans) {
  Coding23(0, individuals, 0, snps, unitsPerIndiv, FROMHAPLO);
}



void coding3(Uint *X, Uint start_individual, Uint end_individual,		
	     Uint start_snp, Uint end_snp, Uint Mnrow,				
	     SEXP Ans, double VARIABLE_IS_NOT_USED * G) {/*closing in .cc file*/
  if (start_snp % genuineCodesPerBlock != 0) BUG;			
  Uint *info = GetInfo(Ans),						
    snps = info[SNPS],
    unitsPerIndiv = UPI(snps),
    *ans = Align(Ans, ALIGN_23) +
         start_snp * UnitsPerBlock / genuineCodesPerBlock;
  									
Coding23(start_individual, end_individual, start_snp, end_snp, Mnrow,
	 FROMINPUT);
 
} // see bit23intern.h



void crossprod3(Uint * M, Uint snps, Uint individuals, double *A) {
  if (TABLE3 == NULL) Init3();
  Uint blocks = genuineBlocks(snps),
    unitsPerIndiv = UPI(snps);
  table_type *table = TABLE3;
  assert(table != NULL);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) 
#endif	  
  for (Uint i=0; i<individuals; i++) {   
    Uint * Mi = M + i * unitsPerIndiv,
      *Mj = Mi;
    
    double *ans = A + i,
      *ptr = A + i * individuals,
      *res_ptr_r = ptr + i;
    for ( Uint j  =  i; j<individuals; j++, res_ptr_r++, Mj += unitsPerIndiv) {
      Uint sum = 0;
      BlockType *MJ = (BlockType0*) Mj,
	*MI = (BlockType0*) Mi;
      for (Uint s=0; s<blocks; s++, MJ++, MI++) {
	block_compressed x;
	x.x = *MJ & *MI;
	sum+= (table[x.b[0]] + table[x.b[1]] + table[x.b[2]] + table[x.b[3]]);
      }
      ans[j * individuals] = *res_ptr_r = (double) sum;
    }
  }
}


double getValue3(BlockType0 *MC, Uint s) {
  uint16_t x = ((uint16_t*) MC)[s / CodesPerMiniblock];
  Uint idx = s % CodesPerMiniblock; 
  //  printf("y=%d %d \n", genuineCodesPerBlock, CodesPerMiniblock);
  // printf("%10e\n", rev_geno_code[y & CodeMask]);
  //  printf("ok\n");

  return rev_geno_code[(x >> (BitsPerCode * idx)) & CodeMask];
}


SEXP get_matrix3(SEXP SNPxIndiv) {
  if (TABLE3 == NULL) Init3();
  Uint 
    *info = GetInfo(SNPxIndiv),
    individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    unitsPerIndiv = UPI(snps),
    *M = Align(SNPxIndiv, ALIGN_23);
  SEXP Ans;
  PROTECT(Ans = get_matrix23_start(snps, individuals, R_NilValue));
  Uint *ans = (Uint *) INTEGER(Ans);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif	  
  for (Uint i=0; i < individuals; i++) {
    Uint *A = ans + i * snps;
    BlockType *Ma = (BlockType0*) (M + i * unitsPerIndiv);
    for (Uint s=0; s<snps; s++) A[s] = getValue3(Ma, s);
  }
  UNPROTECT(1);
  return Ans;
}




SEXP get_matrixN_3(SEXP SNPxIndiv, SEXP Snps) {
  if (TABLE3 == NULL) Init3();
  Uint *M =  Align(SNPxIndiv, ALIGN_23); 
  GET_MATRIX_N_23(getValue3, UPI);
}



SEXP allele_freq3(SEXP SNPxIndiv) {
 if (TABLE3 == NULL) Init3();
 Uint *M = Align(SNPxIndiv, ALIGN_23);
 ALLELE_FREQ23(getValue3, UPI);
}


Uint *Align3(SEXP Code, Uint nr, bool test) { return AlignTest(Code, nr, test); }

void zeroNthGeno3(SEXP SNPxIndiv, SEXP Snps) {					
  Uint *M = Align(SNPxIndiv, ALIGN_23);
  ZERONTH(UPI);
}
