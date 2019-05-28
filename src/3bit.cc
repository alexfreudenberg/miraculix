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


// gcc -mavx -o hello_avx hello_avx.c

// _MM_ALIGN16 Uint ccc;

#include <stdio.h>
#include "options.h"
#include "Miraculix_aux.h"
#include "verwandtschaft.h"
#include <Basic_utils.h>

#ifdef DO_PARALLEL
#include <omp.h>
#endif

#define bits 3
#define blocklength 5
#define TWO_BITSPERBLOCK 32768 // 2^(bits * blocklength)

static double rev_geno_code[7] = {0, RF_NA, RF_NA, 1, RF_NA, RF_NA, 2}; 
static uint64 geno_code[nr_genotypes] = {0, 3, 6}, // 0, 1, 2
  result_code[nr_results] = {0, 3, 2, 6};
static Uint result_value[nr_results] = {0, 1, 2, 4};

static Uint bitsperblock = bits * blocklength,
  two_bitsperblock = (Uint) POW(2, bitsperblock),
  bitpattern = (1 << bits) - 1, // >>
  nextpower2 = (1 << (Uint) CEIL(LOG(bitsperblock) / LOG2)),
  blocks = sizeof(uint64) * bitsperbyte / nextpower2,
  deltaPower2BitsPB = nextpower2 -  bitsperblock,
  bitspercompressed = nextpower2 * blocks,
  snpsPerCompressed = blocks * blocklength;
			

// Codierung: 0->0; 1->3; 2->6
static table_type *TABLE3 = NULL;
static Uint TABLE3_SIZE = 0;
void initiate_table3() { // hash table 
  if (TABLE3 != NULL) BUG; 
  initiate_tableI(&TABLE3, &TABLE3_SIZE, blocklength, bits, result_code,
		  result_value, nr_results);
}

/*
void static printbits(uint64 x) { //
  printbits(x, sizeof(uint64) * bitsperbyte, bits);//
}
*/

// these defintions must match extactly as.genomicmatrix
SEXP matrix_coding_start3 matrix_coding_start(TABLE3_SIZE)

void matrix_coding3(Uint *M, Uint start_individual, Uint end_individual, 
		    Uint start_snp, Uint end_snp,
		    SEXP Ans, double VARIABLE_IS_NOT_USED * G) {
  // always assuming that SNPxIndiv is snps x individuals
  SEXP Infos = getAttrib(Ans, Information);
  Uint *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
 
  if (start_snp % snpsPerCompressed != 0) BUG;
  Uint 
    *pM = M,
    individuals = info[INDIVIDUALS],
    n_compressed = info[MEM] / individuals;
  uint64
    *ans = (uint64*) REAL(Ans)+ start_snp / snpsPerCompressed;
  double
    *p = REAL(VECTOR_ELT(Infos, INFO_P)) + INFO_P_P;
    

  for (Uint s=start_snp; s<end_snp; p[s++] = 0.0);
  for (Uint a=start_individual; a<end_individual; a++) { // calculate SNP means and compress the data
    Uint shift = 0,
      counter = 0;
    uint64 compressed = 0,
      *pAns = ans + a * n_compressed;

    for (Uint s=start_snp; s<end_snp; s++, pM++) {
      //    printf(">> %d ", *pM);
  
      p[s] += (double) (*pM);
      //printf("%d %d\n", s, shift);
      compressed |= geno_code[*pM] << shift;
      shift += bits;
      if (++counter >= blocklength) {
	shift += deltaPower2BitsPB;
	counter = 0;
      }
      if (shift >= bitspercompressed) {
	//	printbits(compressed);
	//	print("%ld a=%d %d %d\n", (pAns - M), a, n_compressed, s );
	*pAns = compressed;
	pAns++;
	shift = counter = 0;
	compressed = 0;
      }
    }
    if (shift > 0) {
      //          	printbits(compressed);print("\n");
      //  print("%ld \n", (pAns - M) );
      *pAns = compressed;
      pAns++;
    }
    // printbits(compressed);
    //  printf("delta=%ld, a=%d n==%d f=%d (snps=%d bitc=%d bl=%d dP2=%d)\n", 
    //	   pAns - M, a , n_compressed, (a + 1) * n_compressed,
    //	   snps, bitspercompressed, blocklength, deltaPower2BitsPB);
    // assert(pAns - M == (a + 1) * n_compressed);
  }
}



double getValue3(uint64 *Indiv, Uint s) {
  block_compressed x;
  x.x =  Indiv[s / snpsPerCompressed];
  Uint idx = s % snpsPerCompressed,
    y = x.b[idx / blocklength] >> ( (idx % blocklength) * bits); 
  //  printf("y=%d \n", (Uint) y);
  // printf("%10e\n", rev_geno_code[y & bitpattern]);
  //  printf("ok\n");

  return rev_geno_code[y & bitpattern];
}


 
void MmPsq3(Uint individuals, Uint min_row, Uint max_row,
	    Uint min_col, Uint max_col,
	    Uint compressed_start, Uint compressed_end, Uint n_compressed,
	    uint64 *M, double *ergb, bool Parallel) {
  table_type *table = TABLE3;
  assert(//min_row >= 0 && min_col >= 0 &&
	 max_row <= individuals && max_col <= individuals);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) if (Parallel)
#endif	  
  for (Uint c=min_col; c<max_col; c++) {
    Uint r  =  c > min_row ? c : min_row;
    uint64 * Mc = M + c * n_compressed,
      *Mr = M + r * n_compressed;
    double *ptr = ergb + c * individuals,
      *res_ptr_r = ptr + r;
    for ( ; r<max_row; r++, res_ptr_r++, Mr += n_compressed) {
      Uint sum = 0;
      uint64 *MR = Mr + compressed_start,
	*MC = Mc + compressed_start;
      for (Uint s=compressed_start; s<compressed_end; s++, MR++, MC++) {
	block_compressed x;
	x.x = *MR & *MC;
	sum+= table[x.b[0]] + table[x.b[1]] + table[x.b[2]] + table[x.b[3]];
      }
      ergb[c + r * individuals] = *res_ptr_r = (double) sum;
    }
  }
}

SEXP get_matrix3(SEXP SNPxIndiv) {
  SEXP
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  Uint 
    individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    n_compressed = info[MEM] / individuals,
    endfor = individuals * n_compressed;
  SEXP Ans = get_matrix_start(individuals, snps, 0);
  Uint *ans = (Uint *) INTEGER(Ans);
  uint64 *M = (uint64 *) REAL(SNPxIndiv);
  for (Uint a=0; a<endfor; a+=n_compressed) {
    uint64 *Ma = M + a;
    for (Uint s=0; s<snps; s++) *(ans++) = getValue3(Ma, s);
  }
  return Ans;
}
