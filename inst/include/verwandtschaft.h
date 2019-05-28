


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 Martin Schlather

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


#ifndef verwandtschaft_H
#define verwandtschaft_H 1

#include <inttypes.h>
#include "miraculix.h"
#include "AutoMiraculix.h"
#include "intrinsics.h"

#define ORIG_BLOCKS 4

typedef char table_type;
typedef uint16_t blockarray[ORIG_BLOCKS];

typedef union block_compressed {
  blockarray b;
  uint64 x;
} block_compressed;


#define nr_genotypes 3
#define nr_results 4


typedef SEXP (*coding_start_t)(Uint , Uint, SEXP);
typedef void (*coding_main_t)(Uint *, Uint, Uint, Uint, Uint, SEXP, double *);
typedef void (*coding_end_t)();

void printbits(uint64 x, Uint size, Uint bits);
void Init23();
#define Init Init23
SEXP get_matrix_start(Uint individuals, Uint snps, SEXP G);


void initiate_tableI(table_type **TABLE, Uint *TABLE_SIZE,
		     Uint blocklength, Uint bits, uint64 *result_code,
		     Uint *result_value, Uint NrResults);
void initiate_table3();
SEXP matrix_coding_start3(Uint individuals, Uint snps, SEXP G);
void matrix_coding3(Uint *SNPxIndiv, Uint start_individual, Uint end_individual,
		    Uint start_snp, Uint end_snp, SEXP Ans, double *G);
void MmPsq3(Uint individuals,
	    Uint min_row, Uint max_row, Uint min_col, Uint max_col,
	    Uint compressed_start, Uint compressed_end, Uint n_compressed,
	    uint64 *M, double *ergb, bool Parallel);
double getValue3(uint64 *Indiv, Uint s);
SEXP get_matrix3(SEXP SNPxIndiv);

void initiate_table2();
SEXP matrix_coding_start2(Uint individuals, Uint snps, SEXP G);
void matrix_coding2(Uint *SNPxIndiv, Uint start_individual, Uint end_individual,
		    Uint start_snp, Uint end_snp, SEXP Ans, double *G);
void MmPsq2(Uint individuals,
	    Uint min_row, Uint max_row, Uint min_col, Uint max_col,
	    Uint compressed_start, Uint compressed_end, Uint n_compressed,
	    uint64 *M, double *ergb, bool Parallel);
double getValue2(uint64 *Indiv, Uint s);
SEXP get_matrix2(SEXP SNPxIndiv);

#define Information Information23
extern SEXP Information;
#define Method Method23


SEXP CreateCode23(int individuals, int snps, int mem, int blocks,
		    int blocklength, int tablesize, SEXP file);
#define matrix_coding_start(TABLE_SIZE)	\
  (Uint individuals, Uint snps, SEXP file) {	 \
   assert(TWO_BITSPERBLOCK == two_bitsperblock); \
    Uint n_compressed = (Uint) CEIL(snps / (double) snpsPerCompressed),	\
      mem = n_compressed * individuals;					\
    SEXP Code = CreateCode23(individuals, snps, mem, blocks, blocklength, \
			     TABLE_SIZE, file),				\
      Infos = getAttrib(Code, Information);				\
    Uint *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));	\
    info[INFO_BITS] = bits;						\
    info[INFO_BITSPERBLOCK] = bitsperblock;				\
    info[SNPSPERCOMPRESSED] = snpsPerCompressed;			\
    return(Code);							\
  }


#endif
