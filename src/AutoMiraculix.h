
#ifndef AutoMiraculix_H
#define AutoMiraculix_H 1
#include <R.h>

typedef enum relmatrix_methods {Shuffle,  // SSE2 or AVX needed
				TwoBit,   // none needed
				ThreeBit, // none needed
				Hamming2, // SSE2
				Hamming3, // SSE3
				NoSNPcoding,
				AutoCoding,    // none needed
				Haplo
                                } relmatrix_methods;

#define last_usr_meth AutoCoding
#define nr_relship_meth (Haplo + 1)

extern const char *RELSHIP_METH_NAME[nr_relship_meth];

// method shuffle; internal coding
// for haplo the information is finally doubled (Haplo in attr(method) & coding)
#define HAPLO 0
#define GENO 1
#define GENOMATRIX 2


// Coding of the attribute "information"
#define INFO_INFO 0
#define INFO_P 1
#define INFO_CODING 2
#define INFO_LAST INFO_CODING

#define WHAT 0
#define SNPS 1
#define INDIVIDUALS 2
#define ADDR0 3
#define ADDR1 4
#define ADDR2 5
#define ADDR3 6
#define MEM 7
#define SNPxIND 8 // INFO_INDIV_PER_COL, i.e. 'percolumn'
#define INFO_BLOCKS 9
#define INFO_BITS 10
#define INFO_BITSPERBLOCK 11
#define SNPSPERCOMPRESSED 12
#define INFO_HEADER 13
#define INFO_DOUBLEINDIV 14
#define INFO_LEADINGCOL 15
#define INFO_INFO_LAST INFO_LEADINGCOL

#define INFO_P_PPT 0
#define INFO_P_SUMPQ 1
#define INFO_P_SUMP 2
#define INFO_P_P 3

#endif
