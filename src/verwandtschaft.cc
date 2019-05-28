
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


// _MM_ALIGN16 Uint ccc;

#include "verwandtschaft.h"
#include "options.h"
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "general.relmatrix.h"

#include "error.h"

char EOL = '\n';

void initiate_tableI(table_type **TABLE, Uint *TABLE_SIZE,
		     Uint blocklength, Uint bits, uint64 *result_code,
		     Uint *result_value, Uint NrResults) { // hash table

  Uint d,
    *nx;
  nx = (Uint *) CALLOC(blocklength, sizeof(Uint));
  unsigned long int size = (unsigned long int) POW(2, bits * blocklength);
  
  *TABLE = (table_type*) CALLOC(size, sizeof(table_type));
  *TABLE_SIZE = size * sizeof(table_type);
  table_type* table = *TABLE;

  while (true) {
    uint64 value = 0;
    Uint shift,
      sum = 0;
    
    for (d=shift=0; d<blocklength; d++, shift+=bits) {
      value |= result_code[nx[d]] << shift;
      sum += result_value[nx[d]];
    }
      
    assert(value < size && table[value] == 0);
    table[value] = sum;
    d = 0;
    nx[d]++;
    while (nx[d] >= NrResults) {
      nx[d] = 0;
      if (++d >= blocklength) break;
      nx[d]++;
    }
    if (d >= blocklength) break;
  }
  FREE(nx);
}



SEXP Information = R_NilValue,
  Method = R_NilValue;
void Init23() {
  //  assert(BytesPerBlock == sizeof(BlockType0));
  // assert(BytesPerUnit == sizeof(Uint));
  assert(sizeof(int) == sizeof(Uint));
  assert(sizeof(uint64) == sizeof(double)); 
  if (sizeof(uint64) != 8) ERR("system does not support 'long int'");
  
  if (Information != R_NilValue) BUG;
  Information = install("information");
  Method = install("method");

  initiate_table3();
  initiate_table2();
}



#define UPDATE								\
  if (isSNPxInd) {							\
    Uint jj, rP1 = r + 1;						\
    if (rP1 % nrow_matrix == 0) {					\
      /* always delivering matrix as snps x individuals */		\
      coding_main(matrix, 0, individuals, rP1-snpsPerCompressed, rP1, Ans,dG); \
      rowidx = -plusrow;						\
      for (jj=0; jj<matrix_size; matrix[jj++] = 0);			\
    }									\
  } else {								\
    Uint jj;								\
    coding_main(matrix, r, r + 1, 0, snps, Ans, dG);			\
    rowidx = -plusrow;							\
    for (jj=0; jj<matrix_size; matrix[jj++] = 0);			\
  }


#define CLOSE								\
  if (isSNPxInd) {							\
    Uint remainder = r % nrow_matrix;					\
    if (remainder != 0) {						\
      coding_main(matrix, 0, individuals, (Uint) (r - remainder), r, Ans, dG); \
    }									\
  } /* KEIN else  */							\
									\
  fclose(fp);								\
  FREE(matrix);								\
  if (PL > 1) { PRINTF("data have been read.");				\
    if (!true) {							\
      Uint j0,k0,z0=0;							\
      for (k0=0; k0<individuals; k0++) {				\
	for(j0=0; j0<n_compressed; j0++) {				\
	  PRINTF("%lu ", (long unsigned int) REAL(Ans)[z0++]);		\
	}								\
	PRINTF("\n");							\
      }									\
    }									\
  }									\
  return Ans;
  

SEXP file_binary_intern(char *file, 
			char *coding, Uint header, bool isSNPxInd, 
			bool doubledindividuals,
			// output
			Uint individuals,
			Uint snps,
			Uint snpsPerCompressed,
			Uint n_compressed,
			coding_start_t coding_start,
			coding_main_t coding_main,
			SEXP G
			) {
  // for plinkbinary files

  char ch, 
    //   A = coding[0],
    B = coding[1],
    C = coding[2]
    //  NA = coding[3],
    //  SEP = coding[4]
    ;
  bool
    haplo = B == C;

  if (haplo) ERR("haplo not recognized for binary files");
   
  FILE *fp;
  Uint j, nrow_matrix, ncol_matrix, perbyte,
    *matrix, // snps x individuals
    rows, // rowsM1, 
    cols,  r, l, i, idx,
    plusrow, rowidx, hr, hc, haplorow, haplocol;
  
  if (snps <= 0 || individuals <= 0 || snps <= individuals) {
    char msg[200];
    SPRINTF(msg, "matrix is supposed to be a '%.20s' x '%.20s' matrix with %d individuals and %d SNPs, which is odd.\n", 
	    isSNPxInd ? "SNPs" : "individuals",
	    isSNPxInd ? "individuals" : "SNPs",
	    individuals, snps);
    if (snps > 0 && individuals > 0) { warn(msg) } else ERR(msg);
  }
  
  if (isSNPxInd) {
    plusrow = 1;
    nrow_matrix = snpsPerCompressed; 
    rows = snps;
    cols = ncol_matrix = individuals;
  } else {
    // plusrow = snps;
    plusrow = snpsPerCompressed;
    rows = individuals;
    cols = ncol_matrix = snps;
    nrow_matrix = 1;
  }
  if (haplo) {
    if (isSNPxInd xor !doubledindividuals) {
      haplocol = 2;
      haplorow = 1;
    } else {
      haplocol = 1;
      haplorow = 2;    
    } 
  } else {
    haplocol = haplorow = 1;
  }

  SEXP Ans = coding_start(individuals, snps, G);
  Uint matrix_size = nrow_matrix * (ncol_matrix + 3);
  if ((matrix = (Uint*) CALLOC(matrix_size, sizeof(Uint))) == NULL)
    ERR("memory space could not be acquired");
  double *dG = REAL(G);

  cols = (Uint) CEIL(0.25 * cols);
  perbyte = 4 / haplocol;
  char twobitpattern = 3;

  //   printf("plusrow %d %d , %d %d c=%d %d haplorow/col %d %d\n", plusrow, nrow_matrix, individuals, haplo, cols, isSNPxInd, haplorow, haplocol);  //assert(false);

  if ((fp = fopen(file, "r")) == NULL)  {
    ERR1("file '%.50s' could not be opened", file);
  }
  for (i=0; i<header; i++) fgetc(fp);

  for (rowidx=r=0; r<rows; r++, rowidx+=plusrow) { // number of (multiple) rows 
    for (hr=0; hr<haplorow; hr++) { 
      // genuine loop only if haplo types are given rowwise
      idx = rowidx;

      //     printf("\n"); //assert(r < 10);

      for (l=0; l<cols; l++) { // number of (multiple cols 
	Uint k = 0;
	ch = fgetc(fp);
	for (j = 0; j < perbyte; j++, idx+=nrow_matrix) {
	  for (hc=0; hc<haplocol; hc++, k+=2) {
	    char value = (ch >> k) & twobitpattern;
	    //	  for (hc=0; hc<haplocol; hc++, ch >>= 2) {
	    // char value = ch & twobitpattern;
	    // genuine loop only if haplo types are given column wise
	    //  printf("%d", value==3 ? 2 : value==2 ? 1 : value == 0 ? 0 : value);

	    if (value != 0) { // == 0
	      // print("idx=%d %d (%d %d) %d\n", idx, matrix_size,  nrow_matrix, ncol_matrix, r);
	      assert(idx <= matrix_size);
	      if (value == 2) matrix[idx]++;
	      else if (value == 3) matrix[idx] += 2;
	      else if (value == 1) ERR("missing value detected.")
	      else {
		ERR8("unknown character detected (>%c<). Note that missings are not allowed here.\n Debugging information: row=%d of %d rows, col=%d of %d columns; plus=%d; current haplocol=%d of %d\n", ch, r, rows, l, cols, nrow_matrix, hc, haplocol);
	      }
	    //assert(matrix[idx] <=2);
	    } // ch != A
	  } //hc
	} // j
      } // cols l
    } // hr

    UPDATE;    

  } // rows r
 
 CLOSE
  
}



SEXP file_coding_intern(SEXP file,
			coding_start_t coding_start,
			coding_main_t coding_main, SEXP G) {
  
  SEXP Infos = getAttrib(file, Information);
  Uint *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
 			
  int
    header = info[INFO_HEADER];
  Uint
    leadingcols = info[INFO_LEADINGCOL],
    isSNPxInd = info[SNPxIND],
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    n_compressed = info[MEM] / individuals,
    doubledindividuals = info[INFO_DOUBLEINDIV],
    snpsPerCompressed = info[SNPSPERCOMPRESSED];
 
  char *coding = (char*) CHAR(STRING_ELT(VECTOR_ELT(Infos, INFO_CODING), 0));
  char *filename = (char*) CHAR(STRING_ELT(file, 0));

  if (header < 0) {
    // only for plinkbinary: snps & individuals already known
    return file_binary_intern(filename, coding,
			      -header, isSNPxInd, doubledindividuals,
			      individuals,
			      snps, snpsPerCompressed, n_compressed,
		       // output
			      coding_start, coding_main, G
			      );
  } // else snps not known yet
  

  // BEAGLE 3 
  // beides: ani x snps und snps x individ
  // A B A B ?
  // ? als voraussetzung nicht vorhanden
  // 1 header-Zeile
  // 2 erste Spalten weg
  //
  // plink: kein header
  // 4 spalten, die nicht brauchen
  char ch, oldch,
    A = coding[0],
    B = coding[1],
    C = coding[2],
    NA = coding[3],
    SEP = coding[4]
    ;
  bool 
    haplo = B == C;

  // printf("%.50s|\n", *coding);  assert(!haplo);
   
  FILE *fp;
  Uint nrow_matrix, ncol_matrix,
    *matrix = NULL, // snps x individuals
    rows, // rowsM1, 
    cols, colsM1,
    plusrow,  haplorow, haplocol;
  if ((fp = fopen(filename, "r")) == NULL) {
    ERR1("file '%.50s' could not be opened", filename);
  }

  // determine number of rows and columns
  for (int i=0; i<header; i++) while (fgetc(fp) != EOL);

  // determine size of matrix:
  rows = 1;
  cols = 0;
  oldch = SEP;
  while ((ch  = fgetc(fp)) != EOL) {  
    if (ch == SEP && oldch != SEP) cols++;
    oldch = ch;
  }
  if (oldch != SEP) cols++;
  cols -= leadingcols;

  while ((ch  = fgetc(fp)) != EOF) { 
    if (ch == EOL) rows++;
    oldch = ch;
  }
  if (oldch != EOL) rows++;
  fclose(fp);
  
  if (isSNPxInd) {
    individuals = cols;  
    snps = rows;
  } else {  
    individuals =  rows;
    snps = cols;    
  }
  
  if (haplo) {
    if (doubledindividuals) {
      Uint dummy = individuals / 2; 
      if (dummy * 2 != individuals) ERR("odd number of values (individuals)");
      individuals = dummy;
    } else {
      Uint dummy = snps / 2; 
      if (dummy * 2 != snps) ERR("odd number of values (SNPs)");
      snps = dummy;
    }
  }

  info[INDIVIDUALS] = individuals;
  info[SNPS] = snps;
 
  if (snps <= 0 || individuals <= 0 || snps <= individuals) {
    char msg[200];
    SPRINTF(msg, "matrix is supposed to be a '%.20s' x '%.20s' matrix with %d individuals and %d SNPs, which is odd.\n", 
	    isSNPxInd ? "SNPs" : "individuals",
	    isSNPxInd ? "individuals" : "SNPs",
	    individuals, snps);
    if (snps > 0 && individuals > 0) { warn(msg) } else ERR(msg);
  }

  colsM1 = cols - 1;
  //rowsM1 = rows - 1;

  // read in table
  
  if (isSNPxInd) {
    plusrow = 1;
    nrow_matrix = snpsPerCompressed; 
    ncol_matrix = individuals;
  } else {
    // plusrow = snps;
    plusrow = snpsPerCompressed;
    ncol_matrix = snps;
    nrow_matrix = 1;
  }
  if (haplo) {
    if (isSNPxInd xor !doubledindividuals) {
      haplocol = 2;
      haplorow = 1;
    } else {
      haplocol = 1;
      haplorow = 2;    
    } 
  } else {
    haplocol = haplorow = 1;
  }

  // now, code matrix
  SEXP Ans = coding_start(individuals, snps, G);
  Uint matrix_size = nrow_matrix * ncol_matrix;
  if ((matrix = (Uint*) CALLOC(matrix_size, sizeof(Uint))) == NULL)
    ERR("memory space could not be acquired");
  double *dG = REAL(G);

  rows /= haplorow;
  cols /= haplocol;

  //   printf("plusrow %d %d, %d %d c=%d %d haplorow/col %d %d\n", plusrow, nrow_matrix, individuals, haplo, cols, isSNPxInd, haplorow, haplocol);  //assert(false);

  fp = fopen(filename, "r");

  // jump header lines
  for (int i=0; i<header; i++) while ((ch=fgetc(fp)) != EOL);

  Uint r=0;
  for (Uint rowidx=0; r<rows; r++, rowidx+=plusrow) {//# of (multiple) rows
    for (Uint hr=0; hr<haplorow; hr++) { 
      // genuine loop only if haplo types are given rowwise
      Uint idx = rowidx;

      //   printf("0\n"); //assert(r < 10);
      
      // jump unnecessary leading stuff
      while ((ch = fgetc(fp)) == SEP);
      for (Uint l=0; l<leadingcols; l++) {
	while ((ch=fgetc(fp)) != SEP && ch!=EOF);
	if (ch == EOF) ERR("unexpected end of file");
	while ((ch = fgetc(fp)) == SEP);
      }
      
      for (Uint l=0; l<cols; l++, idx+=nrow_matrix) { // # of (multiple cols 
	for (Uint hc=0; hc<haplocol; hc++) {
	  // genuine loop only if haplo types are given column wise
	  
	  if (ch != A) { // == 0
	    //  if (idx>=nrow_matrix * ncol_matrix) 
	    // print("idx=%d %d, nrow_matrix=%d\n", idx, nrow_matrix * ncol_matrix, nrow_matrix);
	    if (ch == B) matrix[idx]++;
	    else if (ch == C) matrix[idx] += 2;
	    else if (ch == NA) ERR("missing value detected.")
	    else {
	      PRINTF(">%c< row=%d of %d rows, col=%d of %d columns; plus=%d hc=%d haplocol=%d\n", ch, r, rows, l, cols, nrow_matrix, hc, haplocol);
	      ERR("Unknown character detected. Note that missings are not allowed here.");
	    }
	    //assert(matrix[idx] <=2);
	  }
	  if (l<colsM1 || hc!=haplocol-1)
	    while((ch = fgetc(fp)) == SEP); // printf("%c", ch);
	}
	
	//print("%d", matrix[idx]);
      } // cols l
      
      // search for end of line
      if (ch != EOL && ch != EOF) while ((ch=fgetc(fp)) != EOL && ch!=EOF);
    } // hr
   
    UPDATE;
    
    // printf("\n");
  } // rows r
  
  CLOSE
}



SEXP CreateCode23(int individuals, int snps, int mem, int blocks,
		  int blocklength, int tablesize, SEXP file) {
  if (Information == R_NilValue) Init();
  SEXP Code =  CreateEmptyCodeVector(snps, individuals, mem, REALSXP,
				     Information, true),
    Infos = getAttrib(Code, Information);
  Uint *pi  = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  if (file != R_NilValue) {
    SEXP fileInfos = getAttrib(file, Information);
    Uint *fileinfo = (Uint *) INTEGER(VECTOR_ELT(fileInfos, INFO_INFO));
    pi[SNPxIND] = fileinfo[SNPxIND];
  }
  pi[INFO_BLOCKS] = blocks;
 
 if (PL > 1) {
    Uint
      sizeM = (Uint) (sizeof(uint64) * mem / 1048576);  
    PRINTF("Data: %d individuals and %d SNPs\nStorage mode: %d block(s) of length %d\nSize of M: %d MB\nSize of table: %d kB\n", 
	   individuals, snps,
	   blocks, blocklength, 
	   sizeM < 1 ? 1 : sizeM,
	   (Uint) (tablesize / 1024));
  }

 return(Code);
}

  

SEXP file_coding(SEXP file) {
  bool twobit = GLOBAL.relationship.method == TwoBit;
  assert(twobit || GLOBAL.relationship.method == ThreeBit);
  return file_coding_intern(file, 
			    twobit ? matrix_coding_start2 :matrix_coding_start3,
			    twobit ? matrix_coding2 : matrix_coding3,
			    file);
}


SEXP failure(SEXP VARIABLE_IS_NOT_USED file) {
  ERR("method is unknown or does not work with SNP matrices given by file.");
  return R_NilValue;
}


SEXP matrix_coding(SEXP SNPxIndiv){
  Uint
    method = GLOBAL.relationship.method,
    snps = nrows(SNPxIndiv),
    individuals = ncols(SNPxIndiv);
  SEXP Ans = (method == TwoBit
	      ? matrix_coding_start2(individuals, snps, R_NilValue)
	      : matrix_coding_start3(individuals, snps, R_NilValue));

  //  return(Ans);
  
  SEXP
    Infos = getAttrib(Ans, Information);
  Uint
    //    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
    *M = (Uint*) INTEGER(SNPxIndiv);

  if (method == TwoBit) {
    matrix_coding2(M, 0, individuals, 0, snps, Ans, NULL);
  } else if (method == ThreeBit) {
    matrix_coding3(M, 0, individuals, 0, snps, Ans, NULL);
  } else BUG;
 
#define FORPS(A)					\
  for (Uint s=0; s<snps; s++) {				\
    /* //printf("ps = %10g\n" , p[s]);	*/		\
    double ps = p[s] = A,				\
      p2 = ps * ps;					\
    sump += ps;						\
    ppt += p2;						\
    sum_pi_qi += ps - 0.5 * p2;				\
}
  
  //  Uint n_compressed = info[MEM] / individuals;
  double ppt, sum_pi_qi, sump, 
    *P = REAL(VECTOR_ELT(Infos, INFO_P)),
    *p = P + INFO_P_P,
    inv_individuals = 1.0 / (double) individuals,
    digits = GLOBAL.relationship.digits;

  ppt = sum_pi_qi = sump = 0.0;
  if (digits >= 0) {
    double factor = POW(10, digits) * 0.5, // 27.11.2014 " * 0.5" 
      factor_inv = factor * inv_individuals;
    FORPS(ROUND(p[s] * factor_inv) / factor)
  } else { FORPS(p[s] * inv_individuals) }
  P[INFO_P_PPT] = ppt;
  P[INFO_P_SUMPQ] = sum_pi_qi;
  P[INFO_P_SUMP] = sump * individuals;

 
  return Ans;
}


SEXP do_centering_etc(SEXP SNPxIndiv, SEXP Ans) {
  SEXP
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
    individuals = info[INDIVIDUALS];
  bool centered = GLOBAL.relationship.centered == True,
    normalized = GLOBAL.relationship.centered == True;

  DoCentering(REAL(Ans), individuals, centered, normalized,
	      REAL(VECTOR_ELT(Infos, INFO_P))[INFO_P_SUMP]);
  return Ans;
}


SEXP getRelmatrixIndividuals(SEXP Z, SEXP SNPxIndiv) {
  SEXP Ans,
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
    *zz = (Uint*) INTEGER(Z),
    individuals = info[INDIVIDUALS];
  PROTECT(Ans=allocMatrix(REALSXP, individuals, individuals));
  double
    *ans = REAL(Ans);
  
  for (Uint c=0; c<individuals; c++) {
    for (Uint r=c; r<individuals; r++) {
      ans[r * individuals + c] = ans[r + individuals * c] =(double) *(zz++);
    }
  }

  UNPROTECT(1);
  return do_centering_etc(SNPxIndiv, Ans);
}


SEXP getRelmatrixSNP(SEXP Z, SEXP SNPxIndiv) {
  if (TYPEOF(Z) != VECSXP) ERR("unknown type in 'getmatrixSNP'");
  SEXP Ans,
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
    individuals = info[INDIVIDUALS],
    isq = individuals * individuals,
    len = length(Z);
  PROTECT(Ans=allocMatrix(REALSXP, individuals, individuals));
  double
    *ans = REAL(Ans);
   
  for (Uint i=0; i<isq; ans[i++] = 0.0);
  for (Uint i=0; i<len; i++) {
    double *e = REAL(VECTOR_ELT(Z, i));
    for (Uint r=0; r<isq; r++) ans[r] += e[r];
  }
  UNPROTECT(1);
  return do_centering_etc(SNPxIndiv, Ans);
}


void printbits(uint64 x, Uint size, Uint bits) {
  uint64 mask = INT64_C(1);  
  for (Uint d=0, zaehler=0; d<size; d++) {
    if (zaehler++ == bits) {
      PRINTF(".");
      zaehler = 1;
    }
    PRINTF("%d", (x & mask) != 0);
    x >>= 1;
  }
  PRINTF(" ");
}
//erwandtschaft.cc:657:74: error: invalid conversion from ‘void (*)(unsigned int*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, uint64_t*, unsigned int*, bool) {aka void (*)(unsigned int*, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, long unsigned int*, unsigned int*, bool)}’ to  ‘scalarproduct {aka void (*)(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, long unsigned int*, unsigned int*, bool)}’ [-fpermissive]
 


typedef void (*scalarproduct)(Uint individuals, 
			      Uint min_row, Uint max_row,
			      Uint min_col, Uint max_col,
			      Uint compressed_start, Uint compressed_end,
			      Uint n_compressed,
			      uint64 *M, double *ergb, bool Parallel);

SEXP MmPsqSNPs(SEXP SNPxIndiv, SEXP Compressed_start, SEXP Compressed_end) {
  SEXP
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
    individuals = info[INDIVIDUALS],
    n_compressed = info[MEM] / individuals,
    compressed_start = INTEGER(Compressed_start)[0],
    compressed_end = INTEGER(Compressed_end)[0]
    ;
  uint64 *M = (uint64*) REAL(SNPxIndiv);
  scalarproduct scalar = INTEGER(getAttrib(SNPxIndiv, Method))[0] == TwoBit
    ? MmPsq2 : MmPsq3;
  SEXP Ans;
  PROTECT(Ans=allocVector(REALSXP, individuals * individuals));
  double *ans = REAL(Ans);
    
  if (true) {
    scalar(individuals, 0, individuals, 0, individuals, compressed_start,
	   compressed_end, n_compressed, (uint64*) REAL(SNPxIndiv), ans, false);
  } else {
    Uint
      i = 0,
      delta = compressed_end - compressed_start + 1,
      bytes = individuals * delta * sizeof(uint64);
    uint64 *Mpt = NULL,
      *MX = NULL;
    
    if ((MX = (uint64*) MALLOC(bytes)) == NULL)
      ERR("memory allocation error!");
    for (Uint c=0; c<individuals; c++) {
      Mpt = M + c * n_compressed;
      for (Uint r=compressed_start-1; r< compressed_end; ) {
	MX[i++] = Mpt[r++];
      }
    }
    scalar(individuals, 0, individuals, 0, individuals, 0, delta, delta, MX,
	   ans, false);
    FREE(MX);
  }

  UNPROTECT(1);
  return Ans;
}


SEXP MmPsqIndividualsParallel(SEXP SNPxIndiv) {  
  SEXP
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  Uint 
    individuals = info[INDIVIDUALS],
    n_compressed = info[MEM] / individuals;
  scalarproduct scalar = INTEGER(getAttrib(SNPxIndiv, Method))[0] == TwoBit
				 ? MmPsq2 : MmPsq3;

  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, individuals, individuals));

  scalar(individuals, 0, individuals, 0, individuals,
	 0, n_compressed, n_compressed,
	 (uint64*) REAL(SNPxIndiv), REAL(Ans), true);

  UNPROTECT(1);
  return do_centering_etc(SNPxIndiv, Ans);
}




SEXP get_matrix_start(Uint individuals, Uint snps,
		      SEXP VARIABLE_IS_NOT_USED G) {
  SEXP Ans;
  PROTECT(Ans=allocMatrix(INTSXP, snps, individuals));
  //if (PL > 1) PRINTF("Data: %d individuals and %d SNPs\n", individuals, snps);
  UNPROTECT(1);
  return Ans;
}

SEXP XXXget_matrix(SEXP SNPxIndiv) {
  typedef double (*get_type)(uint64 *Indiv, Uint s);
  SEXP
    Infos = getAttrib(SNPxIndiv, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),  
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    n_compressed = info[MEM] / individuals,
    endfor = individuals * n_compressed;
  uint64 *M = (uint64*) REAL(SNPxIndiv);
  SEXP Ans = get_matrix_start(individuals, snps, R_NilValue);
  double *ans = REAL(Ans);

  //printf("n_compressed = %d\n", n_compressed); BUG;

  get_type value = INTEGER(getAttrib(SNPxIndiv, Method))[0] == TwoBit
    ? getValue2: getValue3;
  for (Uint a=0; a<endfor; a+=n_compressed) {    
    for (Uint s=0; s<snps; s++, ans++) *ans = value(M + a, s);
  }
  return Ans;
}

 
SEXP matrix_get(SEXP SNPxIndiv) {
  if (Information == R_NilValue) Init();
  return INTEGER(getAttrib(SNPxIndiv, Method))[0] == TwoBit
    ? get_matrix2(SNPxIndiv) : get_matrix3(SNPxIndiv);
}



void get_file(Uint *M, Uint start_individual, Uint end_individual, 
		Uint start_snp, Uint end_snp, SEXP Ans,
		double VARIABLE_IS_NOT_USED *G) {
  SEXP
    Infos = getAttrib(Ans, Information);
  Uint
    *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO));
  Uint snps = info[SNPS],
    *ans = (Uint *) INTEGER(Ans);
  for (Uint a=start_individual; a<end_individual; a++) {
    Uint *Mptr = ans + a * snps;
    for (Uint s=start_snp; s<end_snp; M++) {
      Mptr[s++] = *M;
    }
  }
}

SEXP file_get(SEXP file) {
  SEXP Ans = file_coding_intern(file, get_matrix_start,	get_file,
				R_NilValue);
  return Ans;
}



SEXP file_dot_start(Uint individuals, Uint snps, SEXP G) {
  SEXP MTdot;
  if ((Uint) length(G) !=  individuals)
    ERR("vector must have length equal to number of individuals");
  if (PL > 1) {PRINTF("Data: %d individuals and %d SNPs\n", individuals, snps);}
  PROTECT(MTdot=allocVector(REALSXP, snps));
  for (Uint i=0; i<snps; REAL(MTdot)[i++] = 0.0);
  UNPROTECT(1);
  return MTdot;
}

void file_dot_do(Uint *M, Uint start_individual, Uint end_individual, 
		Uint start_snp, Uint end_snp, SEXP Ans, double *G) {
  double *Fdot = REAL(Ans);
  for (Uint a=start_individual; a<end_individual; a++) {
    double dummy = G[a];
    for (Uint s=start_snp; s<end_snp; M++) {
      Fdot[s++] += dummy * (double) *M;
    }
  }
}

SEXP file_dot(SEXP file, SEXP G) {
  return file_coding_intern(file, file_dot_start, file_dot_do, G);
  //  SEXP  Infos = getAttrib(SNPxIndiv, Information);
  //  Uint *info = (Uint *) INTEGER(VECTOR_ELT(Infos, INFO_INFO)),
  
  /*
  usrBoolean centered = GLOBAL.relationship.centered;
  if (centered != False) {
    double *mtdot = REAL(MTdot),
      *centered = GLOBAL.relationship.pcentered,
      g = REAL(G),
      sumg = 0.0;
    Uint individuals = length(G),
      snps = length(MTdot),
      len = GLOBAL.relationship.ncentered;
    if (snps != len) ERR("length of 'centered' must equal the number of SNPs.");
    if (GLOBAL.relationship.centered == Nan) 
      for (Uint i=0; i<individuals; i++) sumg += centered[i] * g[i];
    else {
      double *p = 
    }
      
    for (Uint i=0; i<snps; i++) mtdot[i] -= sumg;
  }
  */
  
  // return Ans;
}





SEXP  dot_file_start(Uint individuals, Uint snps, SEXP G) {
  SEXP Ans;
  if ((Uint) length(G) != snps) ERR("vector must have length equal to number of snps");
  if (PL > 1) {PRINTF("Data: %d individuals and %d SNPs\n", individuals, snps);}
  PROTECT(Ans=allocVector(REALSXP, individuals));
  for (Uint i=0; i<individuals; REAL(Ans)[i++] = 0.0);
  UNPROTECT(1);
  return Ans;
}


void dot_file_do(Uint *M, Uint start_individual, Uint end_individual, 
		  Uint start_snp, Uint end_snp, SEXP Ans, double *G) {
  double *dotF = REAL(Ans);
  for (Uint a=start_individual; a<end_individual; a++) {
    for (Uint s=start_snp; s<end_snp; M++) {
      dotF[a] += G[s++] * (double) *M;
    }
  }
}


SEXP dot_file(SEXP file, SEXP G) {
  return file_coding_intern(file, dot_file_start, dot_file_do, G);
}





SEXP substract_centered(SEXP SnpXindiv) {
  Uint indiv = ncols(SnpXindiv),
    snps = nrows(SnpXindiv),
    len = GLOBAL.relationship.ncentered;
  double 
    *centered = GLOBAL.relationship.pcentered;
  assert(centered != NULL && GLOBAL.relationship.centered==Nan);
  if (snps != len) ERR("length of 'centered' must equal the number of SNPs.");
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, snps, indiv));
  double *ans = REAL(Ans),
    *snpXindiv = REAL(SnpXindiv);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
  for (Uint i=0; i<indiv; i++) {
    int i_snps = i * snps;
    double *s = snpXindiv + i_snps,
      *a = ans + i_snps;
    for (Uint j=0; j<snps; j++) a[j] = s[j] - centered[j];
  }
  UNPROTECT(1);
  return Ans;
}


SEXP get_centered() {
  Uint len = GLOBAL.relationship.ncentered ;
  double *centered = GLOBAL.relationship.pcentered;
  assert (centered != NULL && GLOBAL.relationship.centered==Nan);
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, len));
  double *ans = REAL(Ans);
  MEMCOPY(ans, centered, sizeof(double) * len);
  UNPROTECT(1);
  return Ans;
}


