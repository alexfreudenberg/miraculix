
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

// _MM_ALIGN16 Uint ccc;

#include "AutoMiraculix.h"
#include <General_utils.h>
#include "xport_import.h"
#include "Haplo.h"
#include "Bit23.h"
#include "shuffle.h"
#include "packed.h"
#include "multiply.h"
#include "hamming.h"
#include "nocoding.h"
#include "options.h"
#include "miraculix.h"
#include "align.h"
#include "mmagpu.h"
#include "utils.h"

     

snpcoding check_method(snpcoding method);

typedef SEXP (*coding_start_t)(Uint , Uint, SEXP);
typedef void (*coding_main_t)(Uint *, Uint, Uint, Uint, Uint, Uint,
			      SEXP, double *);
typedef Uint (*codes_per_block_t)();
typedef void (*crossprod_t)(Uint*, Uint, Uint, double *);
typedef void (*h2g_t)(Uint *, Uint, Uint, Uint, Uint *);

char EOL = '\n';

void allInfo(Uint *info) { // ok
  for(Uint i=0; i<=INFO_GENUINELY_LAST; i++) {
    if (i==ADDR1 || i==SUMGENO_E9 || i==ALIGNADDR1 || i==ALIGNEDUNITS1 ||
	i==MEMinUNITS1) continue;
    PRINTF("%s=", INFO_NAMES[i]);
    switch(i) {
    case METHOD : PRINTF("%s", SNPCODING_NAMES[info[i]]); break;
    case ADDR0 : PRINTF(" %u %u", info[ADDR1], info[ADDR0]); break;
    case ALIGNADDR0 :
      PRINTF("%u %u", info[ALIGNADDR1], info[ALIGNADDR0]); break;
    case SUMGENO : PRINTF("%lu", SumGeno(info)); break;
    case MEMinUNITS0 : PRINTF("%lu", MemInUnits(info)); break;
    case ALIGNEDUNITS0 : PRINTF("%lu", AlignedInUnits(info)); break;
    default:
      if ((int) info[i] == NA_INTEGER) PRINTF("NA"); else PRINTF("%u", info[i]);
    }
    PRINTF("\n");
  }
}
void allInfo(SEXP M) { // ok
  Uint *info = GetInfoUnchecked(M);
  allInfo(info); // ok
}
 

Uint GetBytesPerBlock(snpcoding method) {
  switch (method) {
  case Shuffle256 :return BytesPerBlockShuffle256();
  case Shuffle :return BytesPerBlockShuffle();
  case TwoBit : return BytesPerBlock2();
  case Packed :return BytesPerBlockPacked();
  case Packed256: case MMAGPU :return BytesPerBlockPacked256();
  case Multiply :return BytesPerBlockMultiply();
  case Multiply256 :return BytesPerBlockMultiply256();
  case ThreeBit : return BytesPerBlock3();    
  case Hamming2 : case Hamming3 : return BytesPerBlockH();
  case NoSNPcoding: case NoSNPcodingR: return BytesPerBlockPlain();
  case AutoCoding : BUG;
  case Haplo : return BytesPerBlockHaplo();
  default : BUG;
  }
}


Uint GetCodesPerBlock(snpcoding method) {
  switch (method) {
  case Shuffle256 :return CodesPerBlockShuffle256();
  case Shuffle :return CodesPerBlockShuffle();
  case TwoBit : return CodesPerBlock2();
  case Packed :return CodesPerBlockPacked();
  case Packed256: case MMAGPU : return CodesPerBlockPacked256();
  case Multiply :return CodesPerBlockMultiply();
  case Multiply256 :return CodesPerBlockMultiply256();
  case ThreeBit : return CodesPerBlock3();    
  case Hamming2 : case Hamming3 : return CodesPerBlockH();
  case NoSNPcoding: case NoSNPcodingR: return CodesPerBlockPlain();
  case AutoCoding : BUG;
  case Haplo : return CodesPerBlockHaplo();
  case UnknownSNPcoding : return 0;
  default : BUG;
  }
}

Uint GetBitsPerCode(snpcoding method) {
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed :  case Packed256 :
  case Haplo : case Multiply : case Multiply256 : case MMAGPU :
    return BitsPerCode2();
  case ThreeBit : return BitsPerCode3();    
  case Hamming2 : case Hamming3 : return BitsPerCodeH();
  case NoSNPcoding: case NoSNPcodingR: return BitsPerCodePlain();
  case AutoCoding : BUG;
  default : BUG;
  }
}

Uint GetUPI(Uint snps, snpcoding method) {
  switch(method) {
  case Shuffle256 : case Shuffle :case TwoBit : case Packed: case Packed256:
  case Haplo : case Multiply :case Multiply256 : case MMAGPU :
    return UnitsPerIndiv256(snps);
  case ThreeBit : return UnitsPerIndiv3(snps); 
  case Hamming2 : return UnitsPerIndivH2(snps);
  case Hamming3 : return UnitsPerIndivH3(snps); 
  case NoSNPcoding : return UnitsPerIndivPlain(snps); 
  default : BUG;
  }
  
  return 0;
}




SEXP get_centered() {
  option_type *global = &(KEYT()->global);
  Uint len = global->genetics.ncentered ;
  double *centered = global->genetics.pcentered;
  assert (centered != NULL && global->genetics.centered==Nan);
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, len));
  double *ans = REAL(Ans);
  MEMCOPY(ans, centered, sizeof(double) * len);
  UNPROTECT(1);
  return Ans;
}
			  

void DoCentering(double *ans, Uint individuals, bool centred, bool normalized,
		 Ulong sumGeno, Ulong snps) {

  // #define maxlong  9223372036854775807
#define maxlong       9223372036854775000.0 // for safty
    
  double
    nd = (double) individuals,
    nSq = nd * nd,
    factor = RF_NAN,
    sum_pi_qi = RF_NAN;
  if (centred || normalized) {
    double
      sum_P = RF_NAN,
      *sums = (double *) MALLOC(sizeof(double) * individuals),
      totalsum = 0.0;

    if (2.0 * 4.0 * nd * nd * (double) snps > maxlong)
      ERR("");
    
    Ulong k = 0,
      nP1 = individuals + 1;
    for (Uint i=0; i<individuals; i++) { // ok
      double dummy = 0.0;
      for (Uint j=0; j<individuals; j++) dummy += ans[k++];// ok
      sums[i] = nd * dummy;
      totalsum += dummy;
    }

    if (centred) {
      for (Ulong i=0; i<individuals; i++) {
	Ulong idx = i * nP1;
	for (Ulong j=i; j<individuals; j++, idx++) { 
	  ans[idx] = nSq * ans[idx] - sums[i] - sums[j] + totalsum;
	}
      }
    }
    
    if (normalized) {
      sum_P = nd *  (double) sumGeno; // ^= sigma^2 * nd^2
      factor = sum_pi_qi = sum_P - 0.5 * totalsum;
      if (!centred) factor /= nSq;
    } else factor = nSq;
    
    if (factor <= 0.0) ERR("strange input matrix?")
    
    for (Ulong i=0; i<individuals; i++) {
      Ulong idx = i * nP1;
      for (Ulong j=i; j<individuals; j++) {
	ans[i + j * individuals] = (ans[idx++] /= factor);
      }
    }
    
    FREE(sums);
  }  
  // return sum_pi_qi /nSq ;
}


#define UPDATE								\
  if (isSNPxInd) {							\
    Uint jj, rP1 = r + 1;						\
    if (rP1 % nrow_matrix == 0) {					\
      /* always delivering matrix as snps x individuals */		\
      coding_main(matrix, 0, individuals, rP1-codesperblock, rP1,	\
		  nrow_matrix,  Ans, dG);				\
      rowidx = -plusrow;						\
      for (jj=0; jj<matrix_size; matrix[jj++] = 0);			\
    }									\
  } else {								\
    Uint jj;								\
    coding_main(matrix, r, r + 1, 0, snps, nrow_matrix, Ans, dG);	\
    rowidx = -plusrow;							\
    for (jj=0; jj<matrix_size; matrix[jj++] = 0);			\
  }


#define CLOSE								\
  if (isSNPxInd) {							\
    Uint remainder = r % nrow_matrix;					\
    if (remainder != 0) {						\
      coding_main(matrix, 0, individuals, (Uint)(r - remainder), r,	\
		  nrow_matrix,  Ans, dG);				\
    }									\
  } /* KEIN else  */							\
									\
  fclose(fp);								\
  FREE(matrix); 
  

SEXP file_binary_intern(char *file, 
			char *coding, Uint header, bool isSNPxInd, 
			bool doubledindividuals,
			// output
			Uint snps,
			Uint individuals,
			Uint codesperblock,
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
    SPRINTF(msg, "matrix is supposed to be a '%.20s' x '%.20s' matrix with %d individuals and %d SNPs, which looks odd.\n", 
	    isSNPxInd ? "SNPs" : "individuals",
	    isSNPxInd ? "individuals" : "SNPs",
	    individuals, snps);
    if (snps > 0 && individuals > 0) { warn(msg) } else ERR(msg);
  }
  
  if (isSNPxInd) {
    plusrow = 1;
    nrow_matrix = codesperblock; 
    rows = snps;
    cols = ncol_matrix = individuals;
  } else {
    // plusrow = snps;
    plusrow = codesperblock;
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

  SEXP Ans = coding_start(snps, individuals, G);
  Uint matrix_size = nrow_matrix * (ncol_matrix + 3);
  Ulong sumgeno = 0;
  double *dG = NULL;
  if (length(G) > 0) dG = REAL(G);

  if ((matrix = (Uint*) CALLOC(matrix_size, sizeof(Uint))) == NULL)
    ERR("memory space could not be acquired");
  cols = (Uint) CEIL(0.25 * cols);
  perbyte = 4 / haplocol;
  char twobitpattern = 3;


  if ((fp = fopen(file, "r")) == NULL)  {
    ERR1("file '%.50s' could not be opened", file);
  }
  for (i=0; i<header; i++) fgetc(fp);

  for (rowidx=r=0; r<rows; r++, rowidx+=plusrow) { // number of (multiple) rows 
    for (hr=0; hr<haplorow; hr++) { 
      // genuine loop only if haplo types are given rowwise
      idx = rowidx;

      for (l=0; l<cols; l++) { // number of (multiple cols 
	Uint k = 0;
	ch = fgetc(fp);
	for (j = 0; j < perbyte; j++, idx+=nrow_matrix) {
	  for (hc=0; hc<haplocol; hc++, k+=2) {
	    char value = (ch >> k) & twobitpattern;
	    //	  for (hc=0; hc<haplocol; hc++, ch >>= 2) {
	    // char value = ch & twobitpattern;
	    // genuine loop only if haplo types are given column wise

	    if (value != 0) { // == 0
	      assert(idx <= matrix_size);
	      if (value == 2) {
		matrix[idx]++;
		sumgeno++;
	      } else if (value == 3) {
		matrix[idx] += 2;
		sumgeno += 2;
	      } else if (value == 1) ERR("missing value detected.")
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
   
  CLOSE;

  Uint *info = GetInfoUnchecked(Ans);						
  StoreSumGeno(sumgeno);			       

  return Ans;
}




SEXP file_intern(SEXP file,
		 coding_start_t coding_start,
		 coding_main_t coding_main,
		 Uint codesperblock,
		 SEXP G) {

  Uint *info = GetInfoUnchecked(file),
    leadingcols = info[LEADINGCOL],
    isSNPxInd = info[SNPxIND],
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    doubledindividuals = info[DOUBLEINDIV];
  Rint header = info[HEADER]; // can be indeed negative !!
  info[CODESPERBLOCK] = codesperblock;

  char *coding = (char*) CHAR(STRING_ELT(getAttrib(file, Coding), 0));
  char *filename = (char*) CHAR(STRING_ELT(file, 0));

  if (header < 0) {
    // only for plinkbinary: snps & individuals already known
    return file_binary_intern(filename, coding,
			      -header, isSNPxInd, doubledindividuals,
			      snps, individuals,
			      codesperblock, 
		       // output
			      coding_start, coding_main, G
			      );
  } // else snps not known yet:

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

  FILE *fp;
  Uint nrow_matrix, ncol_matrix,
    *matrix = NULL, // snps x individuals
    rows, // rowsM1, 
    cols, colsM1,
    plusrow,  haplorow, haplocol;
  Ulong sumgeno = 0;
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
    SPRINTF(msg, "matrix is supposed to be a '%.20s' x '%.20s' matrix with %d individuals and %d SNPs, which looks odd.\n", 
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
    nrow_matrix = codesperblock; 
    ncol_matrix = individuals;
  } else {
    // plusrow = snps;
    plusrow = codesperblock;
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
  SEXP Ans = coding_start( snps, individuals,G);
  Uint matrix_size = nrow_matrix * ncol_matrix;
  if ((matrix = (Uint*) CALLOC(matrix_size, sizeof(Uint))) == NULL)
    ERR("memory space could not be acquired");
   double *dG = NULL;
  if (length(G) > 0) dG = REAL(G);

  rows /= haplorow;
  cols /= haplocol;

 
  fp = fopen(filename, "r");

  // jump header lines
  for (int i=0; i<header; i++) while ((ch=fgetc(fp)) != EOL);

  Uint r=0;
  for (Uint rowidx=0; r<rows; r++, rowidx+=plusrow) {//# of (multiple) rows
    for (Uint hr=0; hr<haplorow; hr++) { 
      // genuine loop only if haplo types are given rowwise
      Ulong idx = rowidx;
      
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
	    if (ch == B) {
	      matrix[idx]++;
	      sumgeno++;
	    } else if (ch == C) {
	      matrix[idx] += 2;
	      sumgeno +=2;
	    } else if (ch == NA) ERR("missing value detected.")
	    else {
	      PRINTF(">%c< row=%d of %d rows, col=%d of %d columns; plus=%d hc=%d haplocol=%d\n", ch, r, rows, l, cols, nrow_matrix, hc, haplocol);
	      ERR("Unknown character detected. Note that missings are not allowed here.");
	    }
	    //assert(matrix[idx] <=2);
	  } 
	  if (l<colsM1 || hc!=haplocol-1)
	    while((ch = fgetc(fp)) == SEP); 
	}
      } // cols l
      
      // search for end of line
      if (ch != EOL && ch != EOF) while ((ch=fgetc(fp)) != EOL && ch!=EOF);
    } // hr
   
    UPDATE;
    
  } // rows r
  
  CLOSE;

  StoreSumGeno(sumgeno);			       

  return Ans;
}




Ulong sumGeno(Uint * SNPxIndiv, Uint snps, Uint individuals, snpcoding method) {
  switch (method) {
  case Shuffle :  case Shuffle256 : case TwoBit : case Packed: case Packed256:
  case Multiply : case Multiply256 : case MMAGPU :
    if (use256()) return sumGeno256(SNPxIndiv, snps, individuals);
    if (use128()) return sumGeno128(SNPxIndiv, snps, individuals); 
    return sumGeno2(SNPxIndiv, snps, individuals);
  case ThreeBit : return sumGeno3(SNPxIndiv, snps, individuals);
  case Hamming2 : return sumGenoH2(SNPxIndiv, snps, individuals);
  case Hamming3 : return sumGenoH3(SNPxIndiv, snps, individuals);
  case NoSNPcoding:
  case NoSNPcodingR: return sumGenoPlain(SNPxIndiv, snps, individuals);
  case AutoCoding :
  case Haplo : BUG;
  default : BUG;
  }

  return 0L;
}


void showInfo(Uint *info) {
  PRINTF("version=%d\n", info[VERSION]);
  PRINTF("snps=%d\n", SNPS);
  PRINTF("indiv=%d\n", INDIVIDUALS);
  PRINTF("Addr=%d %d\n", ADDR0, ADDR1);
  PRINTF("Align=%d %d\n", ALIGNADDR0 , ALIGNADDR1);
  PRINTF("sumgeno=%d %d\n", SUMGENO, SUMGENO);
  PRINTF("method=%d\n", METHOD);
  PRINTF("algn=%d\n", ALIGNMENT);
  PRINTF("snpXind=%d\n", SNPxIND);
  PRINTF("BPC=%d\n", BITSPERCODE);
  PRINTF("BPB=%d\n", BYTESPERBLOCK);
  PRINTF("CPB=%d\n", CODESPERBLOCK);
  PRINTF("header=%d\n", HEADER);
  PRINTF("doubleindiv=%d\n", DOUBLEINDIV);
  PRINTF("leading=%d\n", LEADINGCOL);
  PRINTF("mem=%d %d\n", MEMinUNITS0, MEMinUNITS1);
  PRINTF("allgnedunits=%d %d\n", ALIGNEDUNITS0, ALIGNEDUNITS1);
  PRINTF("UPI=%d\n", UNITSPERINDIV);
}

void showInfoShort(Uint *info) {
  PRINTF("V=%d, ", info[VERSION]);
  PRINTF("S/I=%d %d, ", SNPS, INDIVIDUALS);
  PRINTF("Ad=%d %d, ", ADDR0, ADDR1);
  PRINTF("Al=%d %d, ", ALIGNADDR0 , ALIGNADDR1);
  PRINTF("M=%d, ", METHOD);
  PRINTF("a=%d, ", ALIGNMENT);
  PRINTF("BPC=%d %d %d %d\n",
	 BITSPERCODE, BYTESPERBLOCK, CODESPERBLOCK, UNITSPERINDIV);
}

Uint *AlignBase(SEXP CM, Uint VARIABLE_IS_NOT_USED nr, Uint alignment, bool test) {
  Uint *info;
  if (test) info = GetInfo(CM); else info = GetInfoUnchecked(CM);
  // showInfoShort(info);
  //

  if (info[ALIGNMENT] != alignment) {
    ERR2("currently, data exchange between different kinds of machines (e.g. AVX2/SSE) is not possible. (Got the aligments %d and %d Bytes.)", info[ALIGNMENT], alignment);
  }
 
  Uint *address = (Uint *) INTEGER(CM),
    *algnaddress = (Uint *) algn_general(INTEGER(CM), alignment),
    *infoaddress = (Uint*) GetAddress(info, ADDR0);
  
  
#ifdef SCHLATHERS_MACHINE
  //  if (info[INDIVIDUALS] != 10) {
  ///int units = info[ UNITSPERINDIV ] * info[INDIVIDUALS];
  //    printf("zaehler = %d %d %d %d %d +++++",info[ZAEHLER], units, info[SNPS], info[INDIVIDUALS], alignment);
  //    for (int i=0; i<units; i++) printf("%d ", algnaddress[i]); printf("\n");
  //  }
#endif    
  
  option_type *global = &(KEYT()->global);
  bool aligned=(uintptr_t)infoaddress % alignment == (uintptr_t) address % alignment;
  if (global->messages.warn_address && infoaddress != address) {
    Uint *recentaddress = (Uint*) GetAddress(info, RECENTALIGNADDR0);    
    if (recentaddress != address) {
      PRINTF("Address has changed in a coded object%s (*%ld->*%ld). %s\n",
	     aligned ?  ", but same modulus" : " by 'gc' or disk saving",
	     (Uint) (uintptr_t) infoaddress, (Uint) (uintptr_t) address,
	     recentaddress == NULL
	     ? "[This messages can be suppressed by 'RFoptions(warn_address=FALSE)']"// OK
	     : ""); 
      ADDADDRESS((Rint*) address, RECENTALIGNADDR0);
    }
  }
  
  if (aligned){
    assert((uintptr_t) algnaddress % alignment == 0);
    return algnaddress;
  }

  Uint bytes = MemInUnits(info) * BytesPerUnit;
#ifdef DO_PARALLEL
  // prevent two processes to move data at the same time
  if (CORES > 1) {
    int mypid, sleep = 1;
    Ext_pid(&mypid);
    if (!info[BLOCKEDINFO]) {
      info[BLOCKEDINFO] = mypid;
      Ext_sleepMicro(&sleep); // wait. Other processes may access simulatenously
    }
    if ((Uint) mypid != info[BLOCKEDINFO]) { // other process has been
      // doing the job, maybe simultaneously
      sleep = bytes / 1000 + mypid; // some randomness in sleeping
      for (int i=0; i<10; i++) {
	if (!info[BLOCKEDINFO]) return AlignBase(CM, nr, alignment, test);
 	Ext_sleepMicro(&sleep);
      }
      if (OPTIONS_UTILS->basic.warn_parallel)
	warn("Simultaneous write access: there is some risk of loosing data. You can suppress this message by 'RFoptions(warn_parallel_write=FALSE)', or avoid any risk by 'RFoptions(cores=1)'."); // ok
      info[BLOCKEDINFO] = 0; 
      return AlignBase(CM, nr, alignment, test);
    }
  }
#endif
  
  Uint *oldInNew =
    address + ((Uint*) GetAddress(info, ALIGNADDR0) - infoaddress);
  MEMMOVE(algnaddress, oldInNew, bytes);
  ADDADDRESS((Rint*) address, ADDR0);
  ADDADDRESS((Rint*) algnaddress, ALIGNADDR0);
#ifdef DO_PARALLEL
  info[BLOCKEDINFO] = 0;
#endif

  return algnaddress;
}

  
Uint* DoAlign(SEXP CM, Uint nr, snpcoding method, bool test) {
  switch(method) {
  case Shuffle : case Shuffle256 : case TwoBit : case Packed: case Packed256:  
  case Haplo : case Multiply : case Multiply256 : case MMAGPU :
    return Align256test(CM, nr, test);
  case ThreeBit : return Align3(CM, nr, test);
  case Hamming2 :
  case Hamming3 : return AlignH(CM, nr, test);
  case NoSNPcoding:
  case NoSNPcodingR:  return AlignPlain(CM, nr, test);
  case AutoCoding : BUG;     
  default : ERR("Method is unknown!");
  }
  return NULL;
}

//Uint* DoAlignWithoutTest(SEXP CM, Uint nr, snpcoding method) {
//  // very dangerous command !!!
//  DoAlign(CM, nr, method, false);
//}

Uint* DoAlign(SEXP CM, Uint nr, snpcoding method) {
  return DoAlign(CM, nr, method, true);
}

void crossprod(Uint * SNPxIndiv, Uint snps, Uint individuals,
		 snpcoding method,  bool centred, bool normalized,
		 Ulong SumGeno, double *A) {
  assert(individuals > 0);
  double nSq = (double) individuals * (double) individuals;      
  if ((double) snps * nSq * 4.0 > 4.5035996e+15) // 2^{manissen-bits = 52}
    ERR("matrix too large to calculate the relationship matrix -- pls contact maintainer of 'miraculix'");
  crossprod_t cp = NULL;

  //  printf("method crsosprod = %d\n", method);
  
  switch (method) {
  case Shuffle : case Shuffle256 : case TwoBit : case Packed :  
    if (useShuffle256(method)) cp = crossprod_shuffle256; 
      else if (useShuffle(method)) cp = crossprod_shuffle;
      else if (usePacked(method)) cp = crossprod_packed;
      else cp = crossprod2;
    break;
  case Packed256 :  cp = crossprod_packed256; break;
  case MMAGPU : cp = crossprod_mmagpu; break;
  case Multiply :  cp = crossprod_multiply; break;
  case Multiply256 :  cp = crossprod_multiply256; break;
  case ThreeBit : cp = crossprod3; break;    
  case Hamming2 : cp = crossprod_H2; break;   
  case Hamming3 : cp = crossprod_H3; break;    
  case NoSNPcoding:
  case NoSNPcodingR: cp = crossprod_Plain; break;    
  case AutoCoding :
  case Haplo : BUG;      
  default : BUG;
  }
  cp(SNPxIndiv, snps, individuals, A); 

  assert(A[0] >= 0.0);
  if (centred || normalized) {
    if (SumGeno == 0) // falls die Matrix tatsaechlich nur aus 0-en besteht,
      // wird eben doppelt gerechnet. Praktisch irrelevanter Fall
      SumGeno = sumGeno(SNPxIndiv, snps, individuals, method);
    
    assert(individuals > 0);
    DoCentering(A, individuals, centred, normalized, SumGeno, snps);
    //  assert(A[0] >= 0.0);
  }
  
}


double *crossprod(Uint * SNPxIndiv, Uint snps, Uint individuals,
		    snpcoding method,  bool centred, bool normalized,
		    Ulong SumGeno) {
  double *A = (double *) MALLOC(individuals * individuals * sizeof(double));
  if (A == NULL) ERR("mem allocation");

  crossprod(SNPxIndiv, snps, individuals, method, centred, normalized,
	      SumGeno, A);
 
  return A;
}
 

SEXP crossprod(SEXP SNPxIndiv) {
   Uint
    *info = GetInfo(SNPxIndiv),
     individuals = info[INDIVIDUALS];
   snpcoding method = check_method((snpcoding) info[METHOD]);
  Uint
     snps = info[SNPS],
    *code = DoAlign(SNPxIndiv, ALIGN_RELATION, method);
  if ((int) info[METHOD] == NA_INTEGER)
     ERR("looks like an uninitialised matrix");
   if (info[METHOD] == Haplo)
     ERR("matrix is a haplotype matrix, not a genomic matric");
   if (info[METHOD] == UnknownSNPcoding) ERR("not a coded Z matrix");
  
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, individuals, individuals));
  option_type *global = &(KEYT()->global);

  crossprod(code, snps, individuals, method,
	      global->genetics.centered == True,
	      global->genetics.normalized == True,
	      SumGeno(info), REAL(Ans));
  
 UNPROTECT(1);
  return Ans;
}
    
    
SEXP file_get(SEXP file) {
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  switch (method) {
  case Shuffle256 : case Shuffle :  case TwoBit : case Packed : case MMAGPU:
    if (useMMAGPU(method)) 
      return file_intern(file, matrix_start_mmagpu, coding2,
    			 CodesPerBlock256, R_NilValue);
    if (useShuffle256(method)) 
      return file_intern(file, matrix_start_shuffle256, coding2,
			 CodesPerBlock256, R_NilValue);
    if (useShuffle(method)) 
      return file_intern(file, matrix_start_shuffle, coding2,
			 CodesPerBlock128(), R_NilValue);
    if (usePacked(method))
      return file_intern(file, matrix_start_packed, coding2,
			 CodesPerBlock128(), R_NilValue);
    return file_intern(file, matrix_start2, coding2,
		       CodesPerBlock2(), R_NilValue);
    
  case Packed256 : return file_intern(file, matrix_start_packed256, coding2,
				      CodesPerBlock256, R_NilValue);

  case Multiply : return file_intern(file, matrix_start_multiply, coding2,
				      CodesPerBlock256, R_NilValue);
  case Multiply256 : return file_intern(file, matrix_start_multiply256, coding2,
					CodesPerBlock256, R_NilValue);
    
  case ThreeBit : return file_intern(file, matrix_start3, coding3,
				     CodesPerBlock3(), R_NilValue);
    
   case Hamming2 : return file_intern(file, matrix_startH2, codingH2,
				     CodesPerBlockH(), R_NilValue);
    
  case Hamming3 : return file_intern(file, matrix_startH3, codingH3,
				     CodesPerBlockH(), R_NilValue);

  case NoSNPcoding: 
  case NoSNPcodingR: return file_intern(file, matrix_start_plain, coding_plain,
					CodesPerBlockPlain(), R_NilValue);

  case AutoCoding :// in R abgesichert
  case Haplo : BUG;    
    
  default : BUG;
    
  }
}




SEXP createSNPmatrix(Uint snps, Uint individuals, snpcoding method) {  
  SEXP Code;
  coding_start_t start = NULL;
  switch(method) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed : 
    if (useShuffle256(method)) start = matrix_start_shuffle256;
    else if (useShuffle(method)) start = matrix_start_shuffle;
    else if (usePacked(method)) start = matrix_start_packed; 
    else start = matrix_start2;
    break;
  case Packed256 : start = matrix_start_packed256; break;
  case MMAGPU : start = matrix_start_mmagpu; break;
  case Multiply : start = matrix_start_multiply; break;
  case Multiply256 : start = matrix_start_multiply256; break;
  case ThreeBit : start = matrix_start3; break; 
  case Hamming2 : start = matrix_startH2; break;    
  case Hamming3 : start = matrix_startH3; break;
  case NoSNPcoding :
  case NoSNPcodingR : start = matrix_start_plain; break;    
  case Haplo : start = create_codevectorHaplo; break; 
  case AutoCoding :  BUG;     
  default :  BUG;
  }

  PROTECT(Code = start(snps, individuals, R_NilValue));      
  UNPROTECT(1);
  return Code;
}


SEXP createSNPmatrix(SEXP SNPs, SEXP Individuals) {
  Uint 
    n = Int0(Individuals),
    snps = Int0(SNPs);
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  if (method == AutoCoding) method = getAutoCodingIntern();
  return createSNPmatrix(snps, n, method);
}



SEXP zeroNthGeno(SEXP SNPxIndiv, SEXP Snps) {
  Uint *info = GetInfo(SNPxIndiv),
    snps = info[SNPS],
    *which = (Uint*) INTEGER(Snps),
    len = length(Snps);  
  for (Uint i=0; i<len; i++) {
    if (which[i] >= snps) ERR("values of 'Snps' out of range.");
  }
    
  snpcoding method = check_method((snpcoding) info[METHOD]);
  switch(method) {
  case Shuffle256 :  case Shuffle : case TwoBit : case Packed : case Packed256 :
  case Multiply : case Multiply256 : case MMAGPU :
    zeroNthGeno2(SNPxIndiv, Snps); break;
  case ThreeBit : zeroNthGeno3(SNPxIndiv, Snps);break;
  case Hamming2 :
  case Hamming3 :  zeroNthGenoH(SNPxIndiv, Snps, method);break;
  case Haplo : zeroNthHaplo(SNPxIndiv, Snps);break; // Haplo, not Geno!!
    
  case NoSNPcoding: //return zeroNthGeno(SNPxIndiv, Snps);
  case NoSNPcodingR: 
  case AutoCoding :  BUG;
  default :  BUG;
  }

  return SNPxIndiv;
}


SEXP allele_freq(SEXP SNPxIndiv) {
  Uint *info = GetInfo(SNPxIndiv);
  snpcoding method = check_method((snpcoding) info[METHOD]);
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed : case Packed256 :
  case Multiply : case Multiply256 : case MMAGPU :
    if (use256()) return allele_freq256(SNPxIndiv);
    if (use128()) return allele_freq128(SNPxIndiv);    
    return allele_freq2(SNPxIndiv);
  case ThreeBit : return allele_freq3(SNPxIndiv);    
  case Hamming2 :return allele_freqH2(SNPxIndiv);
  case Hamming3 :return allele_freqH3(SNPxIndiv);
  case Haplo : ERR("decoding of partial matrix not programmed yet");     
  case NoSNPcoding:  //    return allele_freqPlain(SNPxIndiv); 
  case NoSNPcodingR:
  case AutoCoding : BUG;
      
  default : BUG;
  }
}


SEXP get_matrix_N(SEXP SNPxIndiv, SEXP Snps) {
  Uint *info = GetInfo(SNPxIndiv);
  snpcoding method = check_method((snpcoding) info[METHOD]);
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed : case Packed256 :
  case Multiply : case Multiply256 : case MMAGPU :
    return get_matrixN_2(SNPxIndiv, Snps);    
  case ThreeBit : return get_matrixN_3(SNPxIndiv, Snps);    
  case Hamming2 : return get_matrixN_H2(SNPxIndiv, Snps); 
  case Hamming3 : return get_matrixN_H3(SNPxIndiv, Snps);    
  case Haplo : ERR("decoding of partial matrix not programmed yet"); 
  case NoSNPcoding:
  case NoSNPcodingR: // return get_matrixplain(SNPxIndiv, Snps);
  case AutoCoding :
  BUG;
      
  default : BUG;
  }
}



SEXP matrix_get(SEXP SNPxIndiv) {
  Uint *info = GetInfo(SNPxIndiv);
  snpcoding method = check_method((snpcoding) info[METHOD]);
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed : case Packed256 :
  case Multiply : case Multiply256 : case MMAGPU :
    return get_matrix2(SNPxIndiv);    
  case ThreeBit : return get_matrix3(SNPxIndiv);
  case Hamming2 : return get_matrixH2(SNPxIndiv);
  case Hamming3 : return get_matrixH3(SNPxIndiv);
  case NoSNPcoding:
  case NoSNPcodingR:  return get_matrixPlain(SNPxIndiv);
  case AutoCoding :
  case Haplo : BUG;
      
  default : BUG;
  }
}



SEXP matrix_coding(SEXP M) {
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  if (length(M) == 0) ERR("'M' has length 0.");
  Uint
    snps,
    individuals;
  if (isMatrix(M)) {
    snps = nrows(M);
    individuals = ncols(M);
  } else {
    snps = length(M);
    individuals = 1;
  }  
  
  ToInt(M);
  Ulong total = (Ulong) snps * individuals;
  for (Ulong i=0; i<total; i++) {
    if (Mint[i] > 2) ERR("SNP matrix has only the values 0,1,2");
  }

  SEXP Code;

  PROTECT(Code = createSNPmatrix(snps, individuals, method));
  //  allInfo(Code);
 
  coding_main_t	coding = NULL;
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed : case Packed256 :
  case Multiply : case Multiply256 : case MMAGPU :
    coding = coding2;
    break;
  case ThreeBit : coding = coding3; break;
  case Hamming2 : coding = codingH2; break;
  case Hamming3 : coding = codingH3; break;
  case NoSNPcoding:
  case NoSNPcodingR: coding = coding_plain; break;
  case AutoCoding : // in R abgesichert
  case Haplo : BUG;
  default : BUG;
  }

  coding(Mint, 0, individuals, 0, snps, snps, Code, NULL);
  Uint *info = GetInfoUnchecked(Code);
  
  //  Ulong sumgeno = 0,
  //  endfor = (Ulong) snps * individuals;
  // for (Ulong i=0; i<endfor; i++) sumgeno += (Ulong) Mint[i];

   Uint
      *code = DoAlign(Code, ALIGN_RELATION, method);
   Ulong sumgeno = sumGeno(code, snps, individuals, method);
   StoreSumGeno(sumgeno);

  FREEint(M);
  UNPROTECT(1);
  
  return Code;
}



SEXP fillSNPmatrix(SEXP Z, SEXP Idx, SEXP Vector) {
  ToInt(Idx);
  Uint 
    *infoZ = GetInfo(Z),
    *infoVector = GetInfo(Vector),
    snpsZ = infoZ[SNPS],
    snpsVector = infoVector[SNPS],
    individuals = infoZ[INDIVIDUALS];
  snpcoding     
    methZ = check_method((snpcoding) infoZ[METHOD]),
    methVector = check_method((snpcoding) infoVector[METHOD]);
  Uint
    *z = DoAlign(Z, ALIGN_RELATION, methZ),
    *v = DoAlign(Vector, ALIGN_HAPLO, methZ),
     len = length(Idx);
  if (!isGeno(methZ)) ERR("not a geno matrix that is filled");
  if (snpsZ != snpsVector) ERR("number of snps differ.");
  if (methZ != methVector) ERR("methods of 'SNPxIndiv' and 'values' differ.");
  
  Ulong
    unitsPerIndiv = GetUPI(snpsZ, methZ),
    bytes = BytesPerUnit * unitsPerIndiv;
  
  for (Ulong i=0; i<len; i++, v += unitsPerIndiv) {
    Uint idx = Idxint[i] - 1;
    if (idx >= individuals) ERR("Idx out of bound");
    MEMCOPY(z + idx * unitsPerIndiv, v, bytes);
  }
  FREEint(Idx);
  return R_NilValue;
}



SEXP copyGeno(SEXP CM) {
   Uint
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS];
  
  if (!isGeno(info[METHOD])) ERR("not a geno matrix");
  snpcoding method = check_method((snpcoding) info[METHOD]);
   
  SEXP Code;
  PROTECT(Code = createSNPmatrix(snps, individuals, method));

  Uint
    *from = DoAlign(CM, ALIGN_HAPLO, method),
    *to = DoAlign(Code, ALIGN_HAPLO, method),
    *info2 = GetInfo(Code),
    *addr1 = (Uint*) GetAddress(info, ALIGNADDR0),
    *addr2 = (Uint*) GetAddress(info2, ALIGNADDR0);
  Ulong mem = MemInUnits(info);

  assert(INFO_GENUINELY_LAST== 22 && HEADER == 15);//minicheck ob sich die liste
  // geaendert hat, ohne hier zu korriegieren.
  info2[HEADER] = info[HEADER];
  info2[DOUBLEINDIV] = info[DOUBLEINDIV];
  info2[LEADINGCOL] = info[LEADINGCOL];
 
    
  if (mem != MemInUnits(info2)) BUG;

  intptr_t
    lag1 = (intptr_t) addr1 - (intptr_t) CM,
    lag2 = (intptr_t) addr2 - (intptr_t) Code;
  if (lag1 > 100 || lag1 < 0 || lag2 > 100 || lag2 < 0)
    ERR("Alignment problem -- pls contact maintainer.");
  //  if (lag1 - lag2 ) ERR("Memory was disclocated. Copying is not possible.");
  if (from != addr1 || to != addr2)
    ERR("Alignment problems. Pls contact maintainer.");
  
  MEMCOPY(to, from, mem * BytesPerUnit);
  UNPROTECT(1);
  return Code;
}




#define DotBlockSize 100 // arbitrary number


SEXP file_dot_start(Uint snps,Uint individuals,  SEXP G) {
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
		 Uint start_snp, Uint end_snp, Uint Mnrow,
		 SEXP Ans, double *G) {
  double *Fdot = REAL(Ans);
  for (Ulong a=start_individual; a<end_individual; a++) {
    Uint *pM = M + (a - start_individual) * Mnrow;
    double dummy = G[a];
    for (Uint s=start_snp; s<end_snp; pM++) {
      Fdot[s++] += dummy * (double) *pM;
    }
  }
}

SEXP file_dot(SEXP file, SEXP G) {
  return file_intern(file, file_dot_start, file_dot_do,
			    DotBlockSize, G);
}


SEXP  dot_file_start(Uint snps, Uint individuals, SEXP G) {
  SEXP Ans;
  if ((Uint) length(G) != snps) ERR("vector must have length equal to number of snps");
  PROTECT(Ans=allocVector(REALSXP, individuals));
  for (Uint i=0; i<individuals; REAL(Ans)[i++] = 0.0); // ok
  UNPROTECT(1);
  return Ans;
}


void dot_file_do(Uint *M, Uint start_individual, Uint end_individual, 
		 Uint start_snp, Uint end_snp, Uint Mnrow,
		 SEXP Ans, double *G) {
  double *dotF = REAL(Ans);
  for (Ulong a=start_individual; a<end_individual; a++) {
    Uint *pM = M + (a - start_individual) * Mnrow;
    double dummy = 0;
    for (Uint s=start_snp; s<end_snp; pM++) {
      dummy += G[s++] * (double) *pM;
    }
    dotF[a] = dummy;
  }
}


SEXP dot_file(SEXP file, SEXP G) {
  return file_intern(file, dot_file_start, dot_file_do,
			    DotBlockSize, G);
}



SEXP substract_centered(SEXP SnpXindiv) {
  option_type *global = &(KEYT()->global);
  Uint individuals = ncols(SnpXindiv),
    len = global->genetics.ncentered;
  Ulong snps = nrows(SnpXindiv);
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, snps, individuals));
  double *ans = REAL(Ans),
    *snpXindiv = REAL(SnpXindiv),
    *centered = global->genetics.pcentered;
  assert(centered != NULL && global->genetics.centered==Nan);
  if (snps != len) ERR("length of 'centered' must equal the number of SNPs.");

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
  for (Ulong i=0; i<individuals; i++) {
    Ulong i_snps = i * snps;
    double *s = snpXindiv + i_snps,
      *a = ans + i_snps;
    for (Uint j=0; j<snps; j++) a[j] = s[j] - centered[j];
  }
  UNPROTECT(1);
  return Ans;
}



SEXP vectorGeno(SEXP Vector, SEXP Z) {
  Uint
    *info = GetInfo(Z),
    individuals = info[INDIVIDUALS];
  SEXP Ans;

  PROTECT(Ans = allocVector(REALSXP, individuals));
  double *ans = REAL(Ans);
  
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit: case Packed : case Packed256 :
  case Multiply : case Multiply256 : 
    vectorGeno2(Vector, Z, ans); break;
  default:
    UNPROTECT(1);
    BUG;
  }
  
  UNPROTECT(1);
  return(Ans);
}

SEXP genoVector(SEXP Z, SEXP Vector) {
   Uint
    *info = GetInfo(Z),
    snps = info[SNPS];
  SEXP Ans;
 
  PROTECT(Ans = allocVector(REALSXP, snps));
  double *ans = REAL(Ans);
 
  option_type *global = &(KEYT()->global);
  switch(check_method((snpcoding) global->genetics.method)) {
  case Shuffle256 : case Shuffle : case TwoBit : case Packed : case Packed256 :
  case Multiply: case Multiply256 : case MMAGPU :
    genoVector2(Z, Vector, ans); break;
  default:
    UNPROTECT(1);
    BUG;
  }
  
  UNPROTECT(1);
  return(Ans);
}



Ulong calculateAlignedMem(Ulong memInUnits, snpcoding method,
			  Uint bytesPerBlock) {
  if (is256(method))
    return (2L + (memInUnits - 1L) / UnitsPerBlock256) *  UnitsPerBlock256 - 1L;

  Uint UnitsPerBlock =  bytesPerBlock / BytesPerUnit;
  return (2L + (memInUnits - 1L) / UnitsPerBlock) *  UnitsPerBlock - 1L;
}

void ReUseAs(SEXP Code, snpcoding method) {
  
  SEXP Class;
  Uint *info = GetInfoUnchecked(Code),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    bytesPerBlock = GetBytesPerBlock(method),
    unitsPerIndiv = GetUPI(snps, method);


  if (((int) info[UNITSPERINDIV] != NA_INTEGER &&
       info[UNITSPERINDIV] != unitsPerIndiv))
    ERR4("storage length mismatch (%d; %d; %.20s %.20s)", info[UNITSPERINDIV],
	 unitsPerIndiv, SNPCODING_NAMES[info[METHOD]], SNPCODING_NAMES[method]);
  
  Ulong
    memInUnits = (Ulong) individuals * unitsPerIndiv,
    alignedMem = calculateAlignedMem(memInUnits, method, bytesPerBlock),
    origaligned = AlignedInUnits(info);

  if (origaligned != alignedMem) {
    if (origaligned < alignedMem) ERR("BUG: Reuse of memory impossible");
    Uint *code = (Uint *) INTEGER(Code);
    for (Uint i=alignedMem; i<origaligned; code[i++] = 0L);
  }

  
  info[METHOD] = method;
  if (info[VERSION] != CURRENT_VERSION) ERR("versions do not match");
  if (info[BITSPERCODE] != GetBitsPerCode(method)) ERR("bit-per-code mismatch");
  if (!is256(info[METHOD]) || !is256(method)) {
    if (info[CODESPERBLOCK] != GetCodesPerBlock(method))
      ERR("codes-per-block mismatch");
    if (info[BYTESPERBLOCK] != bytesPerBlock) ERR("block size mismatch");
  }
  if (info[ALIGNMENT] != (is256(method) ? BytesPerBlock256 : bytesPerBlock))
    ERR2("alignment mismatch : %d %d", info[ALIGNMENT], is256(method) ? BytesPerBlock256 : bytesPerBlock);

  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0,
		 mkChar(method == Haplo ? HAPLOMATRIX : GENOMICMATRIX));
  SET_CLASS(Code, Class);
  UNPROTECT(1);
}


int  CECV = 0;
SEXP CreateEmptyCodeVector(Uint snps, Uint individuals, snpcoding method) {
  option_type *global = &(KEYT()->global);
  SEXP Code, Info;
  if (method != Haplo && method != (snpcoding) global->genetics.method)
    ERR2("method mismatch (%.20s; %.20s). Pls contact author",
	SNPCODING_NAMES[method], SNPCODING_NAMES[global->genetics.method]);
  Ulong
    unitsPerIndiv = GetUPI(snps, method),
    memInUnits = (Ulong) individuals * unitsPerIndiv;
  if (method == Hamming2 || method == Hamming3) memInUnits *= 2;
  Ulong
    bytesPerBlock = GetBytesPerBlock(method),
    alignedMem = calculateAlignedMem(memInUnits, method, bytesPerBlock);

  //  if (method ==  )
  if (GetBitsPerCode(method) == 32) { 
    if (snps * (Ulong) individuals != memInUnits || alignedMem != memInUnits)
      BUG;
    PROTECT(Code = allocMatrix(INTSXP, snps, individuals));
  } else {
    if (snps * (Ulong) individuals <= memInUnits && snps > 1000) {
      BUG;
    }
    PROTECT(Code = allocVector(INTSXP, alignedMem)); 
  }
 
  Uint* A = (Uint*) INTEGER(Code);
  for (Ulong k=0; k<alignedMem; A[k++] = 0);

  //  UNPROTECT(1);  return Code;
  
  
  PROTECT(Info = allocVector(INTSXP, INFO_LAST + 1));
  
  Uint *info = (Uint*) INTEGER(Info);
  for (Uint i=0; i<=INFO_LAST; i++) info[i] = NA_INTEGER;
  StoreAlignedInUnits(alignedMem);
  ADDADDRESS(INTEGER(Code), ADDR0);
  setAttrib(Code, Information, Info); // install("information")
  
  info[VERSION] = CURRENT_VERSION;
  info[SNPS] = snps;
  info[INDIVIDUALS] = individuals;  
  info[BITSPERCODE] = GetBitsPerCode(method);
  info[CODESPERBLOCK] = GetCodesPerBlock(method);
  info[BYTESPERBLOCK] = bytesPerBlock;
  info[ALIGNMENT] = is256(method) ? BytesPerBlock256 : bytesPerBlock;
  info[SUMGENO] = info[SUMGENO_E9] = 0;
  info[UNITSPERINDIV] = unitsPerIndiv;
  info[BLOCKEDINFO] = 0;
  info[ZAEHLER] = ++CECV;
  StoreMemInUnits(memInUnits);
  ADDADDRESS((Rint*) algn_general(INTEGER(Code), info[ALIGNMENT]), ALIGNADDR0);
  ADDADDRESS(NULL, RECENTALIGNADDR0);
  ReUseAs(Code, method);

  UNPROTECT(2);
  return Code;
}


void start_info(SEXP Code, SEXP file) {
  Uint *info = GetInfoUnchecked(Code);
  if (file != R_NilValue) {
    Uint *finfo = GetInfoUnchecked(file);
    info[SNPxIND] = finfo[SNPxIND];
  }
  if (PL > 1) {								
    Uint
      mem = MemInUnits(info),
      individuals = info[INDIVIDUALS],
      snps = info[SNPS],
      method = (Uint) check_method((snpcoding) info[METHOD]),
      unitsPerBlock = GetUPI(snps, (snpcoding) method),
      codesperblock = info[CODESPERBLOCK];
    bool inMB = mem > 5000000;
    PRINTF("Data: %d individuals and %d SNPs\nStorage mode: %d block(s) of %d codes\nSize of M: %d %sB.", 
	   individuals, info[SNPS], mem / unitsPerBlock, codesperblock,
	   1 + (mem - 1) / (inMB ? 1048576 : 1024),
	   inMB ? "M" : "k");	       
  }
}


Uint *GetInfo(SEXP Code) { // has negative integers included!!
  SEXP Infos = getAttrib(Code, Information);
  Uint *info = (Uint *) INTEGER(Infos);
  snpcoding method = check_method((snpcoding) info[METHOD]);
  if (info[VERSION] != CURRENT_VERSION)
    ERR2("the stored data format (%d) does not match the current format (%d).", info[VERSION], CURRENT_VERSION );

  if (!is256(method)) {
    if (info[CODESPERBLOCK] != (Uint) GetCodesPerBlock(method)) {
      ERR3("The '%.20s' matrix is unreadable: %d codes per block are used in the data while the system works with %d codes per block",
	 SNPCODING_NAMES[method],
	   info[CODESPERBLOCK], GetCodesPerBlock(method) );
    }
    if (info[BYTESPERBLOCK] != GetBytesPerBlock(method)) {
      ERR("currently, data exchange between different kinds of machines (AVX2/SSE) is not possible");
    }
  }

  if (info[ALIGNMENT] !=
      (is256(method) ? BytesPerBlock256 : info[BYTESPERBLOCK])) {
    //    ERR("Alignment mismatch");
    ERR2("Alignment mismatch : %d %d", info[ALIGNMENT], is256(method) ? BytesPerBlock256 : info[BYTESPERBLOCK]);
 }
  if (info[BITSPERCODE] != GetBitsPerCode(method))
    ERR("Bits per code is inconsistent. Please contact author");
 
       
  if (TYPEOF(Infos) != INTSXP) {
    ERR("obsolete storage mode");
  }
  return info;
}



Uint haplo2geno(Uint *H, Uint snps, Uint individuals,
		snpcoding method, Uint unitsPerIndivH, Uint *A) {

  if (H == A && !is256(method))
    ERR("A memory effecient transformation of a haplotype matrix to a genomic matrix, e.g. from within 'MoBPS', relies on a two-bit-coding assuming for genomic matrices. So, only 'Shuffle(256)', 'Multiply(256)', 'Packed(256)' and 'TwoBit' are allowed snp coding methods");

  h2g_t h2g = NULL;
  switch (method) {
  case Shuffle256 : case Shuffle : case TwoBit :case Packed : case Packed256 :
  case Multiply : case Multiply256 : case MMAGPU :
    if (use256()) h2g = haplo2geno256;
    else if (use128()) h2g = haplo2geno128;
    else h2g = haplo2geno2;
    break;
  case ThreeBit : h2g = haplo2geno3; break;    
  case Hamming2 : h2g = haplo2genoH2; break;
  case Hamming3 : h2g = haplo2genoH3; break;
  case NoSNPcoding:
  case NoSNPcodingR: h2g = haplo2genoPlain; break;    
  case AutoCoding :
  case Haplo : BUG;      
  default : BUG;
  }

  h2g(H, snps, individuals, unitsPerIndivH, A);  
  return sumGeno(A, snps, individuals, method);
  
}


 
void haplo2geno(SEXP H, Uint *code, snpcoding method) {
  Uint
    *info = GetInfo(H),
    individuals = (Uint) info[INDIVIDUALS],
    snps =  (Uint) info[SNPS];

  if (info[METHOD] != Haplo) ERR("not a haplo coding");
  Uint *haplo, unitsPerIndivH;
 
  InitGetHaplo(H, &haplo, &unitsPerIndivH);
  //allInfo(H);
  
  Ulong sumgeno = haplo2geno(haplo, snps, // no PROTECT( neeeded
			     individuals, method, unitsPerIndivH, code);

  if (haplo == code) ReUseAs(H, method);
  StoreSumGeno(sumgeno);
  
}


SEXP haplo2geno(SEXP H) {
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  Uint
   *infoH = GetInfo(H),
   individuals = (Uint) infoH[INDIVIDUALS],
   snps =  (Uint) infoH[SNPS];
  if (infoH[METHOD] != Haplo) ERR("not a haplotype matrix");
  SEXP Code;
  PROTECT(Code = createSNPmatrix(snps, individuals, method));
  Uint
    *UintCode = (Uint*) INTEGER(Code),
    *info = GetInfo(Code),
    *code = (Uint*) algn_general(UintCode, info[ALIGNMENT]);
  haplo2geno(H, code, method);
  UNPROTECT(1);
  return Code;  
}

 

SEXP Debug() { debugging = true; return R_NilValue; }
SEXP StopDebug() { debugging = false;  return R_NilValue; }

SEXP is2BitMethod(SEXP Meth) { // 'Method' is global already
  SEXP Ans;
  int method = INTEGER(Meth)[0];
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = is256(method);
  UNPROTECT(1);
  return Ans;
}
