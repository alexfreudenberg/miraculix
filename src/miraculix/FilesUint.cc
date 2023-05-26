
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

#define BitsPerCode 32
#define MY_VARIANT 32


#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "Template.h"
#include "MXinfo.h"
#include "transform.h"

#include "Files.h"
#include "1bit.h"
#include "2bit.h"
#include "3bit.h"
#include "5codes.h"
#include "4Byte.h"
#include "OneByte.h"
#include "plink.h"
#include "Haplo.h"
#include "intrinsics_specific.h"
 
char EOL = '\n';

#define BitsPerCodeFile 8
#define CodeFileMask ((char) (1 << BitsPerCodeFile) - 1)
#define CodesPerByteFile (BitsPerByte / BitsPerCodeFile)



void file_dot_do(unit_t *M,Long ldM,
		 Long start_snp, Long end_snp, 
		 Long start_individual, Long end_individual,
		 basic_options VARIABLE_IS_NOT_USED *opt,
		 double *G, SEXP Ans) {
  double *Fdot = REALX(Ans);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(opt->cores)) schedule(static)
#endif
  for (Long a=start_individual; a<end_individual; a++) {
    unit_t *pM = M + (a - start_individual) * ldM;
    double dummy = G[a];
    for (Long s=start_snp; s<end_snp; pM++) {
      Fdot[s++] += dummy * (double) *pM;
    }
  }
}



void dot_file_do(unit_t *M,Long ldM,
		 Long start_snp, Long end_snp, 
		 Long start_individual, Long end_individual, 
		 basic_options VARIABLE_IS_NOT_USED *opt,
		 double *G, SEXP Ans) {
  double *dotF = REALX(Ans);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(opt->cores)) schedule(static)
#endif
  for (Long a=start_individual; a<end_individual; a++) {
    unit_t *pM = M + (a - start_individual) * ldM;
    double dummy = 0;
    for (Long s=start_snp; s<end_snp; pM++) {
      dummy += G[s++] * (double) *pM;
    }
    dotF[a] = dummy;
  }
}


void CheckSnpsLessIndiv(Long snps, Long individuals, bool isSNPxInd,
			basic_options *opt) {

  if (snps <= 0 || individuals <= 0 ||
      (snps <= individuals && !opt->skipchecks && opt->helpinfo)) {
    char msg[200];
    SPRINTF(msg, "matrix is supposed to be a '%.20s' x '%.20s' matrix with %ld SNPs and %ld individuals and , which looks odd.\n", 
	    isSNPxInd ? "SNPs" : "individuals",
	    isSNPxInd ? "individuals" : "SNPs",
	    snps, individuals);
    if (snps > 0 && individuals > 0) { WARN0(msg); } else ERR0(msg);
  }
}


void codingSexp(unit_t *M, Long ldM,
		Long start_snp, Long end_snp,
	     Long start_individual, Long end_individual,     
		basic_options *opt, double * G, SEXP Ans) {
  Long *info =  GetInfo(Ans),
    lda = info[LDA],
    indiv = info[INDIVIDUALS];
  unit_t *code =  Align(Ans, opt);
  coding_type coding = (coding_type) info[CODING];
  coding_t coding_4Byte;
  switch(coding) {
  case OneBitGeno : coding_4Byte = coding1_4Byte; break;
  case TwoBitGeno : coding_4Byte = coding2_4Byte; break;
  case Plink : coding_4Byte = codingPlink_4Byte; break;
  case OrigPlink : BUG; // coding_4Byte = codingOrigPlink_4Byte; TODO
    break;
  case ThreeBit : coding_4Byte = coding3_4Byte; break;
  case OneByteGeno : coding_4Byte = codingOneByte_4Byte; break;
  case FourByteGeno : coding_4Byte = codingPlain_4Byte; break;
  case OneByteHaplo : coding_4Byte = codingHaplo1Byte_4Byte ; break;
  default : BUG;
  }
  
  coding_4Byte(M,  ldM,
	      start_snp, end_snp,
	      start_individual,  end_individual, GreaterZero(opt->cores),
	      G, code, indiv, lda);
}


coding_sexp_t matrix_coding_sexp(coding_type coding, int variant,
				 bool bigendian) {
  //  printf("coding_sexp %s\n", CODING_NAMES[coding]);
  assert(has_sexp_coding(coding));
  switch(coding) {
  case TwoBitGeno :
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 : // coding_sexp = coding_2v256; break; Fehler!
    case 128 : if (bigendian) return coding_2v128;
      FALLTHROUGH_OK;
    default : {}
    }
    FALLTHROUGH_OK;
  case OneBitGeno :
  case ThreeBit :
  case OneByteGeno : 
  case Plink :
  case FourByteGeno: 
  case OneByteHaplo:   
    return codingSexp;
    
  case OrigPlink : BUG; // todo
    
  case FiveCodes : return coding5;
  case FileDot : return file_dot_do;
  case DotFile :  return dot_file_do;
  case AutoCoding :// in R abgesichert
  default:
    //printf("coding_sexp %s\n", CODING_NAMES[coding]);
    BUG;
  }
  BUG;
  return NULL;
}



SEXP file_intern(SEXP file, coding_type coding, int variant,
		 option_type *global, utilsoption_type *utils,
		 bool only_complete,
		 SEXP G) {
  basic_options *opt = &(utils->basic);
  tuning_options *tuning = &(global->tuning);
  SEXP XX = getAttribPointer(file, Filecoding);
  char *file_coding = (char*) CHAR(STRING_ELT(XX, 0));
  Long n = STRLEN(file_coding); // kann 4 oder 5 sein
  char ch, oldch,
    A = file_coding[0], // coding 0
    B = file_coding[1], // coding 1
    C = file_coding[n-3], // coding 1 or 2
    NA = file_coding[n-2], // always last but one
    SEP = file_coding[n-1]; // always last
  bool haploFile = B == C;

  //  printf("%c.%c.%c.%c.%c.%d %d; %s \n", A, B, C, NA, SEP, n, haploFile,	 CODING_NAMES[coding]);

      
  coding_sexp_t
    coding_sexp = matrix_coding_sexp(coding, variant, opt->bigendian);
  Long codesperblock = GetCodesPerBlock(coding, variant); // OK
  Long *info = GetInfo(file, 2);
  if (info == NULL || info[FILE_LEADINGCOL] == NA_LONG) BUG;
  Long
    leadingcols = info[FILE_LEADINGCOL],
    isSNPxInd = info[FILE_SNPxIND_READ_BYROW],// = TRUE for plinkbinary
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    LDAbitalign = coding > LastUserCoding ? NA_LONG : GetLDAbitalign(coding),
    doubledindividuals = info[FILE_DOUBLEINDIV];
  stopIfNotInt(info[FILE_HEADER]);
  int header = (int) info[FILE_HEADER]; // can be indeed negative !!
  
  char *filename = (char*) CHAR(STRING_ELT(file, 0));

  if (header < 0 && !only_complete) {
    // only for plinkbinary: snps & individuals already known !! 
    return file_binary(filename,
		       file_coding,
		       -header, isSNPxInd,
		       snps, individuals,
		       codesperblock,// OK
		       // output
		       coding, variant, LDAbitalign,
		       coding_sexp, global, utils,
		       G
		       );
  }

  if (tuning->missingsFully0) {
    ERR0("missings not programmed yet\n");
    // unklar wie der Stream am besten programmiert wird
    // getPlinkMissings(plink, file, opt); s. plink4Byte.cc
  }

 

  // BEAGLE 3 
  // beides: ani x snps und snps x individ
  // A B A B ?
  // ? als voraussetzung nicht vorhanden
  // 1 header-Zeile
  // 2 erste Spalten weg
  //
  // plink: kein header
  // 4 spalten, die nicht brauchen
  //
  // IndivPerCol = TRUE for plinkbinary, i,e,
  //   each row of plink contains all SNPS of an individual,
  //   i.e. 'rows' below equals the number of individuals

  FILE *fp;
  unit_t *matrix = NULL; // snps x individuals
  Long
    nrow_matrix, ncol_matrix,
    rows, // rowsM1, 
    cols, colsM1,
    plusrow,  haplorow, haplocol;
  if ((fp = fopen(filename, "r")) == NULL) {
    ERR1("file '%.50s' could not be opened", filename);
  }

  // determine number of rows and columns
  for (Long i=0; i<header; i++) while ((char) fgetc(fp) != EOL);

  // determine size of matrix:
  rows = 1;
  cols = 0;
  oldch = SEP;
  while ((ch  = (char) fgetc(fp)) != EOL) {
    //   printf("%c.%c. (%d %d) %d %d\n", ch, oldch, (int) ch, (int) SEP,  ch == SEP ,  oldch != SEP);
    if (ch == SEP && oldch != SEP) cols++;
    oldch = ch;
  }
  if (oldch != SEP) cols++;
  cols -= leadingcols;

  while ((ch  = (char) fgetc(fp)) != EOF) { 
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
  
  if (haploFile) {
    // printf("infiv =%d %d %d; (%d x %d) %d\n", individuals, doubledindividuals, isSNPxInd, rows , cols, leadingcols);
    if (doubledindividuals) {
      Long dummy = individuals >> 1; 
      if (dummy << 1 != individuals) ERR0("odd number of values (individuals)");
      individuals = dummy;
    } else {
      Long dummy = snps >> 1; 
      if (dummy << 1 != snps) ERR0("odd number of values (SNPs)");
      snps = dummy;
    }
  }
  
  CheckSnpsLessIndiv(snps, individuals, isSNPxInd, opt);

  if (info[CODING] == NA_LONG || info[CODING] == coding) {
    CompleteCodeVector(snps, individuals, coding, variant, LDAbitalign, false,
		       global, utils, NULL, NULL, file);
    if (only_complete) return file;
  }

  SEXP Ans = PROTECT(createSNPmatrix(snps, individuals, coding, variant,
				     LDAbitalign,
				     false, G, global, utils));
 

  colsM1 = cols - 1;
  //rowsM1 = rows - 1;

  // read in table
  
  if (isSNPxInd) {
    plusrow = 1;
    nrow_matrix = codesperblock; // OK
    ncol_matrix = individuals;
  } else {
    // plusrow = snps;
    plusrow = codesperblock;// OK
    ncol_matrix = snps;
    nrow_matrix = 1;
  }
  if (haploFile) {
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
  
  Long
    DeltaMatrix = 0,
    matrixLDA = SnpNatLDA(nrow_matrix, MY_CODING),
    colsLDA = (1L + (ncol_matrix - 1L) / CodesPerByteFile) * CodesPerByteFile,
    idxIncr = isSNPxInd ? matrixLDA : 1,
    r=0,
    matrix_size = matrixLDA * colsLDA;
  double *dG = NULL;
  if (LENGTH(G) > 0) dG = REALX(G);

  rows /= haplorow;
  cols /= haplocol;

  if (isHaplo(coding)) {
    assert(CodesPerByteFile == 1);
    assert(colsLDA == ncol_matrix);
    assert(coding == OneByteHaplo);
    DeltaMatrix=matrix_size;
    matrix_size *= 2;
  }


  // printf("\n########### matrixLDA %d %d %d\n", matrixLDA, DeltaMatrix, isHaplo(coding));
  
  matrix = (unit_t*) CALLOC(matrix_size, BytesPerUnit);

  // jump header lines
  fp = fopen(filename, "r");
  for (Long i=0; i<header; i++) while ((ch=(char) fgetc(fp)) != EOL);
  
  for (Long rowidx=0; r<rows; r++, rowidx+=plusrow) {//# of (multiple) rows
    for (Long hr=0; hr<haplorow; hr++) { 
      unit_t *M = matrix + hr * DeltaMatrix;
      // genuine loop only if haplo types are given rowwise
      Long idx = rowidx;
      
      // jump unnecessary leading stuff
      while ((ch = (char) fgetc(fp)) == SEP);
      for (Long l=0; l<leadingcols; l++) {
	while ((ch=(char) fgetc(fp)) != SEP && ch!=EOF);
	if (ch == EOF) ERR0("unexpected end of file");
	while ((ch = (char) fgetc(fp)) == SEP);
      }

      for (Long l=0; l<cols; l++, idx+=idxIncr) { // # of (multiple cols
	for (Long hc=0; hc<haplocol; hc++) {
	  // genuine loop only if haplo types are given column wise
	  
	  if (ch != A) { // == 0
	    if (ch == B) {
	      M[idx + hc * DeltaMatrix]++;
	    } else if (ch == C) {
	      M[idx+ hc * DeltaMatrix] += 2;
	    } else if (ch == NA) ERR0("missing value detected.")
	    else {
	      PRINTF(">%c< row=%ld of %ld rows, col=%ld of %ld columns; plus=%ld hc=%ld haplocol=%ld\n", ch, r, rows, l, cols, nrow_matrix, hc, haplocol);
	      ERR0("Unknown character detected. Note that missings are not allowed here.");
	    }
	    //assert(matrix[idx] <=2);
	  } 
	  if (l<colsM1 || hc!=haplocol-1)
	    while((ch = (char) fgetc(fp)) == SEP); 
	}
      } // cols l
      
      // search for end of line
      if (ch != EOL && ch != EOF)
	while ((ch=(char) fgetc(fp)) != EOL && ch!=EOF);
    } // hr
   
    UPDATE;
    
  } // rows r
  
  CLOSE;
  
  UNPROTECT(1);

  return Ans;
}


