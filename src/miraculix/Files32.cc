
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


#define BitsPerCode 32 // ownCoding: FourByteGeno in file_binary
#define MY_VARIANT 32
 
#include "Basic_miraculix.h"
#include "intrinsics_specific.h" 
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "Template.h"
#include "MX.h"
//#include "transform.h"

#include "Files.h"
//#include "1bit.h"
#include "2bit.h"
//#include "3bit.h"
//#include "4Byte.h"
//#include "OneByte.h"
//#include "bitUint.h"

#define BitsPerCodeBinaryFile 2
#define CodeBinaryFileMask ((char) ((1 << BitsPerCodeBinaryFile) - 1))
#define CodesPerByteBinaryFile (BitsPerByte / BitsPerCodeBinaryFile)


SEXP file_binary(char *file, char *file_coding, int header, bool isSNPxInd, 
		 // output
		 Long snps,
		 Long individuals,
		 Long codesperblock, // OK
		 coding_type coding, int variant,
		 Long  LDAbitalign,
		 coding_sexp_t coding_sexp,
		 option_type *global, utilsoption_type *utils,
		 SEXP G
		 ) {
  // for plinkbinary files

  basic_options *opt = &(utils->basic);
  tuning_options *tuning = &(global->tuning);
  char 
    //   A = file_coding[0],
    B = file_coding[1],
    C = file_coding[2]
    //  NA = file_coding[3],
    //  SEP = file_coding[4]
    ;
  bool haplo = B == C;
  if (haplo) ERR0("haplo not recognized for binary files");
   
  FILE *fp;
  unit_t *matrix = NULL; // snps x individuals
  Long bytesperblock = GetBytesPerBlock(coding, variant);// OK
  Long nrow_matrix, ncol_matrix,
    rows, // rowsM1, 
    plusrow;
  Long mainvariant = main_variant(variant);

  if (tuning->missingsFully0) {
    ERR0("missings not programmed yet\n");
  }
 

  CheckSnpsLessIndiv(snps, individuals, isSNPxInd, opt);
  
  if (isSNPxInd) {
    plusrow = 1;
    nrow_matrix = codesperblock;  // OK
    rows = snps;
    ncol_matrix = individuals;
  } else {
    // plusrow = snps;
    plusrow = codesperblock;// OK
    rows = individuals;
    ncol_matrix = snps;
    nrow_matrix = 1;
  }
  
  SEXP Ans = PROTECT(createSNPmatrix(snps, individuals, coding, variant,
				     LDAbitalign, false, G, global, utils));
  double *dG = NULL;
  if (LENGTH(G) > 0) dG = REALX(G);

  coding_sexp = coding==TwoBitGeno && !isSNPxInd ? NULL : coding_sexp;
  if (coding_sexp==NULL && nrow_matrix != 1) BUG;

  TrafoFileValues_t trafoFileValues = (coding == FileDot || coding == DotFile)
    ? NULL
    : mainvariant == 512 ? TrafoFileValues_512
    : mainvariant == 256 ? TrafoFileValues_256
    : mainvariant == 128 ? TrafoFileValues_128
    : TrafoFileValues_64; 
  
  unit_t *aligned = Align(Ans, opt);
  Long *info = GetInfo(Ans),
    allUnits = Units(snps),						
    r=0,								
    OneByteUnits = 1L + (ncol_matrix - 1L) / CodesPerByteBinaryFile,	
    /* as TwoBitGeno also codes in 2 bits, in natural numbers: */	
    MaxByteAlign = MAX_LDA_BITALIGN / BitsPerByte,			
    blocks = (1 + (OneByteUnits - 1) / bytesperblock),	// OK		
    colsLDA = OneByteUnits * CodesPerByteBinaryFile,			
    matrixLDA = SnpNatLDA(nrow_matrix, MY_CODING),
    matrix_size = matrixLDA * colsLDA,					
    idxIncr = isSNPxInd ? matrixLDA : 1;				
  Long lda = info[LDA];						
  matrix = (unit_t*) CALLOC(matrix_size, BytesPerUnit);
  char *buffer0 = (char*) CALLOC(blocks * bytesperblock + MaxByteAlign, 1);// OK
  char *buffer =(char*)((1L + (uintptr_t) buffer0 / MaxByteAlign)*MaxByteAlign);

  if ((fp = fopen(file, "r")) == NULL)
    ERR1("file '%.50s' could not be opened", file);
  for (Long i=0; i<header; i++) fgetc(fp);

  for (Long rowidx=0; r<rows; r++, rowidx+=plusrow) {			
    /* nmbr of (multiple) rows */						
    Long buffered = fread(buffer, 1, OneByteUnits, fp);			
    if (buffered < OneByteUnits) ERR0("file size too small");		
    trafoFileValues((unit_t*) buffer, blocks); // OK

    //   printf("\n########### (coding_sexp %d %d %d\n", coding_sexp == NULL, 0, 0);

    if (coding_sexp == NULL) { /* short cut coding2 && !isSNPxInd */	
      unit_t *code = aligned + r * lda;
      if (opt->bigendian) MEMCOPY(code, buffer, allUnits * BytesPerUnit);
      else {
	char *pb = buffer;
	for (Long j=0; j<allUnits; j++) {
	  assert(code[j] == 0);
	  for (int k=0; k<BytesPerUnit; k++, pb++)
	    code[j] |= ((unit_t) *pb) << (k * BitsPerByte);
	}
      }
      rowidx = -plusrow;							
    } else {								
      for (Long l=0; l<OneByteUnits; l++) {/*# of (multiple OneByteUnits*/ 
	Long idx = rowidx;						
	int k = 0;							
	char ch = buffer[l];						
	for (int j = 0; j < (int) CodesPerByteBinaryFile;			
	     j++, idx+=idxIncr, k+=BitsPerCodeBinaryFile) {		
	  matrix[idx] = (ch >> k) & CodeBinaryFileMask;		
	} /* j */							
      } /* OneByteUnits l */						
      UPDATE;								
    }							
  } /* rows r */							
  
   
  CLOSE;

  FREE(buffer0);							
  UNPROTECT(1);
  return Ans;
}


