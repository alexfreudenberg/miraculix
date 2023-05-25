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



/******** API **********

TO DO


 ****** END API *********/

 

#include <stdio.h>
#include <stdlib.h>
#include "Automiraculix.h"
#include "5codes.h"


void get_compressed_freq(void *compressed, double *f) {
  getFreq5(compressed, f);
}


// if value unknown set it to 0.
void setOptions_compressed(int use_gpu, // 0 == use CPU
			   int cores,   // 0 == use global env
			   int floatLoop, // 0 == use doubles
			   int meanSubstract, // 1 == increase numerical precis.
			   int ignore_missings,// 0==respect missings (untested)
			   int do_not_center, // 0 == do centering
			   int do_normalize,  // 0 == do not normalize
			   int use_miraculix_freq, // 0==use Fortran freq 
			   int variant,  // 0 == use internal {0,32,128,256}
			   int print_details) {//1 == get some messages


  //  printf("floatLoop = %d \n", floatLoop); exit(99);
    printf("get started\n");
    getStartedOptions();
    //  printf("GPU setting options in setOptions_compressed\n");
    setOptions5(use_gpu, cores, floatLoop,	     
	      meanSubstract == 1 || meanSubstract == 2,  // meanV
	      meanSubstract == 1 || meanSubstract == 3,  // meanSxI
	      !ignore_missings, // ! 
	      !do_not_center,   // !
	      do_normalize,
	      !use_miraculix_freq,
	      variant,
	      print_details
	      );
  
}


int is(char *trans) {
  if (*trans == 'T' || *trans == 't' || *trans == 'Y' || *trans =='y') return 1;
  if (*trans != 'N' && *trans != 'n') exit(99);
  return 0;
}


void plink2compressed(char *plink,
                      char *plink_transposed, 
		      int snps, int indiv,
                      double *f, 
                      int max_n,
		      void**compressed) {

  // setOptions_compressed(0, 6, 0, 0, 0, 0, 0, 0, 1);

   int nbytes = indiv /4;   //printf("p2c snps=%d indiv=%d max_n=%d bytes/indiv=%d\n", snps, indiv, max_n, nbytes);
  // exit(99);
 //  printf("c byte 1:5 = %d %d %d %d %d\n", *plink, *(plink+1), *(plink+2), *(plink+3), *(plink+4));

  plinkToSxI(plink, plink_transposed, snps, indiv, FiveCodes, f, max_n,
	     compressed);
  return; 
}

void dgemm_compressed(char *trans, // 'N', 'T'
                      void *compressed,
                      int n, // number of columns of the matrix B
                      double *B,	
                      int Ldb, // how should it be size_t ldb
                      double *C, int Ldc) {

  if (is(trans))
    genoVector5api(compressed, B, n, Ldb, C, Ldc);
  else
    vectorGeno5api(compressed, B, n, Ldb, C, Ldc);
  return;
}

void dgemm_plink(char* trans, // 'N', 'T'                     
		 char *plink,      // first all indiv of one snp
		 char *plink_transposed,// first all snps of one indiv
		 int snps, int indiv,
		 double *f, 
		 int n, // number of columns of the matrix B 
		 double *B,                               
		 int Ldb,    //how should it be size_t ldb
		 double *C,
		 int Ldc) {
  
  if (is(trans)) {
    int tmp = snps; snps=indiv; indiv=tmp;
    plink_transposed = plink;
  } 

  vectorGenoPlinkApi(plink_transposed, snps, indiv, f,  B, n, Ldb, C, Ldc); 
  return;
}



// transsparse = 'N', transcompressed = 'N'
void sparse_times_plink(char *transsparse,// N: matrix as is (sparse row format)
			char *transcompressed, // N: sparse * plink 
			char *plink,      // first all indiv of one snp
			char *plink_transposed,// first all snps of one indiv
			int snps, int indiv,
			int nIdx,// 2nd dimension of sparse; =length(rowIdxB)-1
			int *rowIdxB, // cf. CSR/CSR/Yale 
			int *colIdxB, // cf. CSR/CSR/Yale 
			double *B,    // values of sparse
			double *C,    // result (nIdx x snps)
			int Ldc       // typically =nIdx
			) {

  if (is(transcompressed)) {
    int tmp = snps; snps=indiv; indiv=tmp;
    plink = plink_transposed;
  }
  
  sparseTGenoPlinkApi(plink, indiv, snps, 
		      B, nIdx, rowIdxB, colIdxB, is(transsparse), C, Ldc);

  return;
}

void free_compressed(void **compressed){
   free5(compressed);
}

  
