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

#include <inttypes.h>
#include "def.h"
#include "Basic_miraculix.h"
//#if defined compatibility_to_C_h
#include "xport_import.h"
#include "zzz_RFU.h"
#include "kleinkram.h"
#include "transform.h"
#include "options.h"
#include "Template.h"
#include "plink.h"
#include "MXinfo.h"
#include "haplogeno.h"
#include "Vector.matrix.h"
#include "5codes.h"
#include "time.h"
#include "intrinsics_specific.h"

void get_started(option_type **g, utilsoption_type **u);



// SPARSE
// make standalone && make standalone && Wageningen/run_gcc snps=11200 ind=12300 cores=10 coding=4 mode=3 SNPmode=3 repetV=32 trials=10 centered=0 variant=256 sparse=100000 cmp=1 rowMode=2 colMode=2 fillMode=2
// L 361, plinkUint pr a gma for
//   without static:  18/100
//   with static: 18/100
//   with static ordered: 19/110

// make standalone && make standalone && Wageningen/run_gcc snps=112000 ind=123000 cores=20 coding=4 mode=3 SNPmode=3 repetV=32 trials=10 centered=0 variant=256 sparse=1 cmp=1  ## katastrophal schlecht bei cores=20 -- mit cores=1 viel besser!

// make standalone && make standalone && Wageningen/run_gcc snps=112000 ind=123000 cores=20 coding=4 mode=3 SNPmode=3 repetV=32 trials=10 centered=0 variant=256 sparse=1 cmp=1 ## katastrophal schlecht bei cores=20



// make standalone && Wageningen/run_gcc snps=5024 ind=4100 cores=5 coding=5 mode=3 SNPmode=3 repetV=2 trials=1 cmp=9

//make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=10 coding=5   mode=3 SNPmode=3 repetV=32 trials=1 centered=0 variant=256


// make standalone && Wageningen/run_gcc snps=12800 ind=13800 cores=4 coding=4 mode=3 SNPmode=3 center=0 repetV=16 trials=1 variant=257


// make standalone && Wageningen/run_gcc snps=2400 ind=3800 cores=4 coding=4 mode=3 SNPmode=0 center=0 repetV=16 trials=1 variant=257 cmp=10



// touch miraculix/src/5codesUint.cc && make standalone && Wageningen/run_gcc snps=150000 ind=150000 cores=6 coding=5 mode=3 SNPmode=3 repetV=100 trials=3 centered=0 variant=32

  
// make standalone && Wageningen/run_gcc snps=150000 ind=150000 cores=20 coding=5 mode=3 SNPmode=3 repetV=100 trials=3 centered=0 variant=32


// make standalone && Wageningen/run_gcc snps=150000 ind=150000 cores=20 coding=5 mode=3 SNPmode=3 repetV=10 trials=3 centered=0 variant=32

// make standalone && Wageningen/run_gcc snps=150000 ind=150000 cores=20 coding=5 mode=3 SNPmode=3 repetV=100 trials=3 centered=0 variant=32



// make standalone && valgrind --tool=memcheck --leak-check=full --num-callers=20 Wageningen/run_gcc snps=50241 ind=410000 cores=5 coding=5 mode=3 SNPmode=3 repetV=1 trials=5




// make standalone && Wageningen/run_gcc snps=50241 ind=410000 cores=5 coding=5 mode=3 SNPmode=3 repetV=1 trials=5 


// make standalone && Wageningen/run_gcc snps=50241 ind=410000 cores=6 coding=5 mode=3 SNPmode=3 repetV=1 trials=1 variant=128



// make standalone && Wageningen/run_gcc snps=12800 ind=13800 cores=4 coding=5 mode=3 SNPmode=3 repetV=1 trials=1 variant=128




//  make standalone && Wageningen/run_gcc snps=102400 ind=51200 cores=2 SNPmode=5 mode=1 coding=2 cmp=0


//  make standalone && Wageningen/run snps=262144 ind=512 variant=32 mode=1


// make standalone && (cd Wageningen; make rmR; ls -al Wageningen) && cd ~/svn/miraculix

//  make standalone && Wageningen/run snps=2621440 ind=5120 variant=32

// make standalone && (cd Wageningen; make rmR; ls -al Wageningen; rm *.o && make) && cd ~/svn/miraculix && Wageningen/run snps=2621440 ind=5120 variant=32

//  make StandAloneOfficial && (cd Wageningen && make rmR && make clean) && tar czvf Wageningen.tgz Wageningen && ll Wageningen.tgz


/// make StandAloneOfficial && (cd Wageningen; make rmR; ls -al Wageningen; rm *.o && make) && cd ~/svn/miraculix && Wageningen/run snps=262144 ind=512 variant=32

#define CpBplink 4
unit_t *OneByte2origplink(SEXP SxI, utilsoption_type *utils) {
  basic_options *opt = &(utils->basic);
  Long
    *info = GetInfo(SxI),
    snps = info[SNPS],
    indiv = info[INDIVIDUALS],
    lda = info[LDA],
    ldAns = 1L + (indiv - 1L) / CpBplink; // in bytes
  coding_type coding = (coding_type) info[CODING];
  if (coding != OneByteGeno) BUG;
  _1Byte *code = (_1Byte*) Align(SxI, opt);
  _1Byte *Ans = (_1Byte *) CALLOC(ldAns * snps, 1);
    
  for (Long j=0; j<snps; j++) {
    _1Byte *nc = Ans + ldAns * j;
    _1Byte *c = code + j;
    int shift = 0;
    for (Long i=0; i<indiv; i++) {
      _1Byte cc = c[j * lda];
      nc[i / CpBplink] |=  (_1Byte) (HUMAN2PLINK(cc) << shift);
      shift = (shift + 2) % BitsPerByte;      
    }
  }
  return (unit_t*) Ans;
}

unit_t *OneByte2Tplink(SEXP SxI, utilsoption_type *utils) {
  basic_options *opt = &(utils->basic);
  Long
    *info = GetInfo(SxI),
    snps = info[SNPS],
    indiv = info[INDIVIDUALS],
    lda = info[LDA],
    ldAns = DIV_GEQ(snps, CpBplink); // in bytes
  coding_type coding = (coding_type) info[CODING];
  if (coding != OneByteGeno) BUG;
  _1Byte *code = (_1Byte*) Align(SxI, opt);
  Uchar *Ans = (_1Byte *) CALLOC(ldAns * indiv, sizeof(_1Byte));
    
  for (Long i=0; i<indiv; i++) {
     Uchar *nc = Ans + ldAns * i;
     _1Byte *c = code + lda * sizeof(unit_t) * i;
    int shift = 0;
       for (Long j=0; j<snps; j++) {
       _1Byte cc = c[j];
      nc[j / CpBplink] |=  (Uchar) (HUMAN2PLINK(cc) << shift);
      shift = (shift + 2) % BitsPerByte;      
    }
  }
  return (unit_t*) Ans;
}

SEXP CreateRandomSNPmatrix(Long snps, Long indiv,
			   coding_type coding, int variant, int mode,
			   option_type *global,   utilsoption_type *utils) {
  basic_options *opt = &(utils->basic);
#define npr 17
  char pseudo_random[npr] = {0L,1L,1L,0L,2L,  2L,2L,2L,1L,0L,  0L,1L,0L,0L,0L, 1L,2L};
  SEXP SxI = CreateEmptyCodeVector(snps, indiv,/* PROTECT() */
				   coding, variant, 0, false, global, utils);
  int cores = GreaterZero(opt->cores);
  Long *info = GetInfo(SxI),
    lda = info[LDA];
  unit_t
    *code = Align(SxI, opt);
  
  // (pseudo-random) filling of SNP-matrix
  if (coding == OneByteGeno || coding == FourByteGeno) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
    for (Long i=0; i <indiv; i++) {
      Long k = i  * snps;
      if (coding == OneByteGeno) {
	_1Byte *c = (_1Byte*) (code + i * lda);
	switch(mode % 10) {
	case 0 : MEMSET(c, 0, snps); break;
	case 1 : MEMSET(c, 1, snps); break;
	case 2 : MEMSET(c, 2, snps); break;
	case 3 :
	  for (Long s=0; s<snps; s++, k++) {
	    c[s] = (_1Byte) ((pseudo_random[k % npr] + k + i) % 3);
	  }
	  if (i==0){
	    c  = (_1Byte*) code;
	    c[0] = 2;
	  }
	  break;
	default : for (Long s=0; s<snps; s++, k++) {
	    c[s] = pseudo_random[k % npr];
	  }
	}
      } else if (coding == FourByteGeno) {
	unit_t *c = code + i * lda;
	for (Long s=0; s<snps; s++, k++) {
	  c[s] = pseudo_random[k % npr];
	}
      }
    }
  } else {
    MEMSET(INTEGER(SxI), 0, info[ALIGNEDUNITS] * BytesPerUnit); // OK
    PRINTF(" *** empty matrix created\n"); 
    // ERR1("only OneByteGeno and FourByteGeno allowed in CreateRandomCodeVector. Got %s.", CODING_NAMES[coding]);
  }

  
  if (coding == OneByteGeno) {
    mode /= 10;
    switch(mode) {
    case 1 : 
    case 2 : 
    case 3 : 
    case 4 : {
      int value = mode - 2;
      Long n = value >=0 ? 0L : (snps + indiv - 2);
      Long *M = NULL;
      if (n > 0) {
	SEXP Miss = allocMatrixExtra(LONGSXP, 2, n);/* PROTECT() */
	M = LONG(Miss);
 	info[MISSINGS] = n;
	setAttrib(SxI, Missings, Miss);
	value = 0;
      }
      //      printf("value=%d n=%d\n", value, n);
      Long m = 0;
      for (int i = 1; i<indiv; i++) {
	_1Byte *c  = (_1Byte*) (code + i * lda);
	c[0] = (_1Byte) value;
	if (n > 0) {
	  M[m++] = 0;
	  M[m++] = i;
	}
      }
      _1Byte *c  = (_1Byte*) code;
      for (int j = 1; j<snps; j++) {
	c[j] = 0;
	if (n > 0) {
	  M[m++] = j;
	  M[m++] = 0;
	}
      }
    }
      break;
    default : { } 
    }
  }
    
  return SxI;
}

double *RandomVector(Long n, int mode) {
#define npv 13
  double pseudo_v[npv] = {0,1.7,-2.9,2.7, -1.2,0.2,0.3,1.2,0.5, -3,-3.1,1};
  double *ans = (double*) CALLOC(n, sizeof(double));
  switch(mode) {
  case 0 : break;
  case 1 : for (int i=0; i<n; i++) ans[i] = 1.0 ; break;
  case 2 : for (int i=0; i<n; i++) ans[i] = 2.0 ; break;
  case 3 : for (int i=0; i<n; i++) {
      ans[i] = 1.9 * (double) random() / (POW(2.0, 31.0) - 1.0) - 1.0;
      if (ans[i] > 1.0 || ans[i] < -1.0) {PRINTF("%f\n", ans[i]); BUG;}
    }
    break;
    
  case 4 : for (int i=0; i<n; i++) ans[i] = -1.0 + TWOPI * i - SQRT2 * TRUNC((TWOPI * i) / SQRT2); break;
  default : for (Long i=0; i<n; i++) ans[i] = pseudo_v[i % npv];
  }
  return ans;
}


void input(int argc, char *argv[],
	   const char *argList[], int maxarg, Long argdefault[],
	   bool do_print) {
  for (int k=1; k<argc; k++) {
    int i;
    int j = Match((char *)"argument of 'run'", argv[k], argList, maxarg,
		  '=', &i, false);
    assert(j >= 0);
    if (!STRCMP(argList[j], "coding")) {
      switch(argv[k][i + 2]) {
      case 't' : case 'T' : argdefault[j] = TwoBitGeno; continue;
      case 'p' : case 'P' : argdefault[j] = Plink; continue;
      case 'f' : case 'F' : argdefault[j] = FiveCodes; continue;
      case 'u' : case 'U' : argdefault[j] = UnknownSNPcoding; continue;
      default :  {} 
      }
    } else if (!STRCMP(argList[j], "floatLoop")) {
      //  argdefault[Match((char*) "meanSubst", argList, maxarg)] = true;
      //  argdefault[Match((char*) "SxInotVmean", argList, maxarg)] = false;
    }
 
    argdefault[j] = strtoimax(argv[k] + i + 2, NULL, 10);
  }
  if (do_print) {
    int end = argdefault[3] ? maxarg : 3;
    PRINTF("\n%s ", argv[0]);
    for (int k=0; k<end; k++) PRINTF("%s=%ld ", argList[k], argdefault[k]);
    PRINTF("\n");
  }
}


void RandomSparse(int nIdx, int max_ncol,
		  int rowMode, int colMode, int fillMode,
		  double **valueB, int **RowIdxB, int **ColIdxB) {
  //  printf("entering randomsparse\n");
  *RowIdxB = (int*) MALLOC(sizeof(**RowIdxB) * (nIdx + 1));
  int *rowIdxB = *RowIdxB;
  rowIdxB[0] = 0;

  if (rowMode <= 0) for (int i=1; i<=nIdx; i++) rowIdxB[i] = -rowMode;
  else for (int i=1; i<=nIdx; i++) rowIdxB[i] = (int)(random() % (rowMode + 1));
  for (int i=1; i<=nIdx; i++) rowIdxB[i] += rowIdxB[i-1];
  int total = rowIdxB[nIdx];

  int *colIdxB = *ColIdxB = (int*) MALLOC(sizeof(**ColIdxB) * total);
  switch (colMode) {
  case 0 : for (int i=0; i<total; i++) colIdxB[i] = i % max_ncol;
    break;
  case 1 : for (int i=0; i<total; i++)
      colIdxB[i]= (int) ((i* (Long)337733) % max_ncol);
    break;
  case 2 : for (int i=0; i<total; i++)
      colIdxB[i]= (int)(random() % max_ncol);
    break;
  default : 
    for (int i=0; i<total; i++)
      colIdxB[i]= i % 2 == 0 ? max_ncol - 1 : i % 3 == 0 ? 0 : i % max_ncol;
  }
  
  double *B = *valueB = (double*) MALLOC(sizeof(**valueB) * total);
  switch (fillMode) {
  case 0 : for (int i=0; i<total; i++) B[i] = (double) (i + 1);
    break;
  case 1 : case 2: case 3:
    for (int i=0; i<total; i++) B[i] = (double) fillMode;
    break;
  default :
    for (int i=0; i<total; i++)
      B[i] =  1.9 * (double) random() / (POW(2.0, 31.0) - 1.0) - 1.0;
  }
  //  printf("leaving randomsparse\n");
 }

void Summary(double *x, Long len, Long repet) {
  double total = 0.0;
  for (Long r=0; r<repet; r++, x+=len) {
    double sum =0.0;
    for (Long i=0; i<len; i++) {
      sum += x[i];
    }
    PRINTF("%20.4fs ", sum);
    total += sum;
  }
  if (repet>1) PRINTF("\n18.4%ft", total); // else PRINTF("\n");
}


void Summary(const char*txt, double *x, Long len) {
  PRINTF("%s : ", txt);
  Summary(x, len, 1);
}


void Summary(char *x, Long len) {
  double sum =0.0;
  for (Long i=0; i<len; i++) {
    sum += (double) x[i];
  }
  PRINTF("%6.4f (byte basis)\n ", sum);
}

#if __INTEL_COMPILER
#define TRIALS 5
#else
#define TRIALS 1
#endif


// *********************************************************************
// *** MAIN
// *********************************************************************
#define RIGHT 0x55555555
#define LEFT 0xAAAAAAAA



int main(int argc, char *argv[]) {
  // printf("%d %d %d\n", RIGHT, 0x5555, (int) ((short int) RIGHT)); BUG;
  //printf("%ld\n", (Ulong) RF_NA); BUG;
  
  setOptions5(!true, 7, 0, 0, 0, 0,0,0,0,0, 1);

  STARTCLOCK;
  option_type *global;
  utilsoption_type *utils;
  get_started( &global, &utils);
  basic_options *opt = &(utils->basic);
  tuning_options *tuning = &(global->tuning);
 
  get_started( &global, &utils);
 
#define maxarg 19
  const char *argList[maxarg]  =
    {"snps", "indiv", "repetV",
     "cores", "floatLoop", "meanSubst", "missingsFully0", "centered",
     "trials", "variant", "coding", 
     "SxInotVmean", "mode",  "SNPmode", "cmp",
     "sparse", "rowMode", "colMode", "fillMode"};
  Long argdefault[maxarg]=
    { 51200, 102400, 1, 
      4, false, false, true, true,
      TRIALS, 32, AutoCoding, 
      false, -1, -1, 0,
      0, -4, 0, 0};
  input(argc, argv, argList, maxarg, argdefault, true);
  int z = 0;
  Long snps = argdefault[z++];//printf("%s %ld\n",argList[z-1],argdefault[z-1]);
  Long indiv = argdefault[z++];
  Long repetV = argdefault[z++];
  
  opt->cores = (int) argdefault[z++];
  opt->cores = GreaterZero(opt->cores); // immer getrennt!!!!
  tuning->floatLoop = (int) argdefault[z++];
  tuning->meanVsubstract = argdefault[z++];
  tuning->missingsFully0 = argdefault[z++];
  bool centered = argdefault[z++];
  global->genetics.centered = centered ? RowMeans : NoCentering;
  global->genetics.normalized = NoNormalizing;

  int trials = (int) argdefault[z++];
  int variant = (int) argdefault[z++];
  coding_type newcoding = (coding_type) argdefault[z++];
  
  bool SxInotVmean = argdefault[z++];
  tuning->meanSxIsubstract = tuning->meanVsubstract xor SxInotVmean;
  //printf(" tuning->meanSxIsubstract %d = %d xor %d\n",
  //	 tuning->meanSxIsubstract, tuning->meanVsubstract, SxInotVmean);
  int mode = (int) argdefault[z++];
  int SNPmode = (int) argdefault[z++];
  int cmp = (int) argdefault[z++];

  int sparse_nrow = (int) argdefault[z++];
  int rowMode = (int) argdefault[z++];
  int colMode = (int) argdefault[z++];
   int fillMode = (int) argdefault[z++];
 
  opt->efficient = variant >= 512;
  CLOCK("starting");
  
  Long size = indiv * repetV;
  double *V = RandomVector(size, mode),
    *ans = (double*) MALLOC(sizeof(double) * size);
  Long sizeS = snps * repetV;
  double *VS = RandomVector(sizeS, mode);
  CLOCK("creating random vectors");
    

  coding_type true_coding = newcoding <= AutoCoding ? Plink :
    newcoding == UnknownSNPcoding ? Plink : newcoding;
  tuning->variant = variant;
  if (true_coding == OrigPlink && !sparse_nrow) BUG;
  
    
  SEXP SxI = R_NilValue, Coded = R_NilValue;
  Long *info = NULL;
  if (false) {
    PRINTF("snps * indiv = ");
    if ((double) snps * (double) indiv < 1e9) PRINTF("%ld\n", snps * indiv);
    else PRINTF("%1.2f G\n", (double) snps * (double) indiv / 1e9);
  }

  #define cmplimit 15e9
  if ((double) snps * (double) indiv <= cmplimit) { // 15e9
    SxI = CreateRandomSNPmatrix(snps, indiv, /* PROTECT() */
				OneByteGeno, variant, SNPmode,
				global, utils);
    info = GetInfo(SxI);
    CLOCK("creating pseudo-random SNP matrix (byte coding)");

    //printf("true_coding = %d %s\n", true_coding, CODING_NAMES[true_coding]);

    if (true_coding != OrigPlink) {
      Coded = Transform(SxI, NULL, true_coding, NULL, 0, NULL, 0,/* PROTECT() */
			global, utils);
      CLOCK1("transforming into '%s'-coding", CODING_NAMES[true_coding]);
      
      if (newcoding == UnknownSNPcoding) {
	SEXP CodedOrig = Coded;
	true_coding = FiveCodes;
	Coded = Transform(CodedOrig, NULL, true_coding,/* PROTECT() */
			  NULL, 0, NULL, 0,
			  global, utils);
	CLOCK2("transforming '%s' into '%s'-coding",
	       CODING_NAMES[Plink], CODING_NAMES[true_coding]);
	if (false) {
	  SEXP F =  Transform(SxI, NULL, FiveCodes,/* PROTECT() */
			      NULL, 0, NULL, 0,
			      global, utils);
	  CLOCK1("transforming into '%s'-coding", CODING_NAMES[FiveCodes]);
	  Uchar *cF = (Uchar*) Align(F, opt);
	  Uchar *cT = (Uchar*) Align(Coded, opt);
	  for (int i=0; i<5; i++) PRINTF("F=%d T=%d\n", (int) cF[i], cT[i]);
	  FREE_SEXP(&F);
	}
	FREE_SEXP(&CodedOrig);  
      }
    }
  } else {
    if (cmp) PRINTF("\n\n***** NO COMPARISON ***** \nsnps * indiv = %ld > %ld\n\n",  snps * indiv, (Long) cmplimit);
    cmp = false;
    if (true_coding != newcoding && newcoding != AutoCoding) {
      PRINTF("%s %s\n", CODING_NAMES[true_coding], CODING_NAMES[newcoding]);
    }
    assert(true_coding == newcoding || newcoding == AutoCoding);
    Coded = CreateRandomSNPmatrix(snps, indiv,/* PROTECT() */
				  true_coding, variant, SNPmode,
				  global, utils);
     CLOCK("creating true_coding SNP matrix");
    info = GetInfo(Coded);
  }

  Long *infoCoded = true_coding == OrigPlink ? NULL : GetInfo(Coded);
   
  if (sparse_nrow) {
    double *valueB = NULL;
    int *rowIdxB = NULL, *colIdxB = NULL;
    RandomSparse(sparse_nrow, (int) indiv, rowMode, colMode, fillMode,
		 &valueB, &rowIdxB, &colIdxB);
    Uchar* code = NULL, *codeSxI = NULL;
    Long AnsBytes = sizeof(double) * sparse_nrow * snps;
    //printf("AnsBytes = %ld\n", AnsBytes);
    double* AnsSxI = NULL,
      *Ans = (double*) MALLOC(AnsBytes);
    if (true_coding != OrigPlink) code = (Uchar*) Align(Coded, opt);
    else {
      if (SxI == R_NilValue) BUG;
      code = (Uchar*) OneByte2Tplink(SxI, utils);

      
      SEXP PROTECT(CPxI = Transform(SxI, NULL, Plink, NULL, 0, NULL, 0,
				    global, utils));
      Uchar *CP = (Uchar*) Align(CPxI, opt);
      Long *inf = GetInfo(CPxI);
      Long lcp = inf[LDA];

      for (int i=0; i<indiv; i++) {
	bool first = true;
	Uchar *c = code + DIV_GEQ(snps, 4) * i;
	Uchar *cp = CP + lcp * sizeof(unit_t) * i;
	for (int j=0; j<DIV_GEQ(snps, 4); j++) {
	  if (c[j] != cp[j]) {
	    if (first) { printf("\n%d:", i); first = false; }
	    printf(" %d:%u,%u", j, (Uint) c[j], (Uint) cp[j]);
	  }
	}	
	if (!first) {printf("\n"); BUG;}
      }
      CLOCK("transforming to OrigPlink");
      UNPROTECT(1);
    }

    int p_nrow = sparse_nrow < 10 ? sparse_nrow : 10;
    int p_ncol = snps < 20 ? (int) snps : 20;
    int p_indiv = indiv < 10 ? (int) indiv : 10;

#define repetSparse 10
    for (int k=0; k<repetSparse; k++)
      sparseTGeno(code, snps, indiv,
		  true_coding == OrigPlink ? DIV_GEQ(snps,4) : infoCoded[LDA],
		  true_coding == OrigPlink ? OrigPlink
		                   : (coding_type) infoCoded[CODING],
		  valueB, sparse_nrow, rowIdxB, colIdxB, false, global, utils,
		  Ans, sparse_nrow);
    CLOCK1("%d Sparse multiplications", repetSparse);
    //    BUG;

   
    if (cmp && SxI != R_NilValue) {
      Long lda = info[LDA];
      codeSxI = (Uchar*) Align(SxI, opt);
      AnsSxI = (double*) MALLOC(AnsBytes);
      sparseTGeno(codeSxI, snps, indiv, lda, (coding_type) info[CODING],
		  valueB, sparse_nrow, rowIdxB, colIdxB, false, global, utils,
		  AnsSxI, sparse_nrow);
      CLOCK("Sparse multiplication using OneByte");
      
      if (false) {
	for (Long j=0; j<p_ncol; j++) {
	  for (Long i=0; i<p_indiv; i++) {
	    PRINTF("%u\t", (Uint) codeSxI[j + i * lda * sizeof(unit_t)]);
	    }
	  PRINTF("\n");
	}
	PRINTF("\n");
      }
      
      if (false) {
	for (Long j=0; j<p_nrow; j++) {
	PRINTF("%ld ", j);	
	for (Long i=rowIdxB[j]; i<rowIdxB[j+1]; i++) {
	  PRINTF("%d:%3.1f ", colIdxB[i], valueB[i]);
	}
	PRINTF("\n");
	}
	printf("\n");
      }
      
      if (false)
	for (int k=0; k<p_nrow; k++) {
	  for (int i=0; i<p_ncol; i++) {
	    printf("%6.1f", AnsSxI[i * sparse_nrow + k]);
	    if (Ans[i * sparse_nrow + k] ==  AnsSxI[i * sparse_nrow + k])
	      printf("*   ");
	    else printf(":%3.1f ", Ans[i * sparse_nrow + k]);
	    PRINTF("\t");
	  }
	  printf("\n");
	}
      

      for (int k=0; k<sparse_nrow; k++) {
	bool first = true;
	for (int i=0; i<snps; i++) {
	  if (Ans[i * sparse_nrow + k] !=  AnsSxI[i * sparse_nrow + k]) {
	    if (first) { printf("\n%d:", k); first = false; }
	    if (Ans[i * sparse_nrow + k] == 0.0) printf(" %d", i);
	    else printf(" %d:%3.1f,%3.1f", i, Ans[i * sparse_nrow + k], AnsSxI[i * sparse_nrow + k]);
	    BUG;
	  }
	}
      }    	
    }
    

    FREE(valueB);
    FREE(rowIdxB);
    FREE(colIdxB);
    FREE(Ans);
    FREE(AnsSxI);
    if (Coded != SxI) FREE_SEXP(&Coded);  
    FREE_SEXP(&SxI);
    PRINTF("ENDE \n");
     return(0);
  }
 

  if (cmp) {
    unit_t* codeSxI =  Align(SxI, opt);
    unit_t* code = Align(Coded, opt);
    Long lda = info[LDA];
  // PRINTF("cmp: %ld %ld; %s %ld %d\n", (Long) codeSxI, (Long) code, CODING_NAMES[info[CODING]], (Long) codeSxI, lda);
    Long sum = 0, sumI = 0, codeI;
    for (Long j=0; j<snps; j++) {
      for (Long i=0; i<indiv; i++) {
	if (false && j <= 5) {
	  if (i<10 || i >= indiv -16) {
	    if (i % 5 == 0) PRINTF("|");
	    PRINTF("%u\t", (Uint) (((Uchar *) (codeSxI + i * lda))[j]));
	  } else if (i == 20) PRINTF("||");
	}
	if (false && j==6 && mode == 1) { // after the above is printed
	  sumI += ((_1Byte *) (codeSxI + i * lda))[0];//consider FIRST row!
	  codeI += ((_1Byte *) (codeSxI + i * lda))[0] *
	    (Long)(0.1 + POW(3.0, (double) (i % 5)));
	  if (i % 5 == 4 || i == indiv - 1) {
	    sum += sumI;
	    PRINTF("%ld;%ld;%ld  ", sumI, sum,codeI);
	    sumI = 0L;
	    codeI = 0;
	  }
	}
      }
      if (j <= 6 && mode == 1) PRINTF("\n");
    }
    PRINTF("%s %ld %ld\n", CODING_NAMES[infoCoded[CODING]], (Long) Coded, lda);
    for (int j=0; j<5; j++) {
      for (int i=0; i<10; i++) {
      PRINTF("%u\t", (Uint) (((Uchar *) (code + i * lda))[j]));
    }
    PRINTF("\n");
    }
    PRINTF("\n V="); for (int i=0; i<10; i++) PRINTF("%2.4f ", V[i]);
    PRINTF("\n");
    CLOCK("comparing...");
  }



  
  CLOCK1("BEFORE STARTING TRIALS = %d\n", trials); 
  infoCoded[VARIANT] = variant;
  //
  printf("trials = %d %d\n", trials, cmp);
  for (int k=0; k<trials; k++) {
    VectorRelMatrix(Coded, cmp ? SxI : R_NilValue, V, repetV, 
		    cmp,  global, utils, ans);
    if (k > 0) { CLOCK("dito"); continue; }
    if (!true && cmp) Summary(ans, indiv, repetV);
    CLOCK2("calculating V * RelMatrix (%s V%ld)",
	   CODING_NAMES[infoCoded[CODING]], infoCoded[VARIANT]);
    if (k  < trials - 1) PRINTF("\n");
  }
  
  
  CLOCK("saying that the next call may fail");
  // printINFO(Coded);
  if (false && infoCoded[VARIANT] != 256) {
    for (int k=0; k<trials; k++) {
      infoCoded[VARIANT] = 256;
      VectorRelMatrix(Coded, cmp ? SxI : R_NilValue, V, repetV, 
		      cmp, global, utils, ans);
      if (k > 0) { CLOCK("dito"); continue; }
      if (cmp) Summary(ans, indiv, repetV);
      CLOCK2("calculating V * RelMatrix (%s V%ld)",
	     CODING_NAMES[infoCoded[CODING]], infoCoded[VARIANT]);
      if (k  < trials - 1) PRINTF("\n");
    }
  }

  if (!true) {
    infoCoded[VARIANT] = variant;
    for (int k=0; k<trials; k++) {
      VectorRelMatrix(Coded, cmp ? SxI : R_NilValue, V, repetV, 
		      cmp, global, utils, ans);
      if (k > 0) { CLOCK("dito"); continue; }
      if (cmp) Summary(ans, indiv, repetV);
      CLOCK2("calculating V * RelMatrix (%s V%ld)",
	     CODING_NAMES[infoCoded[CODING]], infoCoded[VARIANT]);
    if (k  < trials - 1) PRINTF("\n");
     }
  }
 
  if (cmp) {
    if (cmp % 10 == 2) {
      unit_t *c = Align(SxI, opt);
      if (snps * indiv < 1000000000L && mode ==  1) {
	Summary(V, indiv, repetV);
	Summary(VS, snps, repetV);
	_1Byte *x = (_1Byte *) c;
	Long sum = 0,
	  end = SizeOfInt * info[LDA] * indiv;
	for (Long i=0; i<end ; i++) sum += (Long) x[i];
	PRINTF("\nTOTAL SUM = %lu of %ld x %ld elements\n", sum, snps, indiv);
      }
    
      if (indiv < 2000L) {
	double *A = (double*) MALLOC(sizeof(double) * indiv * indiv);
	crossprod(c, snps, indiv, info[LDA], OneByteGeno, (int) info[VARIANT], 
		  0, 0, 0, opt, A);
	if (snps * indiv < 1000000000L && mode ==  1) {
	  PRINTF("\n");Summary("SUMMARY CP", A, indiv * indiv); PRINTF("\n");
      }
	CLOCK("calculating V * RelMatrix (byte basis & crossprod)");
	FREE(A);
      }
    }
    if ((double) snps * (double) indiv < 1e6) {
      VectorRelMatrix(SxI, SxI, V, repetV, cmp,global, utils, ans);
      Summary(ans, indiv, repetV);
      CLOCK("calculating V * RelMatrix (byte basis -- slow version)");
    }
  }

 
  if (newcoding <= 0) {
    CLOCK("saying that the next part is time consuming (factor 10)");
    FREE_SEXP(&Coded);
    true_coding = FiveCodes;
    Coded = Transform(SxI, NULL, true_coding, NULL, 0, NULL, 0,/* PROTECT() */
		      global, utils);
    infoCoded = GetInfo(Coded);
    CLOCK1("transforming into '%s'-coding", CODING_NAMES[true_coding]);  

    for (int k=0; k<trials; k++) {
      infoCoded[VARIANT] = variant;
      VectorRelMatrix(Coded, cmp ? SxI : R_NilValue, V, repetV, 
		      cmp,  global, utils,ans);
      if (k > 0) { CLOCK("dito"); continue; }
      if (cmp) Summary(ans, indiv, repetV);
      CLOCK1("calculating V * RelMatrix (5codes8 V%ld)", infoCoded[VARIANT]);
      if (k  < trials - 1) PRINTF("\n");
    }
    
  }
  
 /*
  
  VectorCrossMatrix(SxI, VS, repetV,  global, utils, ansS);
  TOTALCLOCK("calculating V * CrossMatrix");
  Summary(ansS, snps, repetV);

  */

  FREE(V);
  FREE(ans);
  if (Coded != SxI) FREE_SEXP(&Coded);  
  FREE_SEXP(&SxI);
  return 0;
}

//#endif
