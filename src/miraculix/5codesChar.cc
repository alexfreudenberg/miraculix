
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

#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "MXinfo.h"
#include "Template.h"
#include "5codes.h"
#include "plink.h"
#include "Vector.matrix.h"
#include "GPUapi.h"
#include "intrinsics_specific.h"
// #include "compa


#define MY_VARIANT 8


const Long CpByte = GetCodesPerUnit(FiveCodes) / BytesPerUnit;  // OK
const Long CpByteplink = GetCodesPerUnit(Plink) / BytesPerUnit; // OK
const Long BitsPerCodeplink = BitsPerByte / CpByteplink;

//#define CpByteplink 4
//#define CpByte 5
//#define BitsPerCodeplink 2


extern Uchar PLINK2FIVE[1024];
extern Uint ALLELESUM[243];
 


#if !defined CUDA
void plink2gpu(char VARIABLE_IS_NOT_USED  *plink, // including three-byte
	       // header, size in bytes: ceiling(indiv/4) * snps + 3
    char VARIABLE_IS_NOT_USED *plink_transposed, // @JV: If I remember correctly, you need to
                            // transpose the matrix anyway? so it should be
                            // easier to have this as an argument, and not
                            // transpose it again on the GPU
    int VARIABLE_IS_NOT_USED snps,
    int VARIABLE_IS_NOT_USED indiv, // actual number of individuals, not N/4 rounded up
	       double VARIABLE_IS_NOT_USED *f,  // vector guaranteed to be of length  dsnps
    int VARIABLE_IS_NOT_USED n, 
    void VARIABLE_IS_NOT_USED **GPU_obj) {
  BUG;
    }
void dgemm_compressed_gpu(
    bool VARIABLE_IS_NOT_USED trans, //
    void VARIABLE_IS_NOT_USED *GPU_obj,
    //    double *f,
    int VARIABLE_IS_NOT_USED n,     // number of columns of matrix B
    double VARIABLE_IS_NOT_USED *B,  // matrix of dimensions (nsnps, n) if trans =
                // "N" and (nindiv,n) if trans = "T"
    int VARIABLE_IS_NOT_USED ldb,
    int VARIABLE_IS_NOT_USED centered,
    int VARIABLE_IS_NOT_USED normalized,
    double VARIABLE_IS_NOT_USED *C,
    int VARIABLE_IS_NOT_USED ldc
			  ) {
  BUG;
}
void freegpu(void VARIABLE_IS_NOT_USED **GPU_obj){
  BUG;
}
#endif


void print_options(option_type *global, utilsoption_type *utils) {
  if (utils->basic.Cprintlevel) {
    if (false) {
      PRINTF("\ndebug info: addresses %ld %ld", (Ulong) global, (Ulong) utils);
    }
    const char *ft [] ={"false", "true", "tba", "tba"};
    int G = (int) global->tuning.gpu << 1;
    PRINTF("\n");     
    PRINTF("gpu      : %s\n", ft[global->tuning.gpu]);
    if (!global->tuning.gpu) {
      if ( global->genetics.centered != NoCentering)
	PRINTF("ignore missings: %s\n", ft[G + global->tuning.missingsFully0]);
      PRINTF("do center: %s\n", ft[G + (global->genetics.centered !=NoCentering)]);
      if (global->tuning.floatLoop >=0)
	PRINTF("float    : %s\n", ft[G + global->tuning.floatLoop > 0]);
      else PRINTF("precise  : true\n");
      PRINTF("meanV    : %s\n", ft[G + global->tuning.meanVsubstract]);
      PRINTF("meanSxI  : %s\n", ft[G + global->tuning.meanSxIsubstract]);
      PRINTF("normalize: %s\n", ft[G + global->genetics.normalized] );
      PRINTF("variant  : %d \n",   global->tuning.variant);
    }
    PRINTF("use fortran freq: %s\n", 
	   ft[G + global->genetics.prefer_external_freq] );
    PRINTF("cores    : %d \n",  utils->basic.cores );
    PRINTF("print lvl: %d \n",  utils->basic.Cprintlevel );
  }
}



void get_started(option_type **g, utilsoption_type **u) {
  static option_type *global = NULL;
  static utilsoption_type *utils = NULL;
  if (global == NULL) {
    startMiraculix(1);
    WhichOptionList(false, &global, &utils);
    
    //    printf(">>>>>>>>> global %ld\n", (Long) global);

    //printf("GPU setting options in get_started\n");
    setOptions5(
#if defined CUDA
		true,
#else
		false,
#endif
		0, // cores 
		false, // floatLoop
		false, // meanV
		false, // meanSxI
		false, // irgnore missings
		RowMeans, // centering
		NoNormalizing, // normalize
		true, // use Fortran freq
		128,
		false // print details
		);
    Ext_startRFU();
    //   printf(" >> "); //, (Long) (&OPTIONS), (Long) utils;
  }
  *g = global;
  *u = utils;
 assert(global != NULL);
  assert(utils != NULL);
  //  printf("getstarted_done %ld\n", (Long) global);
  if (utils->basic.Cprintlevel > 10) {
    PRINTF("<<<<< global %ld\n", (Long) global);
    print_options(global, utils);
  }
}

void getStartedOptions() {
  static option_type  VARIABLE_IS_NOT_USED *global = NULL;
  static utilsoption_type VARIABLE_IS_NOT_USED  *utils = NULL;
  void get_started(option_type **g, utilsoption_type **u);
}


void setOptions5(int gpu, int cores, int floatLoop,
		 int meanVsubstract, int meanSxIsubstract,
		 int missingsFully0,
		 int centered, int normalized,
		 int prefer_external_freq,
		 int variant, int print_details) {
  if (print_details > 10) PRINTF("entering setoptions cores = %d\n", cores);
  option_type *global;
  utilsoption_type *utils;   
  get_started(&global, &utils);
   tuning_options *tuning = &(global->tuning);
  genetics_options *genetics = &(global->genetics);
  if (cores == 0) {
    const char* s = getenv("OMP_NUM_THREADS");
    if (s != NULL) cores = (int) strtoimax(s, NULL, 10);
    if (cores <= 0) {
      if (print_details > 3)
	PRINTF("WARNING: seeing 'OMP_NUM_THREADS=%s' only -- #threads set to 4\n", s);
       cores = 4;
       //PRINTF("number of threads not positive");  BUG;
    }
    //   FREE(s);
  }
  utils->basic.cores = cores;
  tuning->gpu = (bool) gpu;
  //printf("GPU set to %d\n", tuning->gpu);
  //  gpu: centered works
  if (gpu && (!prefer_external_freq || missingsFully0 || normalized))
    ERR0("in case of 'gpu' the Fortran frequency must always be used; missings/centering/normalizing cannot treated.");
  tuning->floatLoop = floatLoop;
  tuning->meanVsubstract = (bool) meanVsubstract;
  tuning->meanSxIsubstract = (bool) meanSxIsubstract;
  tuning->missingsFully0 = (bool) missingsFully0;
  genetics->centered = (centering_type) centered;
  genetics->prefer_external_freq = prefer_external_freq;
  genetics->normalized = (normalizing_type) normalized;
  tuning->variant = variant < 32 ? 32 : variant > 256 ? 256 : variant;
  utils->basic.efficient = false;
  utils->basic.Cprintlevel = print_details * 3;
  tuning->addtransposed = true;
  if (print_details > 10) {
    print_options(global, utils);
    PRINTF("setoptions5 done \n");
  }
}


#define TBM 0x03
void plink2Geno5codestrans32(Uchar *M, Long snps, Long individuals,
			     int VARIABLE_IS_NOT_USED cores,
			     unit_t *Ans, Long ldAns) { 
  // printf("p25T enter\n");
  if (sizeof(Uchar) * BitsPerByte != MY_VARIANT) BUG;
  Uchar *table = PLINK2FIVE;
  
  Long
    ldM = 1L + (individuals - 1L) / CpByteplink, // in bytes!!!
    end_snpsM1_CpB = (snps - 1L) / CpByte;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<=end_snpsM1_CpB; i++) {
    int rest = (int) CpByte;
    Uchar *pM =  M + ldM * CpByte * i;
    Uchar *a0 = ((Uchar*) (Ans + ldAns * i)); 
    Uchar h[BitsPerByte]; // CpByte <= BitsPerByte !
    
    if (i < end_snpsM1_CpB) {} else {
      rest = (int) (snps - end_snpsM1_CpB * rest);
      assert(rest > 0);
      for (int ll=rest; ll<CpByte; ll++) h[ll] = 0;
    }
    
    for (Long j=0; j<ldM; j++, pM++) {
      Uchar *a = a0 + j * CpByteplink;
      for (int z=0; z<rest; z++) h[z] = pM[z * ldM];
      for (int k=0; k<CpByteplink; k++) {
	const Uint m = ((Uint) (h[0] & TBM)) |
	  (((Uint) (h[1] & TBM)) << 2) |
	  (((Uint) (h[2] & TBM)) << 4) |
	  (((Uint) (h[3] & TBM)) << 6) |
	  (((Uint) (h[4] & TBM)) << 8);
	a[k] = table[m];
	//	printf("%d ", (int) a[k]);
	for (Long ll=0; ll<CpByte; ll++) h[ll] >>= BitsPerCodeplink;
      }
    }
  }
  //  printf("p25T done\n"); BUG;
}




#if defined LONG_HANGING_TYPE
#undef LONG_HANGING_TYPE
#endif
//#define LONG_HANGING_TYPE Long
//#define HTT_BITS 64
#define LONG_HANGING_TYPE Ulong
#define HTT_BITS 64
#define ORIGBITS_TRANS_MASK 0x00000000000003FFL
// irgendwo fehler drin -- sollte eigentlich bisschen schneller gehen

void plink2Geno5codes32(Uchar *M, Long snps, Long individuals,
			int VARIABLE_IS_NOT_USED cores,
			unit_t *Ans, Long ldAns) {
  //  printf("p25 enter\n");
  Uchar *table = PLINK2FIVE;
  const Long ldM = 1L + (individuals - 1L) / CpByteplink; // in bytes!!!
  const Long end_indiM1_CpB = (individuals -1L) / CpByte;
  const Long ldaByte = ldAns * BytesPerUnit; // OK
  const Long totalBytes = ldM * snps;

  if (MY_LDABITALIGN_2BIT < sizeof(Long)) BUG;
  if (HTT_BITS != sizeof(LONG_HANGING_TYPE) * BitsPerByte) BUG;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<snps; i++) {
    //    printf("i=%ld %d\n", i, snps);
    int availableBits = 0;
    LONG_HANGING_TYPE hanging = 0;
    LONG_HANGING_TYPE m = 0;
    int hangingBits = 0;
    Uchar *a0 = ((Uchar*) Ans) + i; // !!
    Long byteNr = i * ldM;
    LONG_HANGING_TYPE *pM = (LONG_HANGING_TYPE*) (M + byteNr);
    for(Long j=0; j<=end_indiM1_CpB; j++) {      
      //  printf("i=%ld %d j=%ld %d ldM=%d %d\n", i, snps, j,end_indiM1_CpB, ldM, individuals);
      // i=127 128 j=224 225 ldM=282 1128
      Uchar *a = a0 + ldaByte * j;					
      while (availableBits < ORIGBITSperFIVE) {
	if (hangingBits == 0) {
	  byteNr += sizeof(LONG_HANGING_TYPE); // endpoint after loading
	  // printf("hanging %d <= %d \n", byteNr,totalBytes );
	  if (byteNr <= totalBytes) {
	    // printf(".");
	    // inbetween just (non-annoying) nonsense is additionally loaded
	    // if ldM is passed by *pM
	    hanging = *(pM++);
	    hangingBits = HTT_BITS;
	  } else {
	    // at the very end, instead of nonsense loading, we would have
	    // a segmenation fault
	    byteNr -= sizeof(LONG_HANGING_TYPE);
	    int left =  (int) (totalBytes - byteNr);
	    for (int k=0; k < left; k++) {
	      //     printf("k=%d left = %d %d\n", k, left, (int) (((Uchar*) pM)[k]));
	      //	  printf("k=%d left = %d %d\n", k, left, (int) ( ((Uchar*) (&hanging))[k]));

	      ((Uchar*) (&hanging))[k] = ((Uchar*) pM)[k];
	      
	    }
	    hangingBits = left * BitsPerByte;
	  }
	}
	//	printf("%d < %d\n", availableBits ,  ORIGBITSperFIVE);
	m |= hanging << availableBits;
	int movedBits = hangingBits <= HTT_BITS-availableBits
	  ? hangingBits : HTT_BITS - availableBits;
	//	printf("i=%ld j=%ld %d %ld %d  %d<%d-%d\n", i, j, hangingBits, (Long) M - (Long) pM, movedBits, hangingBits, sizeof(LONG_HANGING_TYPE), availableBits);
	availableBits += movedBits;
	hangingBits -= movedBits;
	hanging >>= movedBits;
      }
      *a = table[m & ORIGBITS_TRANS_MASK];
       m >>= ORIGBITSperFIVE;
       availableBits -= ORIGBITSperFIVE;
    }
  }
  //  printf("p25 done\n");
  // BUG;
}


		   
void plinkToSxI(char *plink, char *plink_tr, 
	      Long snps, Long indiv,
	      int Coding,
	      double *f, 
	      int max_n,
	      void**compressed) {
  assert(CpByte == 5);
  assert(CpByteplink == 4);
  assert(BitsPerCodeplink == 2);

  option_type *global;
  utilsoption_type *utils;
  get_started(&global, &utils);
  
  stopIfNotInt(snps | indiv);
  if (global->tuning.gpu) {
    plink2gpu((char*) plink, (char*) plink_tr, (int) snps, (int) indiv, f,
	      max_n, compressed);
    return;
  } 

  coding_type coding = (coding_type) Coding;
  basic_options *opt = &(utils->basic);

  SEXP SxI = PROTECT(CompleteCodeVector(snps, indiv, coding,
					global->tuning.variant, 0, false,
					global, utils,
					coding == OrigPlink ? (void*) plink_tr
					: NULL,
					coding == OrigPlink ? plink : NULL,
					NULL));
  if (coding == FiveCodes) {  
    *compressed = (void*) SxI;
    int cores = GreaterZero(opt->cores);
    unit_t *code = Align(SxI, opt);
    Long *info = GetInfo(SxI);
    Long ldAns = info[LDA]; 
    
    if (global->tuning.missingsFully0) getPlinkMissings((Uchar*) plink, SxI,
							opt);
    //  printf("coding ...\n");
    plink2Geno5codes32((Uchar*) plink, snps, indiv, cores, code, ldAns);
    
    SEXP next = getAttribPointer(SxI, Next);
    unit_t *codenext = Align(next, opt);
    Long *infonext = GetInfo(next),
      ldanext = infonext[LDA];
    
    //printf("trans coding ...\n");
    plink2Geno5codestrans32((Uchar*)plink, snps, indiv, cores, codenext, ldanext);
  }
  if (global->genetics.prefer_external_freq) {
    getFreq(SxI, f, global, utils);
  }
  UNPROTECT(1);
}


void vectorGeno5api(void *compressed, double *V, Long repetV, Long LdV, 
		    double *ans, Long LdAns) {
  option_type *global;
  utilsoption_type *utils;
  get_started(&global, &utils);
 if (global->tuning.gpu) {
    //
   stopIfNotInt(repetV | LdV | LdAns);
   dgemm_compressed_gpu(false, compressed, (int) repetV, V, (int) LdV,
			 global->genetics.centered, 
			 global->genetics.normalized, 
			ans, (int) LdAns);
 } else {
    vectorGeno_means_double((SEXP)compressed, V, repetV, LdV, global, utils,
			    ans, LdAns);
  }
}


void genoVector5api(void *compressed, double *V, Long repetV, Long LdV, 
		 double *ans, Long LdAns) {
  option_type *global;
  utilsoption_type *utils;
  get_started(&global, &utils);
  if (global->tuning.gpu) {
    stopIfNotInt(repetV | LdV | LdAns)
      dgemm_compressed_gpu(true, compressed, (int) repetV, V, (int) LdV, 
			 global->genetics.centered, 
			 global->genetics.normalized,
			   ans, (int) LdAns);
  } else {
    genoVector_means_double((SEXP)compressed, V, repetV, LdV, global, utils,
			    ans, LdAns);
  }  
} 



void free5(void **compressed) {
  option_type *global;
  utilsoption_type *utils;   
  get_started(&global, &utils);
  if (global->tuning.gpu) {
    freegpu(compressed);
  } else {
    FREE_SEXP((SEXP*) compressed); // will crash?!!!
  }
}


void getFreq5(void *compressed,  double *f) {
  SEXP SxI = (SEXP) compressed;
  //   printINFO(SxI);
  option_type *global;
  utilsoption_type *utils;
  get_started(&global, &utils);
  Long *info = GetInfo(SxI);
  Long
    snps = info[SNPS],
    rows = snps;
  //  Long individuals = info[INDIVIDUALS];
  getFreq(SxI, global, utils);
  LongDouble *freq = getFreq(SxI);
   for (Long j=0; j<rows; j++) {
    f[j] = (double) freq[j];
    //   printf("j=%d %f\n", j, f[j]);
  }
}


void sparseTGenoPlinkApi(char *compressed,
			 int snps, 
			 int indiv,
			 double *B,
			 int nIdx,// 2nd dimension of sparse; == length(rowIdxB)
			 int *rowIdxB,
			 int *colIdxB,
			 int tSparse,
			 double *C,
			 int Ldc) {
  option_type *global;
  utilsoption_type *utils;
  get_started(&global, &utils);
  int CodesPerByte = (int) (BitsPerByte / BitsPerCodeplink);
  sparseTGeno((unsigned char*) compressed, snps, indiv,
	      DIV_GEQ(snps, CodesPerByte), 
	      OrigPlink, 
	      B, nIdx, rowIdxB, colIdxB, (bool) tSparse, global, utils,
	      C, Ldc);
}


 
void vectorGenoPlinkApi(char *compressed,
			 int snps, 
			 int indiv,
			 double *f,
			 double *B,
			 int n,
			 int LdB,
			 double *C,
			 int Ldc) {
  option_type *global;
  utilsoption_type *utils;
  get_started(&global, &utils);
  basic_options *opt = &(utils->basic);
  tuning_options *tuning = &(global->tuning);
  int CodesPerByte = (int) (BitsPerByte / BitsPerCodeplink);
  if (indiv % 32 != 0) ERR0("currently indiv must be a multiple of 32");

  if (f != NULL) {
    BUG;
    //    plinkToSxI + allg VectorGeno-Mechanismus
  } else {  
    vectorGenoPlinkMatrix256_double((unit_t*) compressed, snps, indiv,
				    DIV_GEQ(snps, CodesPerByte),
				    OrigPlink, 256, NULL,
				    B, n, LdB, 
				    opt, tuning,
				    C, Ldc);
  }
}
