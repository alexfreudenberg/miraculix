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


#define BitsPerCode 8 
#define MY_VARIANT 32

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "TemplateUint.h"
#include "MX.h"


// LDA in Units!
Long LdaOneByte(Long snps, Long LDAbitalign){
  return LDAalignBitSnps(snps * BitsPerCode);
}


void OneByte2T(unit_t *tmp, Long snps, Long individuals, Long lda,
	       unit_t *Ans, Long ldAns) {
  if (tmp == Ans) BUG;
  _1Byte *M = (_1Byte*) tmp;
  _1Byte *N = (_1Byte*) Ans;
  for (Long i= 0; i<individuals; i++)
    for (Long j=0; j<snps; j++) N[i + j * ldAns] = M[j + i * lda];
}



#define coding_OneByte_end(UINT)					\
  for (Long i=start_individual; i<end_individual; i++) {	\
    UINT *pM = (UINT*) (M + (i - start_individual) * ldM);	\
    _1Byte *pAns = (_1Byte*) (Ans + i * ldAns);				\
    for (Long s=start_snp; s<end_snp; pM++) {			\
      pAns[s++] = (_1Byte) *pM;					\
    }								\
  }								\
}

coding_header(_4Byte,OneByte) {				     
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding_OneByte_end(_4Byte)	

coding_header(_1Byte,OneByte) { // geht schneller nur mit MEMCOPY -- ToDo
  // bei gleichem LD !
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding_OneByte_end(_1Byte)


get_matrix_header(_4Byte,OneByte) {
  for (Long i=0; i<individuals; i++) {
    _1Byte *m = (_1Byte*) (code + i * lda);
    unit_t *a = Ans + i * ldAns;
    for (Long j=0;  j<snps; j++) a[j] = (unit_t) m[j];
  }
}


get_matrix_header(_1Byte,OneByte) {
  Long bytes = MIN(lda, ldAns) * BytesPerUnit;
  if (lda == ldAns) {
    if (code != (unit_t*) Ans) MEMCOPY(Ans, code, (Long) individuals * bytes);
  } else {
    for (Long i=0; i<individuals; i++) {
      unit_t *m = code + i * lda,
	*a = ((unit_t*) Ans) + i * ldAns;
      MEMCOPY(a, m, bytes);
    }
  }
}

 
 sumGeno_start(OneByte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) {
    // printf("sumGenoOneBute i=%d < %d\n", i, individuals);
    _1Byte *pcode = (_1Byte*) (code + i * lda);
    Long sum = 0L;					  
    for (Long j=0; j<snps; j++, pcode++) sum += *pcode;
    if (sums != NULL) sums[i] = sum;
    total += sum;
  }
 // printf("sumGenoOneByte ends %ld\n", total);
  return total;						  
}

#define vectorGeno_startOneByte(TYPE, SNDTYPE, TRDTYPE, ENDTYPE, NAME)	\
  vectorGeno_start(TYPE, TRDTYPE, ENDTYPE, NAME);			\
  Long floatLoop = abs(tuning->floatLoop); /* // OK */			\
  floatLoop = floatLoop == 0 || floatLoop > rows ? rows : floatLoop;	\
  Long rowsMfloatLoop = (rows - 1) / floatLoop * floatLoop ;		\
  if (false) PRINTF("float = %ld fl=%d rMfl=%ld r=%ld %ld,%ld,%ld,%ld\n", floatLoop, tuning->floatLoop, rowsMfloatLoop, rows, sizeof(TYPE), sizeof(SNDTYPE), sizeof(TRDTYPE), sizeof(ENDTYPE)); \
  if (false) { for (Long j=0; j<rows; j++) { PRINTF(">%f ", (double) V[j]); } PRINTF("\n");}					\


// float = 1 1 5119 5120

#define vectorGeno_fctnOneByte(TYPE, SNDTYPE, TRDTYPE, ENDTYPE)		\
  for (Long r=0; r < repetV; r++) {					\
    for (Long i=0; i<cols; i++) {					\
      TYPE *vv = V + r * ldV;						\
      ENDTYPE *a = Ans + r * ldAns;					\
      TRDTYPE total = 0;						\
      unit_t *C = code + i * lda;						\
      if (pseudoFreq == NULL) {						\
	for (Long k=0; k<=rowsMfloatLoop; k+=floatLoop) {		\
	  if (false) PRINTF("nopseudofreq; i=%ld k=%ld < %ld\n", i, k, rowsMfloatLoop);	\
	  Long end_j = k < rowsMfloatLoop ?  k + floatLoop : rows;	\
	  SNDTYPE sum = 0;						\
	  for (Long j=k; j < end_j; j++) {				\
	    sum += (SNDTYPE) ((TRDTYPE) vv[j] * (TRDTYPE)(((_1Byte*) C)[j])); \
	  }								\
	  total += (TRDTYPE) sum;					\
	}								\
      } else {								\
	for (Long k=0; k<=rowsMfloatLoop; k+=floatLoop) {		\
	  Long end_j = k < rowsMfloatLoop ?  k + floatLoop : rows;	\
	  SNDTYPE sum = 0;						\
	  for (Long j=k; j < end_j; j++)				\
	    sum += (SNDTYPE) ((TRDTYPE) vv[j] *				\
			      ((TRDTYPE) (((_1Byte*) C)[j])		\
			       - (TRDTYPE) 2 * (TRDTYPE) pseudoFreq[j])); \
	  total += (TRDTYPE) sum;					\
 	}								\
      }									\
      a[i] = (ENDTYPE) total;						\
    }									\
  }									\
}
 
#define sndtype double
#define trdtype double
//#define sndtype LongDouble
//#define trdtype LongDouble
 vectorGeno_startOneByte(double, sndtype, trdtype, double, OneByte);
 //   printf("YYY XXXXXXXXX\n"); BUG;
#ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif
vectorGeno_fctnOneByte(double, sndtype, trdtype, double)

  

#define sndtypeLD LongDouble
#define trdtypeLD LongDouble
 vectorGeno_startOneByte(longD, sndtypeLD, trdtypeLD, double, OneByte);
 // printf("XXXXXXXXX\n"); BUG;
#ifdef DO_PARALLEL
 //  #pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif
 vectorGeno_fctnOneByte(longD, sndtypeLD, trdtypeLD, double)


 
#define sndtypeD float
#define trdtypeD double
 vectorGeno_startOneByte(floatD, sndtypeD, trdtypeD, double, OneByte);
#ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif
 vectorGeno_fctnOneByte(floatD, sndtypeD, trdtypeD, double)
 
vectorGeno_startOneByte(LongDouble, LongDouble, LongDouble,LongDouble, OneByte);
#ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif
vectorGeno_fctnOneByte(LongDouble, LongDouble, LongDouble, LongDouble)

vectorGeno_startOneByte(Ulong, Ulong, Ulong, Ulong, OneByte);
#ifdef DO_PARALLEL
  #pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif
vectorGeno_fctnOneByte(Ulong, Ulong, Ulong, Ulong)


  
#define genoVectorOneByte_start(TYPE, SNDTYPE, TRDTYPE, ENDTYPE)	\
  genoVector_header(TYPE, TRDTYPE, ENDTYPE, OneByte) {			\
    int cores = GreaterZero(opt->cores);				\
    Long floatLoop = abs(tuning->floatLoop); /* // OK */		\
    floatLoop = floatLoop == 0 || floatLoop > cols ? cols : floatLoop;	\
    Long colsMfloatLoop = (cols - 1) / floatLoop * floatLoop ;		\
    SNDTYPE   *Sum = (SNDTYPE *) CALLOC(rows * repetV, sizeof(SNDTYPE)); \
    TRDTYPE *Total = (TRDTYPE *) CALLOC(rows * repetV, sizeof(TRDTYPE)); \
    if (false) PRINTF("FLOATLOOOP = %ld (%ld, %ld) %lu\n",  floatLoop, sizeof(SNDTYPE), sizeof(TRDTYPE), colsMfloatLoop); \
 
 
#define genoVector_fctnOneByte(TYPE,  SNDTYPE, TRDTYPE, ENDTYPE)	\
  for (Long r=0; r < repetV; r++) {					\
    for (Long k=0; k<=colsMfloatLoop; k+=floatLoop) {			\
      SNDTYPE *sum = Sum + r * rows;			\
      Long end_i = k < colsMfloatLoop ?  k + floatLoop : cols;		\
      for (Long i=k; i < end_i; i++) {					\
	/*  printf("her\n");	*/					\
	TRDTYPE vv = (TRDTYPE) V[r * ldV + i];				\
        TRDTYPE aPF = 0;						\
	if (pseudoFreq != NULL) { 			\
	  if (false) {PRINTF("pseudofreq\n"); }				\
	  aPF = (TRDTYPE)2 * vv * (TRDTYPE) pseudoFreq[i];		\
	  if (false) PRINTF("%e ", (double) aPF);			\
	}								\
	_1Byte *c = (_1Byte*) (code + lda * i);				\
	for (Long j=0; j < rows; j++) {					\
	  /* printf("r=%ld i=%ld l=%ld\n", r, i, j); */			\
	  sum[j] += (SNDTYPE) (vv * (TRDTYPE) c[j] - aPF);		\
  	}								\
      }									\
      TRDTYPE *total = Total + r * rows;		\
      for (Long j=0; j < rows; j++) {					\
        total[j] += (TRDTYPE) sum[j];					\
     }									\
    }									\
    ENDTYPE  *a = Ans + r * ldAns;					\
    TRDTYPE *total = Total + r * rows;					\
    for (Long j=0; j < rows; j++) {					\
      /* printf("r=%ld i=%ld l=%ld\n", r, i, j); */			\
      a[j] = (ENDTYPE) total[j];					\
    }			 						\
  }									\
  FREE(Sum);								\
  FREE(Total)

 


 
genoVectorOneByte_start(longD, sndtypeLD, trdtypeLD, double)
  if (false) PRINTF("double LONG LONG double\n");
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
genoVector_fctnOneByte(longD, sndtypeLD, trdtypeLD, double)
  }
 
genoVectorOneByte_start(double, double,  double, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
genoVector_fctnOneByte(double, double, double, double)
  }

genoVectorOneByte_start(floatD, sndtypeD, trdtypeD, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
genoVector_fctnOneByte(floatD, sndtypeD, trdtypeD, double)
  }

genoVectorOneByte_start(LongDouble, LongDouble, LongDouble, LongDouble)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores)  schedule(static)
#endif
genoVector_fctnOneByte(LongDouble, LongDouble, LongDouble, LongDouble)
}
  
genoVectorOneByte_start(Ulong, Ulong, Ulong, Ulong)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
genoVector_fctnOneByte(Ulong, Ulong, Ulong, Ulong)
}


#define gVOneByte_Select(TYPE)						\
  gV_vG_start(TYPE, OneByte);						\
  if (false) PRINTF("gVOneByte gV=%d %s %d\n", gV, "genoVectorOneByte_"#TYPE, pseudoFreq == NULL); \
  vD = gV ? genoVectorOneByte_##TYPE : vectorGenoOneByte_##TYPE;	\
  gV_vG_end;								\
  if (false) PRINTF("end gV_vG1B longD=%d\n", 1 );			\
}  


gVOneByte_Select(LongDouble)
gVOneByte_Select(Ulong)


gV_vG_start(double, OneByte);
if (false) PRINTF("++++ gVOneByte gV=%d %s %d loop=%d\n", gV, "genoVectorOneByte_double", pseudoFreq == NULL, tuning->floatLoop); 
if (tuning->floatLoop < 0) {						
  vD = gV ? genoVectorOneByte_longD : vectorGenoOneByte_longD;	
 } else if (tuning->floatLoop == 0)					
   vD = gV ? genoVectorOneByte_double : vectorGenoOneByte_double;	
 else vD = gV ? genoVectorOneByte_floatD : vectorGenoOneByte_floatD;	

gV_vG_end;								
  //  printf("end gV_vG1B longD=%d\n", vD == vectorGenoOneByte_longD || vD == genoVectorOneByte_longD );
} 
  


void allele_sum_OneByte(unit_t *code, Long snps, Long individuals, 
			Long lda, Long VARIABLE_IS_NOT_USED ldabitalign,
			basic_options VARIABLE_IS_NOT_USED  *opt,
			Ulong *ans) {
#define x127 127    
  unit_t *sumByte =  (unit_t*) CALLOC(lda, BytesPerUnit);// OK
  Uchar *s = (Uchar*) sumByte;
  Long units = snps2entities(snps, BytesPerUnit);	
  for (Long i=0; i<individuals; ) {
    unit_t *c = code + lda * i;
    for (Long j=0; j<units; j++) sumByte[j] += c[j]; // SIMD-bar
    i++;
    if (i % x127 == 0) { // somit egal ob Haplo oder Geno
      for (Long j=0; j<snps; j++) ans[j] += s[j];
      MEMSET(sumByte, 0, lda * BytesPerUnit); // OK
    }
  }
  if (individuals % x127 != 0) for (Long j=0; j<snps; j++) ans[j] += s[j];
  FREE(sumByte);
}


sparseTGeno_start(OneByte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif   
  for (int r=0; r<nrow_code; r++) {
    double *a = Ans + r * ldAns;
    _1Byte *c = code + r;
    for (int k=0; k<nIdx; k++) {
      const int jEnd = rowIdxB[k+1];
      double tmp = 0;
      for (int j=rowIdxB[k]; j<jEnd; j++)
	tmp += valueB[j] * (double) c[colIdxB[j] * ldaInByte];
      a[k] = tmp;
    }
  }
}
