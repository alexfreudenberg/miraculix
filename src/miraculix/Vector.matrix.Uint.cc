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
#include "compatibility.SEXP.h"
#include "xport_import.h"
#include "miraculix.h"
#include "MXinfo.h"
#include "options.h"
#include "kleinkram.h"
#include "haplogeno.h"
#include "Template.h"
#include "Files.h"
#include "2bit.h"
#include "5codesDef.h"
#include "5codes.h"
#include "plink.h"
#include "Vector.matrix.h"
#include "OneByte.h"


#define MULTIPLY				\
  for (int j=0; j<c; m += r) {			\
  double sum;					\
  TYPE_INDEP_SCALAR_PROD(v, m, r, sum);		\
  a[j++] = sum;					\
}

/*

SEXP vector012matrixXX(SEXP vector, SEXP matrix) {
  KEY_type *KT = KEYT_M();
  int cores = GreaterZero(KT->global_utils.basic.cores;
  int err = NOERROR,
    n = LENGTH(vector),
    r = nr ows(matrix), // achtung nur int
    c = nc ols(matrix)
    ; 
  if (r != n) ERR0("vector and matrix do not match");
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, c));
  double *m = REAL(matrix),
    *a = REAL(Ans);


#define FOR					\
  for (int i=0; i<n; i++)			\
    if (v[i] == 1) idx1[n1++] = i;		\
    else if (v[i] == 2) idx2[n2++] = i;		\
  //  else if (v[i] != 0) {err=1; goto ErrorHandling;}
      
  if (c < 999999) {
    switch (TYPEOF(vector)) {
    case REALSXP : { xA(REAL(vector), m, r, c, a, cores); } break;
    case INTSXP : { int *v = INTEGER(vector); MULTIPLY } break;
    case LGLSXP : { int *v = LOGICAL(vector); MULTIPLY } break;
    default : ERR0("vector type incompatible");  
    } 
  } else {
    int
      n1 = 0,
      n2 = 0;
    Uint
      *idx1 = (Uint*) MALLOC(siz eof(in t) * n),
      *idx2 = (Uint *) MALLOC(size of(in t) * n);
    switch (TYPEOF(vector)) {
    case REALSXP : { double *v = REAL(vector); FOR } break;
    case INTSXP  : { int *v = INTEGER(vector); FOR } break;
    case LGLSXP : {
      int *v = LOGICAL(vector);
      for (int i=0; i<n; i++)
	if (v[i] == 1) idx1[n1++] = i;
	else if (v[i] != 0) { err = 1; goto ErrorHandling; }
    }
      break;
    default : err = 2; goto ErrorHandling;
    }
    
    // printf("    martin todo : parallelisieren");
    for (Long j=0; j<c; j++) {  
    double *mm = m + j * r,
      sum = 0.0;
      for (int i=0; i<n2; sum += mm[idx2[i++]]);
      sum *= 2;
      for (int i=0; i<n1; sum += mm[idx1[i++]]);
      a[j] = sum;
    }
    
 ErrorHandling:
    FREE(idx1);
    FREE(idx2); 
    
    if (err == 1) { ERR0("only 012 allowed for the vector") }
    else if (err == 2) { ERR0("unknown type of the vector") }
  }
  UNPROTECT(1);
  return Ans;
}
*/


  
void matrixvector012I(double * m, int r, int c,  SEXP vector, double *a) {
  double *b = (double *) CALLOC(r, sizeof(double));
  // Martin : todo : replace by SIMD
#define MULTIPLYv							\
  for (int j=0; j<c; j++) {						\
    if (v[j] == 1) for (int i=0; i < r;  m++) a[i++] += *m;		\
    else if (v[j] == 2) for (int i=0; i < r; m++) b[i++] += *m;		\
    else m += r;							\
  }

  for (int j=0; j<r; a[j++] = 0);
  switch (TYPEOF(vector)) {
  case REALSXP : { double *v = REAL(vector);MULTIPLYv } break;
  case INTSXP : { int *v = INTEGER(vector); MULTIPLYv } break;
  case LGLSXP : { int *v = LOGICAL(vector); MULTIPLYv } break;
  default : ERR0("vector type incompatible");  
  }

  for (int i=0; i < r; i++) a[i] += 2.0 * b[i];
  FREE(b);
}  


 
void vector012matrixI(SEXP vector, double * m, int r, int c,  double *a) {

#define FOR2					\
  for (int i=0; i<r; i++)			\
    if (v[i] == 1) idx1[n1++] = i;		\
    else if (v[i] == 2) idx2[n2++] = i;		\
    else if (v[i] != 0) {err=1; goto ErrorHandling;}	\
  int old = idx1[0];					\
  for (int i=1; i<n1; ) {					\
    int dummy = idx1[i];					\
    idx1[i++] -= old;						\
    old = dummy;						\
  }							\
  old = idx2[0];					\
  for (int i=1; i<n2; ) {				\
     int dummy = idx2[i];				\
    idx2[i++] -= old;					\
     old = dummy;					\
  }							\

      
  int err = NOERROR;
  if (c < 9) {
    switch (TYPEOF(vector)) {
    case REALSXP : { double *v = REAL(vector);MULTIPLY } break;
    case INTSXP : { int *v = INTEGER(vector); MULTIPLY } break;
    case LGLSXP : { int *v = LOGICAL(vector); MULTIPLY } break;
    default : ERR0("vector type incompatible");  
    } 
  } else {
    int
      n1 = 0,
      n2 = 0;
    int
      *idx1 = (int*) MALLOC(sizeof(int) * r),
      *idx2 = (int *) MALLOC(sizeof(int) * r);
    switch (TYPEOF(vector)) {
    case REALSXP : { double *v = REAL(vector); FOR2 } break;
    case INTSXP  : { int *v = INTEGER(vector); FOR2 } break;
    case LGLSXP : {
      int *v = LOGICAL(vector);
      for (int i=0; i<r; i++)
	if (v[i] == 1) idx1[n1] = n1 > 0 ? i - idx1[n1-1] : i;
	else if (v[i] != 0) { err = 1; goto ErrorHandling; }
    }
      break;
    default : err = 2; goto ErrorHandling;
    }
    
    // martin todo : SIMD / parallel
    for (Long j=0; j<c; j++) {  
      double *mm = m + j * r,
	sum = 0.0;
      for (int i=0; i<n2; ) {
	//	print("%d 2: %d %d %d \n", j, i, idx1[i], mm - m);
	mm += idx2[i++];
	sum += *mm;
      }
      sum *= 2.0;
      mm = m + j * r;
      for (int i=0; i<n1; ) {
	//	print("%d 1: %d %d %d \n", j, i, idx1[i], mm - m);
	mm += idx1[i++];
	sum += *mm;
      }
      a[j] = sum;
    }
    
 ErrorHandling:
    FREE(idx1);
    FREE(idx2); 
    
    if (err == 1) { ERR0("only 012 allowed for the vector") }
    else if (err == 2) { ERR0("unknown type of the vector") }
  }
}


#define gV_vG_0(DOUBLE)					\
  gV_vG_header(DOUBLE,) {				\
 if (SubstractMeanSxI) getFreq(SxI, global, utils);	\
     MEMSET(ans, 0, ldAns * repetV * sizeof(*ans)); 	\
     Long *info = GetInfo(SxI, true);			\
    coding_type coding = (coding_type) info[CODING];	\
    if (false) PRINTF("gV_vG memset size=%ld ldAns=%ld %ld %ld %s\n", sizeof(DOUBLE), ldAns, repetV, (Long) ans, CODING_NAMES[coding]); \
    gV_vG_##DOUBLE##_t vD = NULL;					\
    /* printf("gvvg0 :%sv%d %d bytes=%ld %s\n", CODING_NAMES[coding], (int) info[VARIANT], orig_gV, sizeof(DOUBLE), #DOUBLE); */ \
    switch(coding) {					\
    case OneByteGeno : vD = gV_vGOneByte_##DOUBLE; break;	\
    case Plink :						\
    case PlinkTransposed :					\
    case OrigPlink :						\
    case OrigPlinkTransposed :					\
    case TwoBitGeno :						\
    case TwoBitGenoTransposed :					\
      vD = gV_vG2_##DOUBLE; break;				\
    case FiveCodesTransposed:					\
    case FiveCodes :							\
      if (false) PRINTF("gvvg 5\n");					\
      vD = gV_vG5_##DOUBLE; break;					\
    default:								\
      PRINTF("not coded: %s\n", CODING_NAMES[coding]);			\
      BUG;								\
    }									\
    vD(SxI, V, repetV, ldV, orig_gV, SubstractMeanSxI, global, utils,	\
       ans, ldAns);							\
    /*  printf("gvvg0 done\n");		*/				\
  }
gV_vG_0(double)
gV_vG_0(LongDouble)
gV_vG_0(Ulong)


#define gV_raw(DOUBLE) \
  void genoVector_raw_##DOUBLE(SEXP SxI, DOUBLE *V, Long repetV, Long ldV, \
			       option_type *global, utilsoption_type *utils, \
			       DOUBLE *ans, Long ldAns) {		\
    gV_vG_##DOUBLE(SxI, V, repetV, ldV, true, false, global, utils,ans,ldAns); \
  }									\
									\
 void vectorGeno_raw_##DOUBLE(SEXP SxI, DOUBLE *V, Long repetV, Long ldV, \
			      option_type *global, utilsoption_type *utils, \
			      DOUBLE *ans, Long ldAns) {		\
   gV_vG_##DOUBLE(SxI, V, repetV, ldV, false, false, global, utils,ans,ldAns); \
 }
gV_raw(double)
gV_raw(LongDouble)





#define Tolerance 1e-10 
#define Tolerance1 1e-15
#define Tolerance2 1e-10
#define PRECISION_TYPE LongDouble
void VectorRelMatrix( SEXP SxI, SEXP SxI2, double *V, Long repetV, int compare,
		     option_type *global, utilsoption_type *utils,
		     double *ans) {
  const int debug =  false;

  //  if (false) {   printINFO(SxI); BUG; }
 
  
  //  printf("XA %d\n", SxI == SxI2);
  //  v * M^T * M  for  M = snp x indiv  
  //
  Long *info = GetInfo(SxI);
  //  printf("variant=%d\n", (int) info[VARIANT]); assert(info[VARIANT] == 32);
  compare *= (SxI2 != R_NilValue);
  bool sameInter = compare < 0;
  compare = abs(compare); // OK
  int maxLines =          compare > 10 ? 0 : compare;
  bool stopOnError =      compare % 10 < 5;
  int default_floatLoop = compare % 2 == 0 ? 100000000L  // loop in SxI2
					 : compare % 3 == 0
					 ? -100000000L//additionally:LongDouble!
					 : 0;
 
    
  if (SxI2 == R_NilValue) {
    if (debug > 1) PRINTF("entering VRM %s \n", CODING_NAMES[info[CODING]]);
  } else {
    if (debug > 1) PRINTF("entering VRM %s %s\n", CODING_NAMES[info[CODING]],
		      CODING_NAMES[GetInfo(SxI2)[CODING]]);
  }
  
  Long
    snps = info[SNPS],
    individuals = info[INDIVIDUALS];

  double *inter = (double*) MALLOC(sizeof(double) * snps * repetV);
  //
  
  double *inter2 = NULL;
  double sum=0.0, sum1=0.0, sum2=0.0;
  bool meanSxI = global->tuning.meanSxIsubstract,
    meanV = global->tuning.meanVsubstract;
  int floatLoop = global->tuning.floatLoop;

  LongDouble * V_LD = NULL;
  LongDouble *inter2LD = NULL;
  LongDouble *A_LD = NULL;

  //  printf(" **** genoVector *** \n");
  genoVector_means_double(SxI, V, repetV, global, utils, inter);
    
  
  if (compare) {
      
    if (debug) PRINTF("COMPARING; for indiv x repetV = %ld x %ld; loop=%d sameInter=%d\n", individuals, repetV, default_floatLoop, sameInter);

    global->tuning.meanSxIsubstract = false;
    global->tuning.meanVsubstract = false;
    global->tuning.floatLoop = default_floatLoop;
    inter2 = (double*) MALLOC(sizeof(*inter2) * snps * repetV);


    //printf("default_floatLoop=%d\n", default_floatLoop);
    
    if (default_floatLoop < 0 ) {
      V_LD  = (LongDouble*) MALLOC(sizeof(LongDouble) *
					       individuals * repetV);
      for (Long i=0; i<individuals * repetV; i++) V_LD[i] = (LongDouble) V[i];
      inter2LD = (LongDouble*) MALLOC(sizeof(LongDouble) * snps * repetV);
      genoVector_means_LongDouble(SxI2, V_LD, repetV, global, utils, inter2LD);
      Long totalS = snps * repetV;
      for (Long i=0; i<totalS; i++) inter2[i] = (double) inter2LD[i];
    } else {
      genoVector_means_double(SxI2, V, repetV, global, utils, inter2);
    }
  
    global->tuning.meanSxIsubstract = meanSxI;
    global->tuning.meanVsubstract = meanV;
    global->tuning.floatLoop = floatLoop;
   
    //PRINTF("\nXXX snps * repetV=%d; inter=%ld %f %f  inter2=%f %f cores=%d %s,%s\n", snps * repetV, (Long) inter, inter[0], inter[1], inter2[0], inter2[1], utils->basic.cores, CODING_NAMES[info[CODING]], CODING_NAMES[GetInfo(SxI2)[CODING]]); // 512
      
    if (debug) PRINTF("\nsnps * repetV=%ld; inter=%ld %f %f  inter2=%f %f cores=%d %s,%s\n", snps * repetV, (Long) inter, inter[0], inter[1], inter2[0], inter2[1], utils->basic.cores, CODING_NAMES[info[CODING]], CODING_NAMES[GetInfo(SxI2)[CODING]]); // 512
  
    int zaehler = 0;
    double diff = 0.0, maxdiff = 0.0;
    Long limit = snps * repetV < 50 ? 100000 : maxLines;
    
    for (Long r=0; r<repetV; r++) {
      for (Long ii=0; ii<snps; ii++) {
	Long i = r * snps + ii;
	double difference =  inter[i] - inter2[i];
	maxdiff = MAX(maxdiff, FABS(difference));
	if (i < 5||
	    difference != 0.0 || (i == snps * repetV - 1 && limit > 0)
	  || FABS(inter[i]) > 1e6) {
	  //
	  bool special = false;
	  if (false) special = FABS(difference) > 1e4;
	  if (maxLines > 0) {
	    if (limit > 0 || special) {
	      if (debug) {
		if (repetV==1) PRINTF("i="); else PRINTF("%2ld,", r);
		if (!true)
		  PRINTF("%3ld\t%7.3f\t%7.3ft\td=%7.3e\t%7.3e\n",
			 ii, inter[i], inter2[i],  difference,  difference/inter2[i]);
		else
		  PRINTF("%3ld\t%7.3e\t%7.3et\td=%7.3e\t%7.3e\n",
			 i, inter[i], inter2[i], difference,difference/inter2[i]);
	      }
	      special = false;
	    } else if (limit > -5 && debug) PRINTF(".");	    
	    limit--;
	  } else {
	    zaehler++;
	  }
	  if (false && FABS(difference) > 2e-13)
	    PRINTF("I=%ld\t%7.3f\t%7.3f\td=%e\t%e\n",
		   i, inter[i], inter2[i], difference, difference/inter2[i]);
	  diff += FABS(difference);
	}
	sum += inter[i];
	sum1 += FABS(inter[i]);
	sum2 += inter2[i];
      }
    }
    if (maxLines == 0 && debug)
      PRINTF("NOTE! %d positions with average %e difference\n",
	     zaehler, diff / (double) zaehler);
    if (debug) PRINTF("\nsum=%f sum2=%f diff=%5.1e last=%f / %f  max.diff=%e ",
		      sum, sum2,
		      diff, inter[snps * repetV - 1], inter[snps * repetV - 1],
		      maxdiff);
    if (sum != sum2 || diff > sum1 * Tolerance1) {
      PRINTF("\n**** gV sums differ by %e  \n\n",sum - sum2);
	if (FABS(sum - sum2) > FABS(sum) * Tolerance &&
	    FABS(sum - sum2) > Tolerance2)
	  PRINTF("***** %e > max(%e,%e)\n\n", sum - sum2,
		 FABS(sum) * Tolerance, Tolerance2);
 	if (diff > sum1 * Tolerance1)
	  PRINTF("****diff: %e > %e\n", diff, sum1 * Tolerance1);
	if (stopOnError) BUG;
    }
   
    if (debug) PRINTF("\n------------------------------------------------------\n");
  }
  // if (debug) PRINTF("inter[0] %f %f\n", (double) inter2[0], (double) inter[0]); 
  
  
  //printf(" **** vectorGeno *** \n");
  vectorGeno_means_double(SxI, sameInter ? inter2 : inter, repetV,
			  global, utils, ans);

  if (compare) {
    global->tuning.meanSxIsubstract = false;
    global->tuning.meanVsubstract = false;
    global->tuning.floatLoop = default_floatLoop ;
    Long ldAns = individuals;
    double * ans2 = (double*) MALLOC(sizeof(double) * ldAns * repetV);
    sum = sum1 = sum2 = 0.0;
    if (default_floatLoop < 0 ) {
      A_LD  = (LongDouble*) MALLOC(sizeof(*A_LD) * individuals * repetV);

      Long total = snps * repetV;
      if (sameInter)
	for (Long i=0; i<total; i++) inter2LD[i] = (LongDouble)inter2[i];

      vectorGeno_means_LongDouble(SxI2,
				inter2LD, repetV,
				global, utils, A_LD);
      for (Long i=0; i<individuals * repetV; i++) ans2[i] = (double) A_LD[i];
    } else { 
      vectorGeno_means_double(SxI2, inter2, repetV,
			      global, utils, ans2);
    }

    int limit = maxLines;
    int zaehler = 0;
    double diff = 0.0, maxdiff = 0.0;
    double diff2 = 0.0;
    for (Long i=0; i<individuals * repetV; i++) {
      double difference = ans[i]-ans2[i];
      maxdiff = MAX(maxdiff, FABS(difference));
	if (i < 5 ||
	  difference != 0 || (i == snps * repetV - 1 && limit > 0)) {
	//
	if (compare < 10) {
	  if (limit > 0) {
	    if (debug) PRINTF("i=%ld\t%7.3f\t%7.3f\td=%e\t%e\n",
		   i, ans[i], ans2[i], difference, difference/ans[i]);
	  }
	  else if (limit > -5 && debug) PRINTF(".");
	  limit--;
	}
       	if (false && FABS(difference) > 2e-16)
	  PRINTF("I=%ld\t%7.3f\t%7.3f\td=%e\t%e\n",
		 i, ans[i], ans2[i], difference, difference/ans[i]);
   	zaehler++;
	diff += FABS(difference);
 	diff2 += difference * difference;
      }
      sum += ans[i];
      sum1 += FABS(ans[i]);
      sum2 += ans2[i];
    }
    
    if (compare >= 10 && debug)
      PRINTF("NOTE! %d positions with average %e difference\n",
	     zaehler, diff / (double) zaehler);
    if (debug) PRINTF("\nFINAL sum=%5.5f sum2=%5.5f diff=%5.5e, %e same=%d maxdiff=%e",
		      sum, sum2, diff, diff2, sameInter, maxdiff);
    double thresh = sum1 * Tolerance1 * 2; //* SQRT(SQRT(individuals));
    if (sum != sum2 || diff > thresh) {
	PRINTF("\n***** vG differ by %e \n\n", sum - sum2);
      //  printf("%e > %e &&  %e > %e\n", FABS(sum - sum2), FABS(sum) * Tolerance , FABS(sum - sum2), Tolerance2);  
	if (FABS(sum - sum2) > FABS(sum) * Tolerance &&
	    FABS(sum - sum2) > Tolerance2)
	  PRINTF("***** %e > max(%e,%e)\n\n", sum - sum2,
		 FABS(sum) * Tolerance, Tolerance2);
	if (diff > thresh)
	  PRINTF("****diff: %e > %e\n", diff, thresh);
      BUG;
    } else PRINTF("diff %e < %e = thresh\n", diff, thresh);
 
    FREE(inter2);
    FREE(inter2LD);
    FREE(ans2);
    FREE(V_LD);
    FREE(A_LD);

    global->tuning.meanSxIsubstract = meanSxI;
    global->tuning.meanVsubstract = meanV;
    global->tuning.floatLoop = floatLoop;
  
  }
  if (debug > 1) PRINTF("END V * Relmatrix\n");

  FREE(inter);
}


void VectorCrossMatrix(SEXP SxI, double *V, Long repetV, 
		       option_type *global, utilsoption_type *utils,
		       double *ans) {
  // v * M * M^T   for  M = snp x indiv
  Long *info = GetInfo(SxI),
    individuals = info[INDIVIDUALS];
  double *inter = (double*) MALLOC(sizeof(double) * individuals* repetV);

  vectorGeno_means_double(SxI, V, repetV, global, utils, inter);
  genoVector_means_double(SxI, inter, repetV, global, utils, ans);
 
  FREE(inter);
}



