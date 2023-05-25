
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
#include "5codes.h"
#include "plink.h"
#include "Vector.matrix.h"
#include "OneByte.h"
#include "time.h"


void gV_vG_means_double(SEXP SxI, 
			double *V0, Long repetV, Long ldV, 
			bool gV, // else vG
			centering_type centered,
			normalizing_type normalized,
			bool meanV,
			bool meanSxI,
			option_type *global, utilsoption_type *utils,
			double *ans, Long ldAns) {
  STARTCLOCK;

  basic_options *opt = &(utils->basic);
  const int cores = GreaterZero(opt->cores);
  const Long *info = GetInfo(SxI, true),
    missings = info[MISSINGS];
  //  printf("%ld %ld %ld\n", missings, info[MISSINGS], NA_LONG);
  
  const Long
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    rows = gV ? snps : individuals,    
    cols = gV ? individuals : snps;
  bool PL = utils->basic.Cprintlevel > 5;
  ASSERT_LITTLE_ENDIAN;

  if (centered != RowMeans && centered != ColMeans && centered != NoCentering) {
    PRINTF("only centering by rowmeans and colmeans programmed\n");
    BUG;
  }
  
  const LongDouble
    c = (LongDouble) meanSxI,
    nu = (LongDouble) meanV,
    z = (LongDouble) (centered == RowMeans),
    zeta = (LongDouble) (centered == ColMeans),
    cMinusZ = (c - (gV ? zeta : z)) ;
  double *V = V0;
  
  LongDouble *freqV = (LongDouble*) CALLOC(repetV, sizeof(LongDouble));  
  LongDouble *mu = (LongDouble*) CALLOC(repetV, sizeof(LongDouble));

  LongDouble *freq = NULL, *sxiOne = NULL, *pseudofreq = NULL,
    *freqSxI = NULL, sqTotalSum = -1, pseudoSumFreq = -1, m = -1;
  
  //
  if (centered != NoCentering || meanV || meanSxI ||
      normalized != NoNormalizing) {
    CLOCK("before getFreq");
    getFreq(SxI, global, utils);
    freq =        gV ? getFreq(SxI)       : getPseudoFreq(SxI);
    sxiOne =      gV ? getSum(SxI)       : getPseudoSum(SxI);
    pseudofreq =  gV ? getPseudoFreq(SxI) : getFreq(SxI);
    freqSxI =     gV ? getFreqSxI(SxI)    : getPseudoFreqSxI(SxI);
    sqTotalSum =   gV ? getSqTotalSum(SxI) : getPseudoSqTotalSum(SxI);
    pseudoSumFreq =gV ? getPseudoSumFreq(SxI) : getSumFreq(SxI);
    m = getTotalSum(SxI);
    CLOCK("getFreq");
 }

  bool meanOrCZ  = meanV || cMinusZ != 0.0;
  if (meanOrCZ) { 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
     for (Long r=0; r < repetV; r++) {
      double *vv0 = V0 + r * ldV;
      LongDouble f = 0.0;
      for (Long j=0; j<cols; j++) {
	f += pseudofreq[j]  * (LongDouble) vv0[j];	
      }
      freqV[r] = f;
    }
  }
	
  if (meanV) {
    V = (double*) MALLOC(repetV * ldV * sizeof(*V));
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
    for (Long r=0; r < repetV; r++) {
      double *vv = V + r * ldV;
      double *vv0 = V0 + r * ldV;
      LongDouble mu0 = 0;
      for (Long i=0; i<cols; i++) mu0 += freqSxI[i] * (LongDouble) vv0[i];
      mu0 *= (LongDouble) 2.0 * (LongDouble) cols;
      mu0 += (LongDouble) 2.0 * m * c * (c - 2.0) * freqV[r];
      mu0 /= sqTotalSum + c * (c - 2.0) *  m * m / (LongDouble) rows;
      mu[r] = mu0;
      //   double mu1 = (double) mu0;
      for (Long i=0; i<cols; i++) vv[i] = (double) ((LongDouble) vv0[i] - mu0);
    }
  }
  if (PL && meanOrCZ) CLOCK(meanV ? "meanV" : "cMinusZ");
  
  gV_vG_double(SxI, V, repetV, ldV,
	       gV, meanSxI, global, utils, ans, ldAns);
  if (PL)
    CLOCK3("gV_vG (%s; %s; %s)", CODING_NAMES[info[CODING]],
	   gV ? "G * v" : "v * G",
	   z ? "row means" :  zeta ? "col means" : cMinusZ ? "cMinusZ"
	   : "no means");


  if (centered != NoCentering || meanV || meanSxI) {    
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
    for (Long r=0; r < repetV; r++) {
      double *vv0 = V0 + r * ldV;
      double *a = ans + r * ldAns;
      LongDouble sxiOneFactor = nu * mu[r];
      LongDouble freqFactor = 0.0;
      if ( (gV ? z : zeta) != 0.0) {
	for (Long i=0; i<cols; i++) freqFactor += ((LongDouble) vv0[i]);
	freqFactor *= -2.0L;
      }      
      LongDouble onesFactor = cMinusZ * freqV[r];
      onesFactor -= c * nu * mu[r] * pseudoSumFreq;
      onesFactor *= 2.0L;

      //      printf("c-z=%4.3e frqV=%4.3e c=%4.3e nu=%4.3e mu=%4.3e pF=%4.3e factor=%4.3e,%4.3e \n",(double) cMinusZ, (double) freqV[r], (double) c , (double) nu , (double)  mu[r], (double) pseudoSumFreq, (double) freqFactor, (double) onesFactor );
      
      for (Long j=0; j<rows; j++) {
	// hier geht gegenueber LongDouble die Genauigkeit verloren:
	// wenn in *LD.cc  (double) (LongDouble) ((LongDouble) (double) a[j]
	// gesetzt wird ist selbst bei Centered kein Unterschied mehr da!

	a[j] = (double) ((LongDouble) a[j] + sxiOneFactor * sxiOne[j] +
			 freqFactor * freq[j] + onesFactor);
      }
      // printf("a[0]=%4.3f\n", (double) a[0]);
     }
    if (PL)  CLOCK("centering");
  }


  // Missings
  if (centered != NoCentering) {
     SEXP Miss = getAttribPointer(SxI, Missings);
     if (Miss != R_NilValue) {
       Long *mssngs = LONG(Miss);
       assert(missings > 0);
       Long end_i = 2 * missings ;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
      for (Long r=0; r < repetV; r++) { 
	 double *vv0 = V0 + r * ldV;
	 double *a = ans + r * ldAns;
	 LongDouble Z = gV ? z : zeta;
	 LongDouble Zeta = gV ? zeta : z;

	 for (Long i=0; i<end_i; i+=2) {
	   Long alpha = mssngs[i + 1 - gV];
	   Long beta = mssngs[i + gV];
	   //a[alpha] += 2.0 * freq[alpha] * vv0[ mssngs[i+1] ];
	   a[alpha] = (double) ((LongDouble) a[alpha] +
				(Z * freq[alpha] + Zeta * pseudofreq[beta])
				* 2.0L * (LongDouble) vv0[beta]);
	 }
	 //   printf("(double) a[0]=%4.3f\n", (double) a[0]);
       }
      if (PL) CLOCK("missings");
    } else {
       assert(missings == 0);
     }
  }
  
  
   if (normalized == GlobalNormalizing) {
     LongDouble sigmaSq = gV ? getSigmaSq(SxI) : getPseudoSigmaSq(SxI);
     double sigma = (double) sqrtl(sigmaSq);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
     for (Long r=0; r < repetV; r++) {
       double *a = ans + r * ldAns;
       for (Long i=0; i < rows; i++) a[i] /= sigma;
     }
     if (PL) CLOCK("finalizing") else CLOCK("gV_vG");
  } else if (normalized != NoNormalizing) BUG; // not programmed yet

   //   printf("leaving gv_meansx\n");

  if (meanV) {FREE(V);}
  FREE(mu);
  FREE(freqV);
}



void genoVector_means_double(SEXP SxI, 
			     double *V, Long repetV, Long ldV, 
			     option_type *global, utilsoption_type *utils,
			     double *ans, Long ldAns) {
         gV_vG_means_double(SxI, V, repetV, ldV, 
		     true,
		     global->genetics.centered,
		     global->genetics.normalized,
		     global->tuning.meanVsubstract,
		     global->tuning.meanSxIsubstract,
		     global, utils, ans, ldAns);
}



void genoVector_means_double(SEXP SxI, 
			     double *V, Long repetV,
			     option_type *global, utilsoption_type *utils,
			     double *ans) {
  Long *info = GetInfo(SxI);
  genoVector_means_double(SxI, V, repetV, info[INDIVIDUALS], 
			  global, utils, ans, info[SNPS]);
}

void vectorGeno_means_double(SEXP SxI, 
			     double *V, Long repetV, Long ldV, 
			     option_type *global, utilsoption_type *utils,
			     double *ans, Long ldAns) {
  gV_vG_means_double(SxI, V, repetV, ldV,
		     false, 
		     global->genetics.centered,
		     global->genetics.normalized,
		     global->tuning.meanVsubstract,
		     global->tuning.meanSxIsubstract,
		     global, utils, ans, ldAns);  
}


void vectorGeno_means_double(SEXP SxI, 
			     double *V, Long repetV, 
			     option_type *global, utilsoption_type *utils,
			     double *ans) {
  Long *info = GetInfo(SxI);
  vectorGeno_means_double(SxI, V, repetV, info[SNPS],
			  global, utils, ans, info[INDIVIDUALS]);  
}

