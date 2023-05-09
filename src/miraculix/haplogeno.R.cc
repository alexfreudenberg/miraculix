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

#if defined compatibility_to_R_h

#include "miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "utils_miraculix.h"
#include "extern.h"
#include "mmagpu.h"
#include "transform.h"

#include "haplogeno.h"
#include "Template.h"
#include "intrinsicsCheck.h"
#include "intrinsics_specific.h"

#include "MXinfo.h"
#include "Haplo.h"
#include "Files.h"
#include "1bit.h"
#include "2bit.h"
#include "3bit.h"
//#include "OneByte.h"



SEXP get_centered() {
  option_type *global = &(KEYT_M()->global);
  int len = global->genetics.ncentered ;
  double *centered = global->genetics.pcentered;
  assert (centered != NULL && global->genetics.centered==User);
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, len));
  MEMCOPY(REAL(Ans), centered, sizeof(double) * len);
  UNPROTECT(1);
  return Ans;
}


SEXP createSNPmatrix(SEXP SNPs, SEXP Individuals) {
  int 
    n = Int0(Individuals),
    snps = Int0(SNPs);
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
  coding_type coding = KT->global.genetics.coding;
  int
    variant = check_variant(coding, KT->global.tuning.variant, opt->efficient);
  return createSNPmatrix(snps, n, coding, variant, 0, global, utils);
}


Long crossprodI(unit_t * Code, Long snps, Long individuals, Long lda,
		 coding_type coding, int variant,
		 centering_type centered, normalizing_type normalized,
		 bool squared,
		 LongDouble *Prec,
		 basic_options *opt,
		 tuning_options *tuning,
		 double *A);


Long t_crossprodI(unit_t * Code, Long snps, Long individuals, Long lda,
		 coding_type coding, int variant,
		 centering_type centered, normalizing_type normalized,
		 bool squared,
		  LongDouble *Prec,
		   basic_options *opt,
		 tuning_options *tuning,
		 double *A);


SEXP crossprod(SEXP SxI) {
  KEY_type *KT = KEYT_M();
  option_type * global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
  Long *info = GetInfo(SxI);
  Long
     individuals = info[INDIVIDUALS],
     lda = info[LDA];
   if (info[CODING] == NA_LONG) ERR0("looks like an uninitialised matrix");
   coding_type coding = (coding_type) ;
   if (isHaplo(coding))
     ERR0("matrix is a haplotype matrix, not a genomic matric");
   if (coding == UnknownSNPcoding) ERR0("not a coded Z matrix");
   
  int variant = check_variant(coding, (int) info[VARIANT], opt->efficient),
    snps = info[SNPS];
  unit_t *code = Align(SxI, opt);
  SEXP Ans = PROTECT(allocMatrix(REALSXP, individuals, individuals));
  double *ans = REAL(Ans);

  getFreq(SxI, global, utils);
  BUG;
  crossprodI(code, snps, individuals, lda, coding, variant, 
	     KT->global.genetics.centered,
	     KT->global.genetics.normalized,
	     KT->global.genetics.squared,
	     LONGREAL(getAttribPointer(SxI, Precise)),
	     opt,
	     &(KT->global.tuning),
	     ans);
 
  UNPROTECT(1);
  return Ans;
}

SEXP LD(SEXP SxI) {
  KEY_type *KT = KEYT_M();
  option_type * global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
  Long  *info = GetInfo(SxI);
  if ((int) coding == NA_LONG) ERR0("looks like an uninitialised matrix");
  coding_type coding = (coding_type) info[CODING];
  if (isHaplo(coding))
     ERR0("matrix is a haplotype matrix, not a genomic matric");
  if (coding == UnknownSNPcoding) ERR0("not a coded Z matrix");
  if (uprightlyPacked(coding))
    ERR1("'%s' does not allow LD calculation yet", CODING_NAMES[coding]);

  int variant = check_variant(coding, (int) info[VARIANT], opt->efficient),
    snps = info[SNPS],
    lda = info[LDA],
    individuals = info[INDIVIDUALS];

  SEXP Ans = PROTECT(allocMatrix(REALSXP, snps, snps));
  double *ans = REAL(Ans);
  
  SEXP next = getAttribPointer(SxI, Next);
  if (next != R_NilValue) {
    crossprodI(Align(next, opt),
	       individuals, snps, // reverse ordering of snps, individuals
	       lda, switch_transposed(coding), variant, 
	       RowMeans, AccordingCentering,
	       false,
	       LONGREAL(getAttribPointer(SxI, Precise)),
	       opt,
	       &(KT->global.tuning),
	       ans);
  }

  getFreq(SxI, global, utils);
  t_crossprodI(Align(SxI, opt), snps, individuals, lda, coding, variant, 
	       ColMeans, AccordingCentering,
	       false,
	       LONGREAL(getAttribPointer(SxI, Precise)),
	       opt,
	       &(KT->global.tuning),
	       ans);
 
  UNPROTECT(1);
  return Ans;
}




SEXP zeroGeno(SEXP SxI, SEXP Snps, SEXP Individuals, SEXP Copy) {
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
  Long *info = GetInfo(SxI);
  Uint
    *snps = (Uint*) INTEGER(Snps),
    lenSnps = LENGTH(Snps),
    *Indiv =(Uint*) INTEGER(Individuals),
    lenIndiv = LENGTH(Individuals);
  unit_t
    *code = Align(SxI, opt);
  coding_type coding = (coding_type) info[CODING];
  ASSERT_LITTLE_ENDIAN;
 
  SEXP Ans = SxI;

  if (LOGICAL(Copy)[0]) {
    Ans = PROTECT(createSNPmatrix(info[SNPS], info[INDIVIDUALS], coding,
				  KT->global.tuning.variant, 0,
				  false, R_NilValue, global, utils));
    info = GetInfo(Ans);
    MEMCOPY(Align(Ans, opt), code, info[MEMinUNITS] * BytesPerUnit); // OK
  }
    
  switch(coding) {
  case OneBitGeno :
  case OneBitHaplo :
    zeroGeno1(Ans, snps, lenSnps, Indiv, lenIndiv, opt);break;
  case TwoBitGeno :case TwoBitHaplo :
    zeroGeno2(Ans, snps, lenSnps, Indiv, lenIndiv, opt);break;
  case ThreeBit :
    zeroGeno3(Ans, snps, lenSnps, Indiv, lenIndiv, opt);break;
  default :  ERR0("zero-ing does not work (yet) for given coding");
  }

  if (LOGICAL(Copy)[0]) UNPROTECT(1);

  return Ans;
}



SEXP allele_freq(SEXP SxI) {
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  Long *info = GetInfo(SxI),
    snps = info[SNPS];
  SEXP Ans = PROTECT(allocVector(REALSXP, snps));
  allele_freq_double(SxI, LONG(getAttribPointer(SxI, Missings)),
		     info[MISSINGS],
		     false, global, utils, NULL, REAL(Ans));
  UNPROTECT(1);
  return Ans;
}


SEXP matrix_get(SEXP SxI) {
  KEY_type *KT = KEYT_M();
  basic_options *opt = &(KT->global_utils.basic);
  Long *info = GetInfo(SxI),
    individuals = info[INDIVIDUALS],
    snps =  info[SNPS];
  ASSERT_LITTLE_ENDIAN;
  SEXP Ans =  PROTECT(allocMatrix(INTSXP, snps, individuals));
  matrix_get_4Byte(Align(SxI, opt), snps, individuals, info[LDA],
		   (coding_type) info[CODING], GreaterZero(opt->cores),
		   (unit_t*) INTEGER(Ans), SnpNatLDA(snps, FourByteGeno));
  UNPROTECT(1);
  return Ans;
}


SEXP fillSNPmatrix(SEXP Z, SEXP Idx, SEXP Vector) {
  KEY_type *KT = KEYT_M();
  basic_options *opt = &(KT->global_utils.basic);
  ToIntLocal(Idx);
  Long *infoZ = GetInfo(Z),
    *infoVector = GetInfo(Vector),
    snpsZ = infoZ[SNPS],
    snpsVector = infoVector[SNPS],
    individuals = infoZ[INDIVIDUALS];
  if (snpsZ != snpsVector) ERR0("number of snps differ.");
  coding_type     
    codingZ =  (coding_type) infoZ[CODING],
    codingVector = (coding_type) infoVector[CODING];
  if (!isRoughlyGeno(codingZ)) ERR0("not a geno matrix that is filled");
  if (codingZ != codingVector) ERR0("codings of 'SxI' and 'values' differ.");
  if (infoZ[VARIANT] != infoVector[VARIANT]) ERR0("variants differ");
  if (infoZ[BIGENDIAN] != infoVector[BIGENDIAN]) ERR0("endians differ");
  
  unit_t			     
    *z = Align(Z, opt),
    *v = Align(Vector, opt);
  Long  len = LENGTH(Idx);
   
  Long
    ldZ = GetLDA(snpsZ, individuals, codingZ, infoZ[LDABITALIGN]),
    ldV = GetLDA(snpsVector, individuals, codingVector,
		  infoVector[LDABITALIGN]),   
    bytes = BytesPerUnit * MIN(ldZ, ldV);// OK
  
  for (Long i=0; i<len; i++, v += ldV) {
    Long idx = Idxint[i] - 1;
    if (idx >= individuals) ERR0("Idx out of bound");
    MEMCOPY(z + idx * ldZ, v, bytes);
  }
  FREElocal(Idx);
  return R_NilValue;
}



SEXP substract_centered(SEXP SnpXindiv) {
  option_type *global = &(KEYT_M()->global);
  Long individuals = ncols(SnpXindiv),
    len = global->genetics.ncentered;
  Long snps = nrows(SnpXindiv);
  SEXP Ans;
  PROTECT(Ans = allocMatrix(REALSXP, snps, individuals));
  double *ans = REAL(Ans),
    *snpXindiv = REAL(SnpXindiv),
    *centered = global->genetics.pcentered;
  assert(centered != NULL && global->genetics.centered==User);
  if (snps != len) ERR0("length of 'centered' must equal the number of SNPs.");

#ifdef DO_PARALLEL // omp -- never delete this comment //
 KEY_type *KT = KEYT_M();
 int cores = GreaterZero(KT->global_utils.basic.cores);
#endif
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) {
    Long i_snps = i * snps;
    double *s = i_snps + snpXindiv,
      *a = ans + i_snps;
    for (Long j=0; j<snps; j++) a[j] = s[j] - centered[j];
  }
  UNPROTECT(1);
  return Ans;
}


bool debugging = false;
SEXP Debug() { debugging = true; return R_NilValue; }
SEXP StopDebug() { debugging = false;  return R_NilValue; }

SEXP is2BitMethod(SEXP Meth) { // 'Method' is global already
  SEXP Ans;
  int coding = INTEGER(Meth)[0];
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = is2Bit(coding);
  UNPROTECT(1);
  return Ans;
}


#endif
