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


#define BitsPerCode 2
#define MY_VARIANT 32


#include "Basic_miraculix.h"

#if defined compatibility_to_R_h
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "utils_miraculix.h"
#include "MX.h"
#include "Template.h"
#include "Haplo.h"
#include "miraculix.h"

#include "options_RFU.h"
#include "xport_import.h"

SEXP rhaplomatrixPart2(SEXP Freq1, SEXP Freq2, SEXP Code) {
  //  printf("length =%d %d\n", LENGTH(Freq1), LENGTH(Freq2));
  PROTECT(Code);
  KEY_type *KT = KEYT_M();
  basic_options *opt = &(KT->global_utils.basic);
  double *freq1 = REAL(Freq1),
    *freq2 = REAL(Freq2);
  Long
    *info = GetInfo(Code),
    snps = LENGTH(Freq1);
  unit_t
    *C =  Align(Code, opt);
  assert(snps == info[SNPS] && LENGTH(Freq2) == LENGTH(Freq1));
 
  bool random = false;  
  for (Long s=0; s < snps; s++) {
    double f = freq1[s];
    if (f < 0 || f > 1) ERR0("frequencies not in [0,1]");
    random |= f != 0 && f != 1;
    f = freq2[s];
    if (!ISNA(f)) {
      if (f < 0 || f > 1) ERR0("frequencies not in [0,1]");
      random |= f != 0 &&  f != 1;
    }
  }

  if (random) {
    GetRNGstate();
    InnerRandomTwoBitHaplo(freq1, freq2, info[SNPS], info[INDIVIDUALS],
			   (coding_type) info[CODING], C, info[LDA]);
    PutRNGstate();
  } else {
    InnerDetermTwoBitHaplo(freq1, freq2, info[SNPS], info[INDIVIDUALS],
			   (coding_type) info[CODING], C, info[LDA]);
  }
  UNPROTECT(1);
  return Code;
}

SEXP rhaplomatrix(SEXP Freq1, SEXP Freq2, SEXP Individuals, SEXP Coding) {
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  Long individuals = INTEGER(Individuals)[0],
    snps = LENGTH(Freq1);
  coding_type coding = (coding_type) INTEGER(Coding)[0];
  assert(coding == OneBitHaplo || coding == TwoBitHaplo);
  if (snps != (Long) LENGTH(Freq2))
    ERR0("'freq' and 'freq2' do not have the same length");
  SEXP Ans = createSNPmatrix(snps, individuals, coding, // no PROTECT( needed;
			     MY_VARIANT, MY_LDABITALIGN, global, utils);
  return rhaplomatrixPart2(Freq1, Freq2, Ans);
}

#endif
