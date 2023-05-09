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



// HERE:  TwoBitHaplo-Codes


#define BitsPerCode 2
#define MY_VARIANT 32

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "Template.h"
#include "utils_miraculix.h"
#include "MX.h"
#include "Haplo.h"
#include "bitUint.h"


void assert_haplo(SEXP Code) {
  Long *info = GetInfo(Code);
  if (isHaplo(info[CODING])) ERR0("not a haplotype matrix");
}


unit_t GetTwoBitHaplo(unit_t *code, Long snp) {
  static unit_t rev_haplo_code[4] = {0U, 1U, 1U, 2U};
  return rev_haplo_code[ExtractCode(code, snp)];
}

// AUF KEINEN FALL "PARALLEL" WEGEN RANDOM
#define Inner(mm12) (double *freq1, double *freq2,		\
		     Long snps,Long individuals, coding_type coding, \
		     unit_t *ans, Long ldAns) {				\
    Long bitsPerCodeCompr, shiftIncrCompr, deltaCompressed, CpU;	\
    getHaploIncr(snps, individuals, ldAns, UnknownSNPcoding,		\
		 0, coding,				\
		 0, NULL, NULL, NULL, NULL,				\
		 &bitsPerCodeCompr, &shiftIncrCompr, &deltaCompressed,	\
		 &CpU);							\
    int bPCC = (int) bitsPerCodeCompr;					\
    Long allUnitsM1 = snps2entities(snps, CpU) - 1;			\
    int rest = (int) ((snps - allUnitsM1 * CpU) * (1 + shiftIncrCompr)); \
    for (Long i=0; i<individuals; i++, ans += ldAns) {	/* ok */	\
      Long mm = 0;							\
      for (Long j=0; j<=allUnitsM1; j++) {				\
	unit_t dummy1 = 0, dummy2 = 0;					\
	int end = j == allUnitsM1 ? rest : BitsPerUnit;		\
	for (int shft=0; shft<end && mm < snps; mm++) {		\
	  double f2 = freq2[mm];					\
	  mm12;								\
	  dummy2 |= mm2 << (shft + shiftIncrCompr);			\
	  dummy1 |= mm1 << shft;					\
	  shft += bPCC;					\
	}								\
	if (shiftIncrCompr) ans[j] = dummy1 | dummy2;			\
	else {								\
	  ans[j] = dummy1;						\
	  ans[j + deltaCompressed] = dummy2;				\
	}								\
      }									\
    }									\
  }

void InnerDetermTwoBitHaplo Inner(unit_t mm1 = freq1[mm] == 1.0;	 
			    unit_t mm2 = ISNA(f2) ? mm1 : (f2  == 1.0))

void InnerRandomTwoBitHaplo Inner(unit_t mm1 = UNIFORM_RANDOM <= freq1[mm];
			    unit_t mm2 = ISNA(f2) ? mm1 : UNIFORM_RANDOM <= f2)
  
  
