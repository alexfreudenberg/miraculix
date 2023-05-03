
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




#ifndef miraculix_vector_matrix_templateUint_H
#define miraculix_vector_matrix_templateUint_H 1

#include "Template.h"


#define coding_start(UINT,NAME)						\
  coding_header(UINT,NAME) {						\
  if (start_snp % CodesPerUnit != 0) BUG;				\
  Long  cur_snps = end_snp - start_snp,					\
    allUnits = Units(cur_snps),						\
    allUnitsM1 = allUnits - 1;						\
  unit_t *ans = Ans + start_snp / CodesPerUnit + start_individual * ldAns; \
  Long deltaIndiv = end_individual - start_individual;			\
  Ulong total = 0UL;							\
  int rest = (int) ((cur_snps - allUnitsM1 * CodesPerUnit) * BitsPerCode);    


#define coding_end(UINT,TRAFO)						\
  for (Long i=0; i<deltaIndiv; i++) {					\
    unit_t *code = ans + i * ldAns;					\
    UINT *mm = (UINT*) (M + i * ldM);			\
    bool w = false;							\
    for (Long j=0; j<allUnits; j++) {					\
      unit_t C = 0;							\
      int end = j == allUnitsM1 ? rest : BitsPerUnit;			\
      for (int shft = 0; shft < end; mm++, shft+=BitsPerCode) {	\
	w |= *mm > 2;							\
	if (w > 0) {							\
	  PRINTF("i=%ld<%ld %ld<%ld %d<%d : %u %ld\n",		\
		 i, deltaIndiv, j,allUnits, shft,end, (Uint) *mm, sizeof(UINT)); \
	  BUG;								\
	}								\
	C |= (unit_t) ((TRAFO) << shft);				\
      }									\
      code[j] = C;							\
    }									\
    total += w;								\
  }									\
  if (total > 0) ERR1("2bit#32: %ld values outside 0,1,2", total);	\
  }






#define get_matrix_start(UINT,NAME)					\
  get_matrix_header(UINT,NAME) {					\
  Long allUnitsM1 = Units(snps) - 1;					\
  int VARIABLE_IS_NOT_USED rest_k = (int) (snps - allUnitsM1 * CodesPerUnit);	    

#define get_matrix_end(UINT,TRAFO)				\
  for (Long i=0; i<individuals; i++) {				\
    UINT *a = (UINT*) (Ans + i * ldAns);	\
    unit_t *C = code + i * lda;				\
    for (Long j=0; j<=allUnitsM1; j++, C++) {		\
      unit_t twobit = *C TRAFO;						\
      int end_k = j==allUnitsM1 ? rest_k : CodesPerUnit;		\
      for (int k=0; k < end_k; k++) {					\
	*(a++) = twobit & CodeMask;					\
	twobit >>= BitsPerCode;				\
      }							\
    }							\
  }							\
  }


#define vectorGeno_start(TYPE, TRDTYPE, ENDTYPE, NAME)	\
  vectorGeno_header(TYPE,TRDTYPE,ENDTYPE,NAME)	 {	\
  int VARIABLE_IS_NOT_USED cores = GreaterZero(opt->cores);	\
  Long unitsM1 = Units(rows) - 1;				\
  int VARIABLE_IS_NOT_USED rest = (int) (rows - unitsM1 * CodesPerUnit)	\


#define vectorGeno_fctn(TYPE, TRDTYPE, ENDTYPE, TRAFO)			\
  for (Long r=0; r < repetV; r++) {					\
    for (Long i=0; i<cols; i++) {	\
      TYPE *vv = V + r * ldV;						\
      unit_t *pC = code + i * lda;					\
      ENDTYPE *a = Ans + r * ldAns /* cols */;			\
      LongDouble *aPF = pseudoFreq;				\
      TRDTYPE sum0 = 0;						\
      for (Long j=0; j<=unitsM1; j++, pC++) {				\
	unit_t twobit = *pC TRAFO;					\
	int end = j == unitsM1 ? rest : CodesPerUnit;			\
	if (aPF == NULL) {						\
	  for (int k=0; k<end; k++, twobit >>= BitsPerCode, vv++)	\
	    sum0 +=  (TRDTYPE) (*vv) * (TRDTYPE) (twobit & CodeMask);	\
	} else {							\
 	  for (int k=0; k<end; k++, twobit >>= BitsPerCode, vv++, aPF++) \
	    sum0 += (TRDTYPE) (*vv) *					\
	      ((TRDTYPE) (twobit & CodeMask) - (TRDTYPE) 2 * (TRDTYPE)(*aPF)); \
	}								\
      }									\
      a[i] = (ENDTYPE) sum0;						\
    }									\
  }									\
  }


#if defined sumGeno_start
#undef  sumGeno_start
#endif
#define sumGeno_start(NAME)			\
  sumGeno_header(NAME) {			\
  Ulong VARIABLE_IS_NOT_USED total = 0UL;		\
  Long VARIABLE_IS_NOT_USED units = Units(snps); 


#define sumGeno_fctn(TRAFO)						\
  for (Long i=0; i<individuals; i++) { /* ok */			\
    unit_t *pcode=code + lda * i;						\
    Ulong sum = 0UL;							\
    for (Long j=0; j<units; j++) {					\
      unit_t s = pcode[j];						\
      for (Uint u=0; u<CodesPerUnit; u++) {				\
	const unit_t value = s & CodeMask;				\
	sum += TRAFO(value);						\
	s >>= BitsPerCode;						\
      }									\
    }									\
    if (sums != NULL) sums[i] = sum;					\
    total += sum;							\
  }									\
  return total;								\
  }



#define gV_vG_start(TYPE, NAME)						\
  gV_vG_header(TYPE, NAME) {						\
  if (false) PRINTF("** vG_vG%s_%s\n", ""#NAME, ""#TYPE);		\
  basic_options *opt = &(utils->basic);					\
  tuning_options *tuning=&(global->tuning);				\
  SEXP next = getAttribPointer(SxI, Next);				\
  Long *info = GetInfo(SxI, true);					\
  coding_type coding = (coding_type) info[CODING];			\
  int variant = (int) info[VARIANT];					\
   bool gV = orig_gV;							\
   bool use_transposed =						\
    (gV xor uprightlyPacked(coding,variant)) && next!=R_NilValue;	\
   if (use_transposed) {						\
     if (false) PRINTF("USE TRANSPOSED\n");				\
    SxI = next;								\
    gV = !gV;								\
    info = GetInfo(SxI, true);						\
    coding = (coding_type) info[CODING];				\
  }									\
  bool is_transposed = transposed(coding);				\
  Long rows = info[is_transposed ? INDIVIDUALS : SNPS],		\
    cols = info[is_transposed ? SNPS : INDIVIDUALS];			\
 unit_t *code = Align(SxI, opt);					\
  if (false) PRINTF("%s %s ldAns=%ld >= %ld;  ldV=%ld >= %lu  trans=%d/%d meanSxi/V=%d/%d gV=%d/%d r/c=%ld/%ld\n", CODING_NAMES[coding],""#TYPE, ldAns, (gV ? rows : cols), ldV, (gV ? cols : rows), use_transposed, is_transposed, tuning->meanSxIsubstract, tuning->meanVsubstract, orig_gV, gV, rows, cols); \
  assert(ldAns >= (gV ? rows : cols));					\
  assert(ldV >= (gV ? cols : rows));					\
  if (false) PRINTF("end\n");						\
  LongDouble *pseudoFreq = NULL;					\
  if (SubstractMeanSxI)	{						\
    pseudoFreq = orig_gV ? getPseudoFreq(SxI) : getFreq(SxI);		\
  }									\
  vector_##TYPE##_t vD = NULL

  
#define gV_vG_end							\
  vD(code, rows, cols, info[LDA], coding, variant, pseudoFreq,		\
     V, repetV, ldV, opt, tuning, ans, ldAns)				\
  


#endif
