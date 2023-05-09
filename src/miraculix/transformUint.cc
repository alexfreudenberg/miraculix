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
#include "options.h"
#include "haplogeno.h"
#include "Template.h"
#include "utils_miraculix.h"
#include "MXinfo.h"
#include "intrinsicsCheck.h"
#include "Files.h"
#include "Haplo.h"
#include "1bit.h"
#include "2bit.h"
#include "3bit.h"
#include "5codes.h"
#include "4Byte.h"
#include "OneByte.h"
#include "plink.h"
#include "intrinsics_specific.h"


coding_type OneSet2G(coding_type coding) {
  coding_type tmp;
  if      (coding == OneByteHaplo) tmp = OneByteGeno;
  else if (coding == EightByteHaplo ||
	   coding == FourByteHaplo ||
	   coding == FourByteSingleBit) tmp = FourByteGeno;
  else if (coding == FourByteSingleBit_Tr ||	     
	   coding == EightByteHaplo_Tr) tmp = FourByte_Tr;
  else BUG;
  return tmp;
} 

coding_type swapHaploGeno(coding_type coding, bool ToH, bool ToG,
			  usr_bool stopOnError) {
  switch(coding) {
  case OneBitGeno:   if (ToH) coding = OneBitHaplo;
    break;
  case TwoBitGeno:   if (ToH) coding = TwoBitHaplo;
    break;
  case OneByteGeno:  if (ToH) coding = OneByteHaplo;
    break;
  case FourByteGeno: if (ToH) coding = FourByteHaplo;
    break;
  case OneBitHaplo:  if (ToG) coding = OneBitGeno;
    break;
  case TwoBitHaplo:  if (ToG) coding = TwoBitGeno;
    break;
  case OneByteHaplo: if (ToG) coding = OneByteGeno;
    break;
  case FourByteHaplo:if (ToG) coding = FourByteGeno;
    break;
  default :
    if (stopOnError == True || 
	(stopOnError == Nan && ( coding == AutoCoding ||
				 coding == UnknownSNPcoding ||
				 (ToH && !isHaplo(coding)) ||
				 (ToG && isHaplo(coding)) )))
      ERR2("no counterpart within %.20s-coding found for '%.20s'.",
	   ToH ? "haplo" : "geno", CODING_NAMES[coding]);
  }
  return coding;
}

coding_type swapHaploGeno(coding_type coding, bool ToH) {
  return swapHaploGeno(coding, ToH, !ToH, Nan);
}


typedef void (*get_1Byte_t )(unit_t *, Long snps, Long indiv, Long lda,
			     int cores, unit_t *ans, Long ldAns);

void matrix_get_1Byte(unit_t *Code, Long snps, Long indiv, Long lda,
		      coding_type coding, int cores,
		      unit_t *ans, Long ldAns) {
  // BUG;
  MEMSET(ans, 0, totalMem(snps, indiv, lda, coding, false) * BytesPerUnit);
  get_1Byte_t get = NULL;
  switch (coding) {
  case OneBitGeno : get = get_matrix1_1Byte; break;
  case TwoBitGeno : get = get_matrix2_1Byte ; break;
  case ThreeBit :   get = get_matrix3_1Byte; break;
  case Plink :   get = get_matrixPlink_1Byte; break;
  case FourByteGeno:get = get_matrixPlain_1Byte; break;
  case OneByteGeno: get = get_matrixOneByte_1Byte; break;
  default : BUG;
  }
  get(Code, snps, indiv, lda, cores, ans, ldAns);
}


coding_t  matrix_coding_4Byte(coding_type coding, int variant,
				  bool bigendian) {
  coding_t matrix_coding;
  switch (coding) {
  case OneBitGeno : matrix_coding = coding1_4Byte; break;
  case TwoBitGeno :
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 : // matrix_coding = coding_2v256; break; -- Fehler
    case 128 :
      if (bigendian) {
	matrix_coding = coding_2v128_4Byte;
	break;
      }
      FALLTHROUGH_OK;
    default:
       matrix_coding = coding2_4Byte;
    }
    break;
  case TwoBitGenoTransposed : matrix_coding = coding2trans_4Byte; break;
  case PlinkTransposed : matrix_coding = codingPlinktrans_4Byte; break;
  case ThreeBit : matrix_coding = coding3_4Byte; break;
  case OneByteGeno : matrix_coding = codingOneByte_4Byte; break;
  case FourByteGeno: matrix_coding = codingPlain_4Byte; break;
  case Plink :
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 : // matrix_coding = coding_2v256; break; -- Fehler
    case 128 :
      FALLTHROUGH_OK;
    default:
      matrix_coding = codingPlink_4Byte;
    }
    break;
  case FiveCodes : matrix_coding = coding5_4Byte; break;
  case FiveCodesTransposed : matrix_coding = coding5trans_4Byte; break;
  default : BUG;
  }
  return matrix_coding;
}


coding_t matrix_coding_1Byte(coding_type coding, int variant,
			     bool VARIABLE_IS_NOT_USED bigendian) {
  coding_t matrix_coding;
  switch (coding) {
  case OneBitGeno : matrix_coding = coding1_1Byte; break;
  case TwoBitGeno :
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 :
      // if (bigendian) {matrix_coding = coding_2v256_1Byte; break; }
      FALLTHROUGH_OK;
    case 128 :
      // if (bigendian) {matrix_coding = coding_2v128_1Byte; break;}
      FALLTHROUGH_OK;
    default: matrix_coding = coding2_1Byte;
    }
    break;
  case TwoBitGenoTransposed : matrix_coding = coding2trans_1Byte; break;
  case PlinkTransposed : matrix_coding = codingPlinktrans_1Byte; break;
  case ThreeBit : matrix_coding = coding3_1Byte; break;
  case FourByteGeno: matrix_coding = codingPlain_1Byte; break;
  case OneByteGeno : matrix_coding = codingOneByte_1Byte; break;
    //  case FiveCodes : matrix_coding = coding5_1Byte; break;
  case Plink :
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 :
      // if (bigendian) {matrix_coding = coding_2v256_1Byte; break; }
      FALLTHROUGH_OK;
    case 128 :
      // if (bigendian) {matrix_coding = coding_2v128_1Byte; break;}
      FALLTHROUGH_OK;
    default: matrix_coding = codingPlink_1Byte;
    }
    break;
  case FiveCodes : matrix_coding = coding5_1Byte; break;
  case FiveCodesTransposed : matrix_coding = coding5trans_1Byte; break;
 default : BUG;
  }
  return matrix_coding;
}




#define SelectGeno(UINT,UCHAR)						\
  void SelectGeno##UINT##UCHAR(unit_t *Old, Long OldSnps, Long OldIndiv, \
			       Long OldLDA,				\
			       int *snp_list, Long snp_list_N,		\
			       int *indiv_list, Long indiv_list_N,	\
			       unit_t *Ans, Long ldAns) {		\
    /*//    printf("select %u %u (%lu %lu)\n", OldLDA,ldAns, sizeof(UINT), sizeof(UCHAR)	); */ \
    Long indiv_N = indiv_list_N > 0 ? indiv_list_N : OldIndiv;		\
    Long nIndiv = 0;							\
    for (Long i=0; i<indiv_N; i++) {					\
      Long ii=i;							\
      if (indiv_list_N > 0) {						\
	ii = indiv_list[i];						\
	if (ii >= OldIndiv) continue;					\
      }									\
      UCHAR *neu = (UCHAR*) (Ans + (nIndiv++) * ldAns);	\
      MEMSET(neu, 0, BytesPerUnit * ldAns); /* avoids undefined memory */ \
      if (ii < 0)  continue;						\
      UINT *old = (UINT*) (Old + OldLDA * ii);		\
      if (snp_list_N == 0) {						\
	if (sizeof(UCHAR) == sizeof(UINT)) /* lda may differ */		\
	  MEMCOPY(neu, old, sizeof(UINT) * OldSnps);			\
	else for (Long j=0; j<OldSnps; j++) neu[j]=(UCHAR) old[j];	\
	continue;							\
      }									\
      for (Long j=0; j<snp_list_N; j++) {				\
	if (snp_list[j] >= OldSnps) continue;				\
	/*//printf("i=%d %d: %d %d; %lu %u n=%d\n", i, j, ii, snp_list[j], \
	  (Long) neu - (Long) Ans, ldAns, nIndiv); */			\
	if (snp_list[j] >= 0) neu[j] = (UCHAR) old[snp_list[j]];	\
      }									\
    }									\
  }

SelectGeno(_4Byte,_4Byte)
SelectGeno(_1Byte,_1Byte)
SelectGeno(_4Byte,_1Byte)
SelectGeno(_1Byte,_4Byte)



#define SelectHaplo(UINT,UCHAR)						\
  void SelectHaplo##UINT##UCHAR(unit_t *Old, Long OldSnps, Long OldIndiv, \
				    Long OldLDA, coding_type OldCoding,	\
				    int *snp_list, Long snp_list_N,	\
				    int *indiv_list, Long indiv_list_N,	\
				    usr_bool set1,			\
				    unit_t *Ans, Long AnsSnps, \
				    Long ldAns, coding_type AnsCoding) { \
    /* currently, always indiv_list=NULL as indiv already selected */	\
    const bool both = set1 == Nan,					\
      doubled = both & UnCompressedDoubleCols(OldCoding, true);		\
    assert(UnCompressedDoubleCols(OldCoding, true) ==			\
	   UnCompressedDoubleCols(AnsCoding, true));			\
    /*//    printf("yX %s %s \n",CODING_NAMES[OldCoding], CODING_NAMES[AnsCoding]); */ \
    assert(isSomeByteHaplo(OldCoding, !both) && isSomeByteHaplo(AnsCoding, !both)); \
    assert(!transposedHaplo(OldCoding) && !transposedHaplo(AnsCoding));	\
    Long oldIncr, oldSndIncr, ansIncr, ansSndIncr, oldNext, ansNext,	\
      indiv_N = indiv_list_N > 0 ? indiv_list_N : OldIndiv,		\
      snps_N = snp_list_N > 0 ? AnsSnps : OldSnps;			\
    if ((sizeof(UINT)==1) xor (OldCoding == OneByteHaplo)) BUG;		\
    if ((sizeof(UCHAR)==1) xor						\
	(AnsCoding == OneByteHaplo || AnsCoding == OneByteGeno)) BUG;	\
    getHaploIncr(OldIndiv, OldLDA, OldCoding, 2,  /* in UINT ! */	\
		 &oldNext, /* in UINT */				\
		 &oldSndIncr, /* in unit_t */				\
		 &oldIncr /* in unit_t */					\
		 );							\
    getHaploIncr(indiv_N, ldAns, AnsCoding, 1 + both,			\
		 &ansNext, &ansSndIncr, &ansIncr);			\
    if (set1 == False) {						\
      Old += oldSndIncr;						\
      ansSndIncr = 0;							\
    }									\
     									\
    Long min = MIN(OldLDA, ldAns) * BytesPerUnit,			\
      nIndiv = 0;							\
    for (Long i=0; i<indiv_N; i++) {					\
      Long ii = i;							\
      if (indiv_list_N > 0) {						\
	ii = indiv_list[i];						\
	if (ii >= OldIndiv) continue;					\
      }									\
      UCHAR *neu = (UCHAR*) (Ans + (nIndiv++) * ansIncr);	\
      UCHAR *neu2 = (UCHAR*) (((unit_t*) neu) +  ansSndIncr)	;	\
      MEMSET(neu, 0, ldAns * BytesPerUnit);				\
      if (doubled) MEMSET(neu2, 0, ldAns * BytesPerUnit);		\
      if (ii < 0) continue;						\
      UINT *old = (UINT*) (Old + oldIncr * ii);		\
      UINT *old2 = (UINT*) (((unit_t*) old) + oldSndIncr); /* // ((( OK */	\
      /*//   printf("%d %d %d %d %d ==%d %d=%d AnsSize=%d %lu D=%d\n", i, ii, snp_list_N, snps_N, indiv_N, sizeof(UCHAR) \
	== sizeof(UINT), ansIncr, ldAns, sizeof(UCHAR), Ans, (Long) neu2 - (Long) neu); */ \
      if (snp_list_N == 0) {						\
	if (sizeof(UCHAR) == sizeof(UINT)) {				\
	  MEMCOPY(neu, old, min);					\
	  if (doubled) MEMCOPY(neu2, old2, min);			\
	} else for (Long j=0; j<snps_N; j++) {				\
	    /*// printf("j=%d %d %d \n", j, old[j], old2[j]);	*/	\
	    neu[j] = (UCHAR) old[j];					\
	    if (both) neu2[j] = (UCHAR)	old2[j];			\
	  }								\
	continue;							\
      }									\
      for (Long j=0; j<snp_list_N; j++) {				\
	if (snp_list[j] >= OldSnps) continue;				\
	if (snp_list[j] >= 0) {						\
	  Long idx = oldNext * snp_list[j];				\
	  *neu =  (UCHAR) old[idx];					\
	  if (both) *neu2 = (UCHAR) old2[idx];				\
	} 								\
	neu += ansNext;							\
 	neu2 += ansNext;						\
      }									\
    }									\
  }

SelectHaplo(_1Byte,_1Byte)
SelectHaplo(_4Byte,_1Byte)
SelectHaplo(_1Byte,_4Byte)
SelectHaplo(_4Byte,_4Byte)


coding_type haplo2geno(unit_t *X, Long snps, Long individuals, Long lda,
		       coding_type coding,
		       int main_var, int cores,
		       unit_t *ans, Long ldAns, coding_type anscoding) {
  // endian need not be considered
  switch (coding) {
  case OneBitHaplo :
    switch (main_var) {
    case 512 :
    case 256 :
      OneBithaplo2geno256(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    case 128 :
      OneBithaplo2geno128(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    case 64:
      OneBithaplo2geno64(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    case 32:
      OneBithaplo2geno32(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    default: BUG;
    }
    return OneBitGeno;
	
  case TwoBitHaplo :
    if (anscoding == ThreeBit) {
      TwoBithaplo2geno3(X, snps, individuals, lda, cores, ans, ldAns);
      return ThreeBit;
    }  
    switch (main_var) {
    case 512 :
    case 256 :      
      TwoBithaplo2geno256(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    case 128 :
      TwoBithaplo2geno128(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    case 64:
    case 32:
      TwoBithaplo2geno2(X, snps, individuals, lda, cores, ans, ldAns);
      break;
    default: BUG;
    }
    return TwoBitGeno;
  default: BUG;
  }
}


void haplo2geno(unit_t *code, Long snps, Long individuals, Long lda,
		coding_type oldCoding, int variant,int cores) {//MoBPS
  haplo2geno(code, RoundUpSnps(snps, oldCoding, variant),
	     individuals, lda, oldCoding, main_variant(variant),
	     cores,
	     code, lda, oldCoding==TwoBitHaplo ? TwoBitGeno : UnknownSNPcoding);
}


void haplo2geno(SEXP Code, basic_options *opt) { //MoBPS  
  Long *info = GetInfo(Code);// may have been changed in createSNPmatrix
  haplo2geno(Align(Code, opt), info[SNPS], info[INDIVIDUALS], info[LDA],
	     (coding_type) info[CODING], (int) info[VARIANT],
	     GreaterZero(opt->cores));
}


void transform(unit_t *Old, Long OldSnps, Long OldIndiv, Long OldLDA,
	       coding_type OldCoding,  
	       int *snp_list, int snp_list_N,
	       int *indiv_list, int indiv_list_N,
	       usr_bool set1, basic_options *opt, 
	       bool VARIABLE_IS_NOT_USED ordered,
	       unit_t *Ans, Long AnsSnps, Long AnsIndiv, Long ldAns,
	       coding_type AnsCoding, int AnsVariant) {
 
  // idea is that for matrix multiplication the ordering
  // need not be kept if the reordering is the same for each column
  // then faster algorithms might be possible

 
  const bool OneHaploOnly = set1 != Nan;
  const coding_type selCoding = (!OneHaploOnly && isHaplo(OldCoding)) ||
    userHaplo(OldCoding) ? OneByteHaplo : OneByteGeno;
  const Long selSnps = snp_list_N > 0 ? AnsSnps : OldSnps,
    selIndiv = indiv_list_N > 0 ? AnsIndiv : OldIndiv,
    selLDA = SnpNatLDA(selSnps, selCoding),
    safety = 3 * MAX_LDA_BITALIGN / BitsPerUnit;
  const int main_ansvar = main_variant(AnsVariant);
  const Long selMemInUnits= totalMem(selSnps, selIndiv, selLDA,
					 selCoding, OneHaploOnly);

  //  printf("LDA = %d %d sel=%s:%d %d\n",  OldLDA, ldAns,  CODING_NAMES[selCoding], selSnps, selLDA );

  coding_type tmpCoding = OneByteGeno;
  Long tmpLDA = SnpNatLDA(OldSnps, tmpCoding);
  Long mem1ByteInUnits = totalMem(OldSnps, OldIndiv, tmpLDA, tmpCoding,
				      OneHaploOnly);
  coding_type oldCoding = OldCoding;
  int cores = GreaterZero(opt->cores);
  unit_t *old = Old,
    *neu = Ans;
  Long
    oldIndiv = OldIndiv,
    oldSnps = OldSnps;
  Long oldLDA = OldLDA;
  
  const Long tmp_units = mem1ByteInUnits;
  unit_t *tmp_unaligned = NULL,  *tmp = NULL,
    *sel_unaligned = NULL, *sel = NULL;

  bool explicit_selection_necessary = false, // dummy!
    any_selection = (snp_list_N > 0 || indiv_list_N > 0 ||
		     (OneHaploOnly && isHaplo(OldCoding)));

  assert(FirstUserHaploCoding == 20 && LastUserHaploCoding == 24);

  if (Old == Ans) {
    if (OldLDA != ldAns) 
      ERR0("transformation on place impossible (different storing length)");
    if (any_selection) ERR0("selection on place never possible");
    if (set1 != Nan) ERR0("haplo selection on place never possible");
  }

  // printf("%s -> %s\n", CODING_NAMES[OldCoding], CODING_NAMES[AnsCoding]);  
  if (AnsCoding == OldCoding && !any_selection) {
    Long cols = transposedHaplo(OldCoding) ? OldSnps : OldIndiv;
    if (doubledCols(OldCoding, false)) cols *= 2;
    Long bytes = MIN(ldAns, OldLDA) * BytesPerUnit;
   if (ldAns != OldLDA) {
     for (Long i=0; i<cols; i++)
	MEMCOPY(Ans + i * ldAns, Old + i * OldLDA, bytes);
    } else if (Old != Ans)
      MEMCOPY(Ans, Old, cols * bytes);
    goto ende;
  }

 
#define ALIGN_MALLOC(tmp, inUnits)					\
  tmp = (unit_t*) algn_generalL( (int*)(tmp##_unaligned = (unit_t*)	\
				      CALLOC((inUnits)+safety, BytesPerUnit)), \
			       MAX_LDA_BITALIGN )
 

  if (isHaplo(OldCoding)) {
    // printf("!!haplo\n");
     
    // Special cases, already coded directly
    if (!any_selection &&
	((OldCoding == OneBitHaplo && AnsCoding == OneBitGeno) ||
	 (OldCoding == TwoBitHaplo && ( AnsCoding == TwoBitGeno ||
					AnsCoding == ThreeBit)))) {
      haplo2geno(old, OldSnps, OldIndiv, oldLDA, OldCoding, main_ansvar, 
		 cores, Ans, ldAns, AnsCoding);
      goto ende;
    }

    if ((UnCompressedDoubleCols(OldCoding, true) ||
	 UnCompressedDoubleRows(OldCoding)) // range of SelectHaplo
	&&
	userHaplo(AnsCoding) && ! transposedHaplo(AnsCoding)// !OneByteHaplo
	&&
	UnCompressedDoubleCols(OldCoding, true) ==
	UnCompressedDoubleCols(AnsCoding, true) ) {
      if (isOneByteHaplo(OldCoding, false))
	SelectHaplo_1Byte_4Byte(Old, OldSnps, OldIndiv, OldLDA, OldCoding, 
			       snp_list, snp_list_N, indiv_list, indiv_list_N,
			       set1, Ans, AnsSnps, ldAns, AnsCoding);
      else 
	SelectHaplo_4Byte_4Byte(Old, OldSnps, OldIndiv, OldLDA, OldCoding,
			      snp_list, snp_list_N, indiv_list, indiv_list_N,
			      set1, Ans, AnsSnps, ldAns, AnsCoding);
      goto ende;   
    }

    if (OldCoding == OneByteHaplo) {
      //printf("EQUAL\n");
      
      tmp = Old;        // corresponds toAns Coding
      tmpCoding = OldCoding;
      tmpLDA = OldLDA;
    } else {
      // general case
      // printf("UNEQUAL userhaplo %d %lu\n",userHaplo(OldCoding), (Long) Ans);
      tmp = Ans;        // corresponds to AnsCoding
      tmpCoding = AnsCoding;
      tmpLDA = ldAns;
      if (userHaplo(OldCoding)) {
	// printf("UnCompressedDoubleCols(OldCoding, true)=%d\n", UnCompressedDoubleCols(OldCoding, true));
	if (UnCompressedDoubleCols(OldCoding, true) && any_selection) {
	  if (!isOneByteHaplo(tmpCoding, OneHaploOnly)) {
	    tmpCoding = selCoding;
	    tmpLDA = selLDA;
	    mem1ByteInUnits = selMemInUnits;
	    ALIGN_MALLOC(tmp, mem1ByteInUnits);
	  }
	  //	  printf("tmp = %lu\n", tmp);	  
	  SelectHaplo_4Byte_1Byte(Old,
				 OldSnps, OldIndiv, OldLDA, OldCoding,
				 snp_list, snp_list_N, indiv_list, indiv_list_N,
				 set1,
				  tmp, AnsSnps, tmpLDA,
				 tmpCoding);
	  explicit_selection_necessary = false;
	  set1 = set1 == False ? True : set1;
	  snp_list_N = 0;
	} else {
	  explicit_selection_necessary = snp_list_N > 0 || OneHaploOnly;
	  if (explicit_selection_necessary || !isCompressedHaplo(AnsCoding)) {
	    tmpCoding = selCoding;
	    tmpLDA = SnpNatLDA(OldSnps, tmpCoding); // not selSnps!
	    mem1ByteInUnits = totalMem(OldSnps, OldIndiv, tmpLDA, tmpCoding,
				       OneHaploOnly);
	    ALIGN_MALLOC(tmp, mem1ByteInUnits); ////// !!!!!!!!
	  }
	  codeHaplo_4Byte(Old, OldSnps, OldIndiv, OldLDA, OldCoding, 
			  indiv_list, indiv_list_N,
			  opt, tmp, AnsIndiv, tmpLDA, tmpCoding);
	  //
	  if (false && isOneByte(tmpCoding)) {
	    BUG;
	    _1Byte *t = (_1Byte*) tmp;
	    // printf("OldLDA=%d::\n", OldLDA); // R print
	    Long S = opt->bigendian ?  (1 + (OldSnps -1) / 4) * 4 : OldSnps;
	    for (Long jj=0; jj<S; jj++){
	      for (Long ii=0; ii<2 * OldIndiv; ii++)
	       PRINTF("%ld ", (Long) (t[jj + tmpLDA * ii * BytesPerUnit])); 
	      PRINTF("\n");
	    }
	  }
	}
	
      } else  { // !userHaplo(OldCoding)
	//	printf("trafo:not OneByteHaplo %lu\n", (Long) Old);
	explicit_selection_necessary = snp_list_N > 0;
 	if ((userHaplo(AnsCoding) || (OneHaploOnly && isOneSet2G(AnsCoding)))
	    && !explicit_selection_necessary) {
	  decodeHaplo_4Byte(Old, OldSnps, OldIndiv, OldLDA, OldCoding,
			   indiv_list, indiv_list_N, set1, cores,
			    tmp, AnsIndiv, tmpLDA, tmpCoding);
	} else {// !userHaplo(OldCoding) && (snp_list_N > 0 || ...)
	  tmpCoding = selCoding;
	  tmpLDA = SnpNatLDA(OldSnps, tmpCoding); // not selSnps!
	  mem1ByteInUnits = totalMem(OldSnps, OldIndiv, tmpLDA, tmpCoding,
				     OneHaploOnly); 	  
	  ALIGN_MALLOC(tmp, mem1ByteInUnits);	////// !!!!!!!!  
	  decodeHaplo_1Byte(Old, OldSnps, OldIndiv, OldLDA, OldCoding,
			    indiv_list, indiv_list_N, set1, cores,
			     tmp, AnsIndiv, tmpLDA, tmpCoding);
	}
	set1 = set1 == False ? True : set1;
      } // !userHaplo(OldCoding)
      indiv_list_N = 0;
    } // OldCoding != OneByteHaplo
    
    if (tmp == Ans) {
      assert(!explicit_selection_necessary);
      goto ende;
    }
    
    if (!OneHaploOnly) { // OneHaploOnly may or may not be true
      assert(tmpCoding == OneByteHaplo);
  
      if (isFourByteHaplo(AnsCoding, OneHaploOnly)) {
	SelectHaplo_1Byte_4Byte(tmp,
			       OldSnps, selIndiv, tmpLDA, tmpCoding,
			       snp_list, snp_list_N, indiv_list, indiv_list_N,
			       set1, Ans, selSnps, selLDA, selCoding);
	goto ende;
      }

      if (explicit_selection_necessary) { // i.e., snp_list_N > 0
	if (isOneByteHaplo(AnsCoding, OneByteHaplo)) sel = Ans;
	else ALIGN_MALLOC(sel, selMemInUnits);
	
	SelectHaplo_1Byte_1Byte(tmp,
				OldSnps, selIndiv, tmpLDA, tmpCoding,
				snp_list, snp_list_N, indiv_list, indiv_list_N,
				set1,
				sel, selSnps,  selLDA,
				selCoding);
	if (sel == Ans) goto ende;
	assert(set1 == Nan);
 	snp_list_N = indiv_list_N = 0;
        explicit_selection_necessary = false;
	FREE(tmp_unaligned);
	tmp_unaligned = sel_unaligned;
	tmp = sel;
	sel_unaligned = NULL;
	tmpCoding = selCoding;
	tmpLDA = selLDA;	
      }
      
      // printf("here %d %d %s -> %s -> %s\n", tmp == Ans, tmp == Old, CODING_NAMES[OldCoding], CODING_NAMES[tmpCoding],  CODING_NAMES[AnsCoding]);
      
      if (isHaplo(AnsCoding)) {
	if (userHaplo(AnsCoding))
	  decodeHaplo_4Byte(tmp, selSnps, selIndiv, tmpLDA,  tmpCoding, 
			   indiv_list, indiv_list_N, set1, cores,
			   Ans, AnsIndiv, ldAns, AnsCoding);
	else codeHaplo_1Byte(tmp, selSnps, selIndiv, tmpLDA, tmpCoding,
			     indiv_list, indiv_list_N, opt,
			     Ans, AnsIndiv, ldAns, AnsCoding);
	goto ende;
      }

      // halplo --> geno !
      OneByteHaplo2geno64(tmp, selSnps, selIndiv, tmpLDA, cores, tmp, tmpLDA);
      oldSnps = selSnps;
    } else { // OneHaploOnly      
      if (set1==False)
	MEMCOPY(tmp, tmp + tmpLDA * selIndiv, tmpLDA * selIndiv * BytesPerUnit);
      explicit_selection_necessary = snp_list_N > 0; 
      oldSnps = explicit_selection_necessary ? OldSnps : selSnps; // note:
      // snp_list_N may have changed value since definition of selSnps!
    }
    
    old = tmp;
    oldCoding = OneByteGeno;
    oldLDA = SnpNatLDA(oldSnps, oldCoding);
    oldIndiv = selIndiv;
    any_selection = explicit_selection_necessary;
  }
  
  //////////////////// END HAPLO //////////////////////


  //////////////////// START GENO /////////////////////
  assert(!isHaplo(oldCoding));
  assert(sel_unaligned==NULL); 
  if (isHaplo(AnsCoding))
    ERR0("Genotype matrix cannot be transformed into a haplotype matrix");

  
  if (oldCoding == FourByteGeno && AnsCoding == FourByteGeno) {
    // either selection or different lda
    SelectGeno_4Byte_4Byte(old, oldSnps, oldIndiv, oldLDA,
			   snp_list, snp_list_N, indiv_list, indiv_list_N,
			   Ans, ldAns // eigene Zeile wg findallM
			 );
    goto ende;
  }
  
  if (any_selection) {
    // printf("%s --> %s\n", CODING_NAMES[oldCoding], CODING_NAMES[AnsCoding]);
    
    if (AnsCoding == OneByteGeno || AnsCoding == FourByteGeno
	) sel = Ans;
    else ALIGN_MALLOC(sel, selMemInUnits);
    if (oldCoding == FourByteGeno) {
      assert(tmp_unaligned==NULL); // as tmp_unaligned is set only if Haplo
      SelectGeno_4Byte_1Byte(old, oldSnps, oldIndiv, oldLDA,
			     snp_list, snp_list_N, indiv_list, indiv_list_N,
			     sel, selLDA);
    } else {
      coding_type TempCoding= oldCoding;
      Long TempLDA = oldLDA;
      unit_t *Temp_unaligned = tmp_unaligned,
	*Temp = old;
      // printf("TempNULL = %d %s \n", Temp_unaligned==NULL, CODING_NAMES[oldCoding] );
      if (oldCoding == OneByteGeno) {
	assert((tmp_unaligned == NULL) xor (old == sel)); // coming from Haplo?
	tmp_unaligned = NULL;
      } else {
	assert(tmp_unaligned == NULL);
	if (!isCompressed(oldCoding)) BUG;
	TempCoding = OneByteGeno;
	TempLDA = SnpNatLDA(oldSnps, TempCoding);
	Long memInUnitsTemp = totalMem(oldSnps, oldIndiv, TempLDA,
				       TempCoding, false);
	ALIGN_MALLOC(Temp, memInUnitsTemp);
	assert(Temp_unaligned != NULL);
	matrix_get_1Byte(old, oldSnps, oldIndiv, oldLDA, oldCoding, 
			 cores, Temp, TempLDA);
      }
      if (AnsCoding == FourByteGeno) {
	SelectGeno_1Byte_4Byte(Temp, oldSnps, oldIndiv, TempLDA,
			      snp_list, snp_list_N, indiv_list, indiv_list_N,
			      sel, ldAns);
	FREE(Temp_unaligned);
 	goto ende;
      } else
	SelectGeno_1Byte_1Byte(Temp, oldSnps, oldIndiv, TempLDA,
			       snp_list, snp_list_N, indiv_list, indiv_list_N,
			       sel,  selLDA);
      FREE(Temp_unaligned);
    } // oldCoding != FourByteGeno
    if (sel == Ans) goto ende;    
    old = sel;
    oldSnps = selSnps;
    oldIndiv = selIndiv;
    oldCoding = OneByteGeno;
    oldLDA = selLDA;
  }

  ///// no selection ////

  if (anyByteGeno(AnsCoding)) {
    //    printf("ans = 4Byte (%s->%s)\n", CODING_NAMES[oldCoding], CODING_NAMES[AnsCoding]);
    if (AnsCoding == FourByteGeno)
      matrix_get_4Byte(old, oldSnps, oldIndiv, oldLDA, oldCoding, cores,
		      Ans, SnpNatLDA(oldSnps, AnsCoding));
    else matrix_get_1Byte(old, oldSnps, oldIndiv, oldLDA, oldCoding, 
			  cores, Ans, SnpNatLDA(oldSnps, AnsCoding));
    goto ende;
  }

  // printf("reaching switch OldCoding\n");
  assert(tmp_unaligned == NULL);
  
  switch(oldCoding) {
  case OneBitGeno:{
    switch (AnsCoding) {
    case TwoBitGeno:
      switch (main_ansvar) {
      case 512 : case 256 : case 128 :
	if (opt->bigendian) {
	  if (old == Ans) {
	    ALIGN_MALLOC(tmp, tmp_units); 
	    neu = tmp;
	  }
	  trafo1Geno2Geno128(old, oldSnps, oldIndiv, oldLDA, neu, ldAns);
	  goto copy;
 	}
	FALLTHROUGH_OK;
      case 64: case 32: 
	MEMCOPY(Ans, old, oldLDA * oldIndiv * BytesPerUnit);
	old = Ans;
	goto via_1Byte; 	
      default: BUG;
      }
      break;
    case ThreeBit:
      if (ldAns < oldLDA) ERR0("transformation to 3bit not possible");
      MEMCOPY(Ans, old, oldLDA * oldIndiv * BytesPerUnit);  
      old = Ans;
      goto via_1Byte;
    default: BUG;
    }
  }
    break;
    
  case Plink:
  case TwoBitGeno:{
    switch (AnsCoding) {
    case OneBitGeno:
      if (oldCoding == Plink) BUG;
      switch (main_ansvar) {
      case 512 : case 256 : case 128 :
	if (opt->bigendian) {
	  if (old == Ans) {
	    ALIGN_MALLOC(tmp, tmp_units);    
	    neu = tmp;
	  }
	  trafo2Geno1Geno128(old, oldSnps, oldIndiv, oldLDA, neu, ldAns);
	  goto copy;
	}
	FALLTHROUGH_OK;
      case 64: case 32:
	MEMCOPY(Ans, old, oldLDA * oldIndiv * BytesPerUnit);
	old = Ans;
	goto via_1Byte;
      default: BUG;
      }
      break;
    case FiveCodes:
      //  printf("bigendian = %d\n", opt->bigendian);
      if (opt->bigendian) BUG;
      if (old == Ans) BUG;
      switch (main_ansvar) {
      case 512 : case 256 : case 128 :
	trafo2Geno5codes256(old, oldSnps, oldIndiv, oldLDA, oldCoding,
			    cores,
			    neu, ldAns);
	goto copy;
      case 64: case 32:
	trafo2Geno5codes32(old, oldSnps, oldIndiv, oldLDA, oldCoding,
			    cores, neu, ldAns);
	goto copy;
      default: BUG;
      }
      break;
    case FiveCodesTransposed: 
      if (opt->bigendian) BUG;
      if (old == Ans) BUG;
      switch (main_ansvar) {
      case 512 : case 256 : case 128 :
	trafo2Geno5codestrans256(old, oldSnps, oldIndiv, oldLDA, oldCoding,
				 cores, 
				 neu, ldAns);
	goto copy;
      case 64: case 32:
	trafo2Geno5codestrans32(old, oldSnps, oldIndiv, oldLDA, oldCoding,
				cores,
				neu, ldAns);
	goto copy;
      default: BUG;
      }
      break;
    case ThreeBit:
      if (oldCoding == Plink) BUG;
      if (ldAns < oldLDA) ERR0("transformation to 3bit not possible");
      MEMCOPY(Ans, old, oldLDA * oldIndiv * BytesPerUnit);
      old = Ans;
      goto via_1Byte;
    default: BUG;
    }
  }
    break;

  case ThreeBit :
  case OneByteGeno:
    //    printf("onebyte via 1Byte");
    goto via_1Byte;
    
  case FourByteGeno : {
    coding_t
      coding_uint = matrix_coding_4Byte(AnsCoding, AnsVariant, opt->bigendian);
    
    coding_uint(old , oldLDA,
		0, oldSnps, 0, oldIndiv, 
		cores, NULL, Ans, oldIndiv, ldAns);
    goto ende;
  }
  
  case FourBit :
  case TwoByte :
  case FiveCodes :
  case FiveCodesTransposed:  BUG; // not programmed
  case FourByte_Tr: BUG; // programmed, but should not appear
	   
  case DotFile:
  case FileDot:
  case CorrespondingGeno :
  case UnknownSNPcoding :
  default: // programming error
    //  printf("%s\n", CODING_NAMES[oldCoding]);
    BUG;
  }

 via_1Byte :
  tmpCoding = OneByteGeno;
  if (oldCoding == tmpCoding) {
    tmp = old;
    tmpLDA = oldLDA;
  } else {
    tmpLDA = SnpNatLDA(oldSnps,tmpCoding);
    if (tmp_unaligned == NULL) ALIGN_MALLOC(tmp, tmp_units);
    matrix_get_1Byte(old, oldSnps, oldIndiv, oldLDA, oldCoding,
		     cores,  tmp, tmpLDA);
  }
  if (tmp == Ans) BUG;
  coding_t matrix_coding;
  matrix_coding = matrix_coding_1Byte(AnsCoding, AnsVariant, opt->bigendian);
  matrix_coding(tmp, tmpLDA,
		0, oldSnps, 0, oldIndiv,
		cores, NULL,
		Ans, oldIndiv, ldAns);
  assert(neu == Ans);
  
 copy :  
  if (neu != Ans) {
    Long memInUnits = totalMem(oldSnps, oldIndiv, ldAns, AnsCoding, set1);
    MEMCOPY(Ans, neu, memInUnits * BytesPerUnit);
  }
			
 ende :
  FREE(sel_unaligned);
  FREE(tmp_unaligned);
}



void transform(unit_t *Old, Long snps, Long individuals, Long OldLDA,
	       coding_type OldCoding, 
	       int *snp_list, int snp_list_N,
	       int *indiv_list, int indiv_list_N,
	       usr_bool set1, basic_options *opt,
	       bool ordered, SEXP Code, Long AnsSnps, Long AnsIndiv) {
#if defined compatibility_to_R_h
  PROTECT(Code);
#endif  
  Long *info = GetInfo(Code);// may have been changed in createSNPmatrix
  unit_t *code = Align(Code, opt);
  ASSERT_LITTLE_ENDIAN;
  //	printf("intermediate trafo %lu\n", (Long) Old);
  
  transform(Old, snps, individuals, OldLDA, OldCoding, 
	    snp_list, snp_list_N, indiv_list, indiv_list_N,
	    set1, opt, ordered,
	    code, AnsSnps, AnsIndiv, info[LDA],
	    (coding_type) info[CODING], (int) info[VARIANT]);
  //printf("done %u %u %u %u\n", code[0], code[1], code[2], code[3]);
  
  UNPROTECT(1);
}




void transposeIntern(unit_t * Code, Long snps, Long individuals,
		     coding_type coding,  Long LDAbitalign,
		     bool VARIABLE_IS_NOT_USED bigendian,
		     unit_t *Ans) {
  if (Code == Ans) ERR0("traasposing cannot be performed on place");
  if (coding == FourByteGeno) {
    Long lda = GetLDA(snps, individuals, coding, LDAbitalign);
    for (Long j=0; j<individuals; j++) {
      Long j0 = j * lda; // == j * snps
      for (Long i=0; i<snps; i++) {
	Ans[i * snps + j] = Code[j0 + i];
      }
    }
    return;
  }

  bool twobit = coding == TwoBitGeno || coding == TwoBitHaplo;
  bool onebit = coding == OneBitGeno || coding == OneBitHaplo;

  if (!(onebit || twobit)) {
    ERR1("transposing for %.50s not programmed yet.", CODING_NAMES[coding]);
    return;
  }
  //if (bigendian)ERR0("transposing for %.50s not programmed for big endian.");
  
  const int
    bitPerCode = 1 + twobit,
    codePerByte_default = BitsPerByte / bitPerCode,
    codePerUnit_default = BytesPerUnit * codePerByte_default;
  Long
    Lda = GetLDA(snps, individuals, coding, LDAbitalign),
    real_units[2] = {snps / codePerUnit_default, 1};
  assert((Long) snps - real_units[0] * codePerUnit_default >= 0);
  Long
    real_CpU[2] = {codePerUnit_default,
                   snps - real_units[0] * codePerUnit_default},
    end_ll = snps % codePerUnit_default > 0;
  Long
    TldaBytes = GetLDA(individuals, snps, coding, LDAbitalign) * BytesPerUnit;
  unsigned char *ans = (unsigned char *) Ans;

#define N6E1 0x03

  if (onebit) {
    ERR0("t() for OneBit Not programmed yet");
    return;
  }
  
  assert(twobit);
   unit_t *zeros = NULL;
   int rest_indiv = (int) (individuals % codePerByte_default);
   Long   indiv = individuals  - rest_indiv;
   int codePerByte = codePerByte_default;
   unit_t
    *Code0 = Code,
    *Code1 = Code0 + Lda,
    *Code2 = Code1 + Lda,
    *Code3 = Code2 + Lda;//, *Code4, *Code5, *Code6, *Code7;
  unsigned char *ans1 = ans;
  for (Long k=0; k<=(rest_indiv > 0); k++) {
    //    printf("k=%d %d\n", k, rest_indiv > 0);
    if (k) {
      zeros = (unit_t*) CALLOC(real_units[0] + 1, BytesPerUnit); //real_units=[0]=0 moeglich!
      Code0 = Code + Lda * indiv ;
      Code1 = codePerByte > 1 ? Code0 + Lda : zeros;
      Code2 = codePerByte > 2 ? Code1 + Lda :  zeros;
      Code3 = zeros;
      ans1 = ans + indiv / codePerByte;
      indiv = rest_indiv > 0;
      codePerByte = rest_indiv;
      //printf("zeros %d %d %d\n", Code1 == zeros, Code2 == zeros,  Code3== zeros);
    }
  
    for (Long j = 0; j<indiv; j+=codePerByte, ans1++) {
      unit_t 
	*S00 = Code0 + j * Lda,
	*S11 = Code1 + j * Lda,
	*S22 = Code2 + j * Lda,
	*S33 = Code3 + j * Lda;
      //printf("%u %u %u %u bPc=%d %d\n", S0, S1, S2, S3, bitPerCode, N6E1);
      unsigned char *ans2 = ans1;
      for (Long ll = 0; ll <= end_ll; ll++) { // steuert S N P S
	Long units = real_units[ll],
	  CpU = real_CpU[ll];
	//	printf("k=%d j=%d ll=%d units=%d CpU=%d rest=%d indiv=%d, codePerB=%d\n", k, j, ll, units, CpU, rest_indiv, indiv, codePerByte);
	
	for (Long m=0; m<units; m++, S00++, S11++, S22++, S33++) {
	  unit_t S0 = *S00,
	    S1 = *S11,
	    S2 = *S22,
	    S3 = *S33;
	  for (int i=0; i<CpU; i++, ans2 += TldaBytes) {
	    *ans2 = (S0 & N6E1) |
	      ((S1 & N6E1) << bitPerCode) |
	      ((S2 & N6E1) << (2 * bitPerCode)) |
	      ((S3 & N6E1) << (3 * bitPerCode));
	    S0 >>= bitPerCode;
	    S1 >>= bitPerCode;
	    S2 >>= bitPerCode;
	    S3 >>= bitPerCode;
	    //	    printf("%u:%u,%u,%u,%u ", *ans2, S0, S1,S2,S3);
	  } // i
	} // m
      } // ll
    } // j
  } // k
  FREE(zeros);
}


SEXP transpose(SEXP SxI, option_type *global, utilsoption_type *utils) { 
  basic_options *opt = &(utils->basic);
  int efficient =  opt->efficient;
  Long *info = GetInfo(SxI);
  Long
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    LDAbitalign = info[LDABITALIGN];
  coding_type coding = (coding_type) info[CODING];
  int variant = check_variant(coding, (int) info[VARIANT], efficient);
  unit_t *C = Align(SxI, opt);

  SEXP Ans = PROTECT(createSNPmatrix(individuals, // snp and indiv swapped!
				     snps, coding, variant, 0, global, utils));
  Long *infoAns = GetInfo(Ans);
  Long LDAbitAns = infoAns[LDABITALIGN];
  coding_type codingAns = (coding_type) infoAns[CODING];
  if (LDAbitalign != LDAbitAns)
    ERR0("transposition not possible when LDAbitalign has been changed.");
  if (coding != codingAns)
    ERR0("transposition not possible when snpcoding has been changed.");
  unit_t *A =  Align(Ans, opt);
   
  transposeIntern(C, snps, individuals, coding, 
		  LDAbitalign, info[BIGENDIAN], A);
  UNPROTECT(1);
  return Ans;
}

 
  
SEXP Transform(SEXP SxI, unit_t* SxIint, coding_type codingInfo,
	       int* selSnps, int lenSnps, int *selIndiv, int lenIndiv,
	       option_type *global, utilsoption_type *utils) {
  //  printf("entering Transform %s\n", CODING_NAMES[codingInfo]);
#if defined compatibility_to_C_h
  if (SxIint != NULL) {
    // SxI contains always all relevant information
    // SxI$code might be real-valued in R
    // Transform expects however integer values -> SxIint
    ERR0("In standalone mode, 'Transform' may be called only with SxIint=NULL.");
  }
#endif
  		 
   basic_options *opt = &(utils->basic);
  const usr_bool set1 =  global->genetics.haplo_set1;
  const bool as_is = global->genetics.interprete_coding_as_is,
    OneHaploOnly = set1 != Nan,
    is_file = TYPEOF(SxI) == STRSXP;
  //  tuning_options *tuning = &(global->tuning);
  
  Long *info = GetInfo(SxI, true);
  
  const bool
    only_complete = codingInfo == UnknownSNPcoding,
    no_select = lenSnps == 0 && lenIndiv == 0,
    copy = !(no_select &&
	     (only_complete || (info != NULL && info[CODING] == codingInfo)));
  int efficient =  opt->efficient;
 
  coding_type
    adjustedCodingInfo = (codingInfo==AutoCoding ||(only_complete && info==NULL)
			  ? global->genetics.coding
			  : only_complete ? (coding_type) info[CODING]
			  : codingInfo); 
 
  assert(LastGenoCoding == 13);
  if (adjustedCodingInfo == FourBit || adjustedCodingInfo == TwoByte)
    ERR0("'FourBit' and 'TwoByte' are not programmed yet.");  
  if (!(isNarrowGeno(adjustedCodingInfo) || isHaplo(adjustedCodingInfo))) BUG;
  
  SEXP
    filecoding = R_NilValue,
    filename = R_NilValue,
    Code = R_NilValue;
  int n_protect = 0;  
  int variant = //  info != NULL ? info[VARIANT] :
    global->tuning.variant;
  unit_t *code = NULL;
  Long
    AnsSnps = NA_LONG,
    AnsIndiv = NA_LONG,
    oldLDA = NA_LONG, 
    OldIndiv = NA_LONG,
    OldSnps = NA_LONG;

  //  printf("info[%d] variant %d %d -> %d\n",  info != NULL, (int) info[VARIANT] , global->tuning.variant, variant );
 
  coding_type oldCoding = UnknownSNPcoding,
    ansCoding = UnknownSNPcoding;

  
  if (is_file) {
    // printf("check has_sexp_coding !!! \n");
    ansCoding = adjustedCodingInfo;
    if (isCompressed(ansCoding) && !as_is)
      ansCoding = swapHaploGeno(ansCoding, false);
    filename = PROTECT(COPY_SEXP_STRING(SxI)); n_protect++;
    bool early_ret = (no_select && has_sexp_coding(ansCoding)) || only_complete;
    coding_type tmpCoding = early_ret ? ansCoding :
      isHaplo(ansCoding) ? OneByteHaplo : OneByteGeno;
    int tmpVariant = check_variant(tmpCoding, variant, efficient);
    
    // printf("is_file************* ret=%d %s --> %s (as_is=%d)\n", early_ret, CODING_NAMES[tmpCoding], CODING_NAMES[ansCoding], as_is);

    PROTECT(Code = file_intern(SxI, // ToDo: hier gleich auswahl!
			       // dann muss tmpCoding nicht auf OneByteGeno
			       // gesetzt werden.
			       tmpCoding, tmpVariant,
			       global, utils, only_complete, R_NilValue));

    //printf("************* early_ret=%d complete=%d\n", early_ret, only_complete);
    
    n_protect++;
    if (early_ret) goto ende;
    FREE_SEXP(&SxI);
    SxI = Code;
    info = GetInfoUnchecked(SxI); // OK
    oldCoding = (coding_type) info[CODING];
    assert(info != NULL && ansCoding != CorrespondingGeno);
  } else {
    //   printf("is matrix ***** %s %s %s\n", CODING_NAMES[oldCoding], CODING_NAMES[adjustedCodingInfo], CODING_NAMES[codingInfo]);
    if (copy) {
      SEXP tmp = getAttribPointer(SxI, Filecoding);
      if (tmp != R_NilValue) {
	filecoding = PROTECT(COPY_SEXP_STRING(tmp)); n_protect++;
      }
      tmp = getAttribPointer(SxI, Filename);
      if (tmp != R_NilValue) {
	filename = PROTECT(COPY_SEXP_STRING(tmp)); n_protect++;
      }
    }

    if (info == NULL) { //  from User To Compressed
      oldCoding = (coding_type) adjustedCodingInfo;
      if (oldCoding==AutoCoding) oldCoding = FourByteGeno;
      //    OldVariant = 32;
      ansCoding = global->genetics.coding;
      if (!as_is) ansCoding = swapHaploGeno(ansCoding, isHaplo(oldCoding));
    } else {
      oldCoding = (coding_type) info[CODING];
      //     OldVariant = (int) info[VARIANT];      
      ansCoding = adjustedCodingInfo;
      ASSERT_LITTLE_ENDIAN;
    }

    if (ansCoding == CorrespondingGeno)
      ansCoding = swapHaploGeno(oldCoding, false);
  }
  
  //  printf("******** %s coding: task: decoded=%s; %s glbl= %s -> %s  isHa=%d\n", info == NULL ? "InfoIsNull" : "Internal",   CODING_NAMES[adjustedCodingInfo],CODING_NAMES[oldCoding],	    CODING_NAMES[ global->genetics.coding],  CODING_NAMES[ansCoding],	   isHaplo(oldCoding));
  

  if (info == NULL) {
    //   printf("NO INFO\n");
   code = SxIint;

    //   printf("code = %d %d  %d %d\n", code[0], code[1], code[2], code[3]);
    
    if (LENGTH(SxI) == 0) ERR0("'SxI' has length 0.");
    if (isMatrix(SxI)) {
      OldIndiv = ncols(SxI); // OK
      OldSnps = nrows(SxI); // OK
    } else {
      OldIndiv = 1;
      OldSnps = LENGTH(SxI);
    }
    if (isOneByte(oldCoding)) OldSnps *= 4;
    if (userHaplo(oldCoding)) {
      if (transposedHaplo(oldCoding)) {
	Long tmp = OldSnps;
	OldSnps = OldIndiv;
	OldIndiv = tmp;
      }
      if (doubledSnps(oldCoding)) {
	if (OldSnps % 2 != 0) ERR0("snps are not doubled");
	OldSnps /= 2;
      } else {
	if (OldIndiv % 2 != 0) ERR0("snps are not doubled");
	OldIndiv /= 2;
      }
    }

    //  printf("snps=%d, %s/%s: %s -> %s\n", OldSnps, CODING_NAMES[codingInfo],CODING_NAMES[adjustedCodingInfo],CODING_NAMES[oldCoding],  CODING_NAMES[ansCoding]);
    
    oldLDA = GetNatLDA(OldSnps, OldIndiv, oldCoding);
    
    bool check_included = ansCoding == TwoBitGeno || ansCoding == OneBitGeno
      || userHaplo(oldCoding);
    if (!opt->skipchecks && !check_included) {
      assert(oldCoding == FourByteGeno); 
      bool ok = true;
      unit_t *SxIU = (unit_t*) SxIint;
      Long total = OldSnps * OldIndiv;
      for (Long i=0; i<total; i++) ok &= SxIU[i] <= 2;
      if (!ok)  ERR0("SNP matrix nay have only the values 0,1,2");
    }
  } else {
    //   printf(" INFO\n");
    assert(info != NULL);
    code = Align(SxI, opt); // nach info

    //   printf("alignment difference %ld\n", (Long) code - (Long) INTEGER(SxI));
      
    OldIndiv = info[INDIVIDUALS];
    OldSnps =  info[SNPS];
    oldLDA = info[LDA];
  }

  // printf("preprocessing\n");
  
  // preprocessing the selection options
  AnsSnps = OldSnps;
  AnsIndiv = OldIndiv;
  if (lenSnps > 0) {
    // comfortable to give "to the end" by 1:n where n is just
    // as large as than the genuine size!
    // So, at the end, values exceeding the genine size are expected:
    while (lenSnps >= 0 && selSnps[lenSnps -1] >= OldSnps) lenSnps--;    
    AnsSnps = 0;
    for (Long i=0; i<lenSnps; i++) AnsSnps += selSnps[i] < OldSnps;
    if (AnsSnps == 0) ERR0("SNP selection out of range");
  }
  if (lenIndiv > 0) {
    while (lenIndiv >= 0 && selIndiv[lenIndiv -1] >= OldIndiv) lenIndiv--;
    AnsIndiv = 0;
    for (Long i=0; i<lenIndiv; i++) AnsIndiv += selIndiv[i] < OldIndiv; 
    if (AnsIndiv == 0) ERR0("individual selection out of range");
  }

  //  printf("check_variant %s %d\n", CODING_NAMES[ansCoding], variant);
  variant = check_variant(ansCoding, variant, efficient);
  //  printf("variant checked %d\n", variant);
  PROTECT(Code = CompleteCodeVector(AnsSnps, AnsIndiv, ansCoding,
				    variant, 0, OneHaploOnly,
				    global, utils, NULL, NULL,
				    copy ? R_NilValue : SxI));
  // printf("code vector created\n");

  n_protect++;
  // createSNPmatrix sets INFO and Class (if possible)
  // But Filename/Coding is set at the very end, since the part here might be
  // skipped


  assert(OldSnps != NA_LONG &&  OldIndiv!=NA_LONG && oldLDA !=NA_LONG);
  assert(code != NULL);

  if (false && is32BitAligned(oldCoding)) {
    BUG;
    int endI = (1 + isHaplo(oldCoding) * (int) MIN(OldIndiv, 10)  *
		UnCompressedDoubleCols(oldCoding, false)),
      endS = (int) MIN(OldSnps, 15);
    for (int jj=0; jj<endS; jj++){
      for (int ii=0; ii<endI; ii++)
	printf("%d ",(int)  code[jj + oldLDA * ii]);     // Rprintf
      if (endI < OldIndiv) printf("...[%lu]\n", OldIndiv); // Rprintf
      printf("\n"); // Rprintf
    }
    if (endS < OldSnps) printf("... [%lu]", OldSnps); // Rprintf
  }

  // printINFO(Code); BUG;
 

  transform(code, OldSnps, OldIndiv, oldLDA, oldCoding, 
	    selSnps, lenSnps, selIndiv, lenIndiv,
	    set1, opt, true,
	    Code, AnsSnps, AnsIndiv);

  //  printINFO(SxI);
  // printINFO(Code);
  //
  
  SEXP next;
  next = getAttribPointer(Code, Next);
  if (next != R_NilValue) {
    // printINFO(Code);    printINFO(next);
    SEXP N = PROTECT(Transform(SxI, SxIint, (coding_type) GetInfo(next)[CODING],
			       selSnps, lenSnps, selIndiv, lenIndiv,
			       global, utils));
    n_protect++;

    setAttrib(Code, Next, N);
  }

  
 ende:
 
  if (filename != R_NilValue) setAttrib(Code, Filename, filename); // install !
  if (filecoding != R_NilValue) setAttrib(Code, Filecoding, filecoding);
  copySEXPattrib(Code, SxI, false);
  if (info != NULL) {
    Long *infoCode = GetInfo(Code);
    infoCode[MISSINGS] = info[MISSINGS]; 
  }
  if (n_protect) { UNPROTECT(n_protect); }

  // printf("ende transformUint\n");
  return Code;  
}
