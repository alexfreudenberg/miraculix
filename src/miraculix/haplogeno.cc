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


// _MM_ALIGN16 Uint ccc;

#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "Template.h"
#include "utils_miraculix.h"
#include "extern.h"
#include "mmagpu.h"
#include "transform.h"


#include "MXinfo.h"
#include "intrinsicsCheck.h"
#include "intrinsics_specific.h"
#include "Haplo.h"
#include "Files.h"
#include "1bit.h"
#include "bitBase.h"
#include "2bit.h"
#include "3bit.h"
#include "5codes.h"
#include "OneByte.h"
#include "4Byte.h"
#include "plink.h"
#include "Vector.matrix.h"
#include "kleinkram.h"


void allInfo(Long *info) { // ok
  for (int J=0; J<=INFO_GENUINELY_LAST; J++) {
    PRINTF("%s=", INFO_NAMES[J]);
    switch(J) {
    case CODING : PRINTF("%s", CODING_NAMES[info[J]]); break;
    case VARIANT : PRINTF("implementation %ld", info[J]); break;
    case ADDR : PRINTF(" %lu", info[ADDR]); break;
    case ALIGNADDR :
      PRINTF("%lu",  info[ALIGNADDR]); break;
    case MEMinUNITS : PRINTF("%lu", info[MEMinUNITS]); break; // OK
    case ALIGNEDUNITS : PRINTF("%lu", info[ALIGNEDUNITS]); break;// OK
    default:
      if (info[J] == NA_LONG) PRINTF("NA"); else PRINTF("%lu", info[J]);
    }
    PRINTF("\n");
  }
}
void allInfo(SEXP M) { // ok
  Long *info = GetInfoUnchecked(M);
  if (info == NULL) BUG;
  allInfo(info); // ok
}
 

Long GetBytesPerBlock(coding_type coding, int variant) {
  variant = main_variant(variant);
  if (variant <= 1024) return variant >> 3;
  if (variant == VARIANT_GPU)
    return coding == FourByteGeno ? 4 : (256 >> 3);
  if (variant == VARIANT_R) return 4; 
  BUG; 
  return 0;
}

Long GetCodesPerBlock(coding_type coding, int variant) {
#define twobits 2  
  assert(FirstUserHaploCoding == 20 && LastUserHaploCoding == 24);
  Long BpB = GetBytesPerBlock(coding, variant);
  Long codes;
  switch (coding) {
  case OneBitGeno: case OneBitHaplo: codes = BpB * BitsPerByte; break;
  case TwoBitGeno: case TwoBitHaplo: case Plink: case OrigPlink:
    codes = BpB * BitsPerByte / twobits; break;
  case ThreeBit : codes = CodesPerBlock3(); break;
  case FourBit : BUG; break;
  case OneByteGeno :  codes = BpB; break;
  case TwoByte : BUG; break;
  case FiveCodesTransposed :
  case FiveCodes : codes = 5 * BpB; break;
  case FourByteGeno : case FourByte_Tr : codes = BpB / 4; break;
  case OneByteHaplo: codes = BpB; break;
  case FourByteHaplo: codes = BpB / 4; break;
  case FourByteSingleBit: codes = BpB / 4; break;
  case EightByteHaplo: codes = BpB / 8; break;
  case FourByteSingleBit_Tr: codes = BpB / 8; break;
  case EightByteHaplo_Tr: codes = BpB / 4; break;    
  case FileDot:
  case DotFile : codes = 100; break; // arbitrary number
  default : BUG;
  }
  return codes;
}

Long GetCodesPerUnit(coding_type coding) { // OK
  // ACHTUNG: EightHaplo braucht 2 U n i t  per Code
  int variant = 64;
  Long zaehler = GetCodesPerBlock(coding, variant) * BytesPerUnit, // OK
    nenner = GetBytesPerBlock(coding, variant);
  //printf("zaehler /nenneer %ld %ld %ld\n", zaehler, nenner, GetCodesPerBlock(coding, variant));
  if (zaehler % nenner) BUG;
  return  zaehler / nenner; // OK
}

Long GetBlocks(Long snps, coding_type coding, int variant) {
  return snps2entities(snps, GetCodesPerBlock(coding, variant));
}


Long GetLDA(Long snps, Long individuals, coding_type coding, Long LDAbitalign) {
  // LDA in Einheiten von unit_t
  assert(snps == 1 || LDAbitalign > 0);
  if (LDAbitalign % BitsPerUnit != 0) // OK
    ERR1("bit alignment not multiple of %d", BitsPerUnit); // OK
  Long lda;
  switch(coding) {
  case OneBitGeno : case OneBitHaplo : lda = Lda1Bit(snps, LDAbitalign); break;
  case TwoBitGeno : case TwoBitHaplo : lda = Lda2Bit(snps, LDAbitalign); break;
  case TwoBitGenoTransposed : lda = Lda2Bit(individuals, LDAbitalign); break;
  case ThreeBit : lda = Lda3(snps, LDAbitalign); break;
  case FourBit : BUG;
  case OneByteGeno :
  case OneByteHaplo : lda = LdaOneByte(snps, LDAbitalign); break;
  case TwoByte: BUG;
  case FourByteGeno : lda = LdaPlain(snps, LDAbitalign); break;
  case FourByte_Tr : lda = LdaPlain(individuals, LDAbitalign); break;
  case FourByteHaplo : lda = LdaPlain(snps, LDAbitalign); break;
  case Plink : lda = LdaPlink(snps, LDAbitalign); break;
  case PlinkTransposed : lda = LdaPlink(individuals, LDAbitalign); break;
  case OrigPlink : lda = LdaOrigPlink(snps, LDAbitalign); break; // ldaInByte
  case OrigPlinkTransposed : lda = LdaOrigPlink(individuals, LDAbitalign);break;
  case FiveCodes : lda = LdaFiveCodes(snps, LDAbitalign); break;
  case FiveCodesTransposed: lda = LdaFiveCodes(individuals, LDAbitalign);
    break;
  case FourByteSingleBit : lda = LdaPlain(snps, LDAbitalign); break;
  case EightByteHaplo : lda = LdaPlainDoubled(snps, LDAbitalign); break;
  case FourByteSingleBit_Tr :
    lda = LdaPlainDoubled(individuals, LDAbitalign);
    break;
  case EightByteHaplo_Tr :  lda = LdaPlain(individuals, LDAbitalign); break;
  default :
    PRINTF("unknown '%s' %d\n", CODING_NAMES[coding],coding);
    BUG;
  }
  //  printf("%s %s %ld %d = %d\n", LDAbitalign ? "Getlda" : "Getlda(bitalign)", CODING_NAMES[coding], snps, individuals, lda);
  if (lda == 0) BUG;
  return lda;
}


Long SnpLDA(Long snps, coding_type coding, Long LDAbitalign) { // LDA in unit_t
  return GetLDA(snps, NA_LONG, coding, LDAbitalign);
}

  
Long GetLDAbitalign(coding_type coding) {
  // assuming that snps=1 behaves adequately
  return GetLDA(1L, 1L, coding, 0) * BitsPerUnit;  // OK
}

Long GetNatLDA(Long snps, Long individuals, coding_type coding) { // LDA in unit_t
  return GetLDA(snps, individuals, coding, GetLDAbitalign(coding));
}


Long SnpNatLDA(Long snps, coding_type coding) { // LDA in unit_t
  return GetLDA(snps, NA_LONG, coding, GetLDAbitalign(coding));
}


Long totalMem(Long snps, Long individuals,
	      Long lda, // in units  unit_t  except OrigPlink
	      coding_type coding,  bool OneHaploOnly, bool relaxed) {
  // in units unit_t  !!!
  //  printf("toatlMem : lda=%ld %d\n", lda, has_lda_in_Byte(coding));
  assert(lda > 0);
  if (has_lda_in_Byte(coding)) lda = DIV_GEQ(lda, sizeof(unit_t));
  assert(FirstUserHaploCoding == 20 && LastUserHaploCoding == 24);
  if (OneHaploOnly) {
    switch(coding) {
    case TwoBitHaplo : coding = OneBitHaplo; break;
    case EightByteHaplo : coding = FourByteHaplo; break;
    case FourByteSingleBit_Tr : coding = FourByteGeno; break;
    default : {} // none
    }
  }
  if (transposedHaplo(coding))
    return lda * snps * (1 + (coding == EightByteHaplo_Tr) * (!OneHaploOnly));
  // printf("totalmem %s %d %d %d %d \n", CODING_NAMES[coding], UnCompressedDoubleCols(coding, false), doubledCols(coding, false),!isHaplo(coding),!OneHaploOnly);
  if (transposedGeno(coding)) individuals = snps;
  if (coding == FiveCodes || coding == FiveCodesTransposed)
    return lda * (1L + (individuals - 1L) / 20L) << 2;
  if (doubledCols(coding, relaxed) && (!isHaplo(coding) || !OneHaploOnly))
    individuals <<= 1;
  return lda * individuals;
}

Long totalMem(Long snps, Long individuals, Long lda,
	      coding_type coding,  bool OneHaploOnly) {
  // in units unit_t  !!!
  return totalMem(snps, individuals, lda, coding, OneHaploOnly, false);
}


Long calculateAlignedMem(Long memInUnits, coding_type coding) {
  if (is32BitAligned(coding)) {
    // R is 32bit-aligned. So no additional alignment considerations
    // necessary. Do not add safe byte either, since vector/matrix
    // returned in readable form to user
    return memInUnits;
  }

  Long UPB = MAX_LDABITALIGN / (sizeof(unit_t) * BitsPerByte); // MAX_LDABITALIGN always save
  Long mem = ROUND_GEQ(memInUnits, UPB);
   
 // + 1 UPB : starting point for alignedmem
  // + 1 UBP : safely at the end
  mem += UPB << 1;
  
  return mem;
}

Long RoundUpSnps(Long snps, coding_type coding, int variant) {
  Long CpB = GetCodesPerBlock(coding, variant);
  return   (1 + (snps - 1) / CpB) * CpB;
}


Long  CECV = 0;
inline static void swap(Long *A, Long *B) { Long tmp = *A; *A = *B; *B=tmp; }
SEXP CompleteCodeVector(Long snps, Long individuals, coding_type coding,
			int variant, Long LDAbitalign, bool OneHaploOnly,
			option_type *global, utilsoption_type *utils,
			void* pointer, // to some existing memeory
			void* nextpointer, // to some existing memeory
			SEXP Code) {
  //  printf(" coding = %.50s %d %ld\n", CODING_NAMES[coding], variant, LDAbitalign);
 
  basic_options *opt = &(utils->basic);
  SEXP Info = R_NilValue;
  Long *info;
  Long
    lda = 0,
    memInUnits =0,
    alignedMem = 0;
  bool CodeGiven = Code != R_NilValue,
    IsHaplo = false;
  if (CodeGiven && pointer != NULL) BUG;
  if (pointer == NULL && nextpointer != NULL) BUG;
  coding_type newcoding = UnknownSNPcoding;
  if (CodeGiven) {
#if defined compatibility_to_R_h
    PROTECT(Code);
#endif  
    Info = getAttribPointer(Code, Information);
  }
  bool info_new = Info == R_NilValue;
  if (info_new) {
    PROTECT(Info = allocVectorExtra(LONGSXP, (INFO_LAST + 1) ));
    info = LONG(Info);
    for (int J=0; J<=INFO_LAST; J++) info[J] = NA_LONG;
    info[MISSINGS] = 0;
  } else {
#if defined compatibility_to_R_h
    PROTECT(Info);
#endif  
    info = LONG(Info);
  }
  
  if (info[BIGENDIAN] != opt->bigendian) { 
    if (info[BIGENDIAN] == NA_LONG)// happens also in case of filename  
      info[BIGENDIAN] = opt->bigendian;
    else ERR2("big/little endian do not match (info=%ld; param=%d)",
	      info[BIGENDIAN], opt->bigendian);
  }

  if (info[IMPLEMENTATION] != CURRENT_IMPLEMENTATION) {
    if (info[IMPLEMENTATION] == NA_LONG)
      info[IMPLEMENTATION] = CURRENT_IMPLEMENTATION;
    else ERR0("SNP matrix in obsolete format cannot be read");
  }
  
  if (info[SNPS] != snps) {
    if (info[SNPS] == NA_LONG) info[SNPS] = snps;
    else ERR0("number of snps has changed");
  }
  
  if (info[INDIVIDUALS] != individuals) {
    if (info[INDIVIDUALS] == NA_LONG) info[INDIVIDUALS] = individuals;
    else ERR0("number of individuals has changed");
  }
  
  if (info[ZAEHLER] == NA_LONG) info[ZAEHLER] = ++CECV;
  
  if (LDAbitalign == 0) LDAbitalign = NA_LONG;
  if ((int) LDAbitalign != NA_LONG && coding == UnknownSNPcoding) BUG;
  if (coding == UnknownSNPcoding) {
    if (!CodeGiven) BUG;
    goto ende;
  }

  assert(coding != AutoCoding);
  IsHaplo = isHaplo(coding);
  if (info[ORIGINALLY_HAPLO] == NA_LONG) info[ORIGINALLY_HAPLO] = IsHaplo;
  //if (OneHaploOnly == NA_IN TEGER && !info[ORIGINALLY_HAPLO])
  // Achtung, info[] koennte auch NA sein -- passt trotzdem
  //  OneHaploOnly = false;
  newcoding = coding;
  if (IsHaplo && OneHaploOnly) newcoding = OneSet2G(coding);
   
  //  printf("CODING %s [%d; %d] %d Bitalgn=(%d %d) NA=%ld\n", CODING_NAMES[coding], coding, info[CODING],  IsHaplo, info[LDABITALIGN], LDAbitalign, NA_LONG);
  
  if (LDAbitalign == NA_LONG) LDAbitalign = GetLDAbitalign(coding);
  assert(LDAbitalign != 0);
  if ((LDAbitalign > MAX_LDA_BITALIGN) &&
      LDAbitalign != NA_LONG) ERR0("alignment out of range");
  if (info[LDABITALIGN] != LDAbitalign) {
    if (info[LDABITALIGN] == NA_LONG) info[LDABITALIGN] = LDAbitalign;
    else ERR0("'LDAbitalign' has changed");
  }

  if (opt->Cprintlevel > 7) 
    PRINTF(" coding = %ld x %ld : %.50s v=%d lda_align=%ld\n",
	   snps, individuals, CODING_NAMES[coding], variant, LDAbitalign);
 
  if (snps == NA_LONG || (int) individuals == NA_LONG ||
      LDAbitalign == NA_LONG // || OneHaploOnly == NA_INTE GER
      ) {
    if (!CodeGiven) BUG;
    goto ende;
  }
 
  lda = GetLDA(snps, individuals, coding, LDAbitalign);
  //  printf("***LDA=%lu %lu\n", lda, (Long)  Code);
  
  if (info[LDA] != lda) {
    //    printf("%d %d; %d %lu\n", info[CODING], coding, info[LDA], lda);
    if (info[LDA] == NA_LONG) {
      info[LDA] = lda;
    }
    else if (info[CODING] != coding) {
      ERR2("'coding' has changed incompatibly from '%.50s' to '%.50s'.",
	   CODING_NAMES[info[CODING]], CODING_NAMES[coding]);
    } else ERR2("'LDA' has changed from %ld to %lu",  info[LDA], lda);
  }
    
  memInUnits = totalMem(snps, individuals, info[LDA], coding, OneHaploOnly);
  if (info[MEMinUNITS] == NA_LONG) info[MEMinUNITS] = memInUnits;
  else if (info[MEMinUNITS] != memInUnits)
    ERR0("memory requirement changed");
  alignedMem = calculateAlignedMem(memInUnits, coding);
  
  if (info[ALIGNEDUNITS] == NA_LONG) {info[ALIGNEDUNITS] = alignedMem;}
  else if (info[ALIGNEDUNITS]!=alignedMem) ERR0("alignment changed");

  {
    Long tmpIndiv = individuals,
      tmpSnps = snps;
    switch(coding) {
    case FourByteSingleBit_Tr :
      swap(&tmpSnps, &tmpIndiv);
      FALLTHROUGH_OK;
    case EightByteHaplo:  
      tmpSnps *= (2- OneHaploOnly);
      if (tmpSnps * tmpIndiv != memInUnits ||
	  alignedMem != memInUnits) BUG;
      break;
    case EightByteHaplo_Tr : 
      swap(&tmpSnps, &tmpIndiv);
      FALLTHROUGH_OK;
    case FourByteSingleBit : case FourByteHaplo :
      tmpIndiv *= (2 - OneHaploOnly);
      FALLTHROUGH_OK;
    case FourByteGeno :
      // printf("%ld * %ld = %ld %ld %ld %ld \n", tmpSnps , tmpIndiv, tmpSnps * tmpIndiv, memInUnits, alignedMem, info[LDA]);
      if (tmpSnps * tmpIndiv != memInUnits || alignedMem != memInUnits) BUG;
      break;
    case OneByteHaplo  :
      FALLTHROUGH_OK;
    default: 
      if (snps > 1000 && snps * individuals <= memInUnits) {
	PRINTF("%s %ld * %ld <= %ld\n",
	      CODING_NAMES[coding], snps, individuals, memInUnits);
	BUG;
      }
      if (Code == R_NilValue) {
	if (pointer == NULL) PROTECT(Code = allocVector(INTSXP, alignedMem));
	else PROTECT(Code = allocVectorPointer(INTSXP, alignedMem, pointer));
      }
    }

  
    //  printf("xxxx=%d \n", Code==R_NilValue);
    if (Code==R_NilValue) { // !!!ACHTUNG: ja nicht CodeGiven abfragen!!!
      assert(tmpSnps * tmpIndiv == alignedMem);
      if (pointer == NULL) PROTECT(Code = allocMatrix(INTSXP,tmpSnps,tmpIndiv));
      else PROTECT(Code = allocMatrixPointer(INTSXP,tmpSnps,tmpIndiv, pointer));
      //   printf("matrix: alignedMem = %lu snps=%d->%lu indiv=%d->%lu; %lu\n", alignedMem, snps, tmpSnps, individuals, tmpIndiv, (Long) INTEGER(Code));
    }
  
   
    if (!CodeGiven) MEMSET(INTEGER(Code), 0, BytesPerUnit * alignedMem); // OK
  
    //  printf("xxxx\n");  return R_NilValue;
    if (info[CODING] != newcoding) {
      if (info[CODING] == NA_LONG) info[CODING] = newcoding;
      else ERR0("internal coding has changed");
    }
    //  printf("Dxxxx %s variant = %d\n", CODING_NAMES[coding], variant);
    //     printf("xxxx\n");  return R_NilValue;
    info[VARIANT] = check_variant(newcoding, variant, opt->efficient);
    if (TYPEOF(Code) != STRSXP) {
      info[ALIGNADDR] = (Long) algn_generalL(INTEGER(Code), LDAbitalign);
      info[RECENTALIGNADDR] = (Long) NULL;
      info[ADDR] = (Long) INTEGER(Code);
    }
  }
 
 ende:

  //  printf("ende: %ld\n", (Long) Code);
  if (Code == R_NilValue) BUG;
  //  printf("info_new %d\n", info_new);  
  if (info_new) setAttrib(Code, Information, Info); // install("information")
  //printf("info_new end %d\n", info_new);  
  
  //  printf("xxxx=%d %d\n", Code==R_NilValue, GetInfo(Code)[100000]);

   if (info[ORIGINALLY_HAPLO] != NA_LONG // && OneHaploOnly != NA_INTE GER
      &&  GET_CLASS(Code) == R_NilValue) {
    //    printf("class: %d %d ", info[ORIGINALLY_HAPLO] , OneHaploOnly);
    SEXP Klasse;  
    PROTECT(Klasse = allocVector(STRSXP, 1));
    SET_STRING_ELT(Klasse, 0,
		   mkChar(info[ORIGINALLY_HAPLO] && !OneHaploOnly
			  ? HAPLOMATRIX : GENOMICMATRIX));
    SET_CLASS(Code, Klasse);
    UNPROTECT(1);
  }
  
  if (hasTransposed(coding) && global->tuning.addtransposed) {
    SEXP NewNext = PROTECT(CompleteCodeVector(snps, individuals,
					      switch_transposed(coding),
					      variant,
					      info[LDABITALIGN],
					      OneHaploOnly,
					      global, utils,
					      nextpointer, NULL,
					      getAttribPointer(Code, Next)));
    setAttrib(Code, Next, NewNext);
    UNPROTECT(1);
  }
      
  //  printf("done Completing\n");
  UNPROTECT(2);
  return Code;
}


/*
SEXP CompleteCodeVector(Long snps, Long individuals, coding_type coding,
			int variant, Long LDAbitalign, bool OneHaploOnly,
			option_type *global, utilsoption_type *utils,
			SEXP Code) {
  CompleteCodeVector(snps, individuals, coding,
		     variant, LDAbitalign, OneHaploOnly,
		     global,utils, NULL, NULL, Code);

}
*/


SEXP CreateEmptyCodeVector(Long snps, Long individuals, coding_type coding,
			   int variant, Long LDAbitalign, bool OneHaploOnly,
			   option_type *global,  utilsoption_type *utils) {
  return CompleteCodeVector(snps, individuals, coding,
			    variant, LDAbitalign, OneHaploOnly,
			    global,utils, NULL, NULL, R_NilValue);
}


void start_info(SEXP Code) {
  Long *info = GetInfoUnchecked(Code);
  coding_type coding =  (coding_type) info[CODING];
  Long
    mem =info[MEMinUNITS],
    individuals = info[INDIVIDUALS],
    snps = info[SNPS],
    lda = GetLDA(snps, individuals, coding, info[LDABITALIGN]); 
  bool inMB = mem > 5000000;
  PRINTF("Data:  %ld SNPs for each of %ld individuals\nStorage mode: (approx) %ld block(s) of %ld codes\nSize of M: %ld %sB.", 
	 info[SNPS], individuals,  mem / lda,
	 GetCodesPerBlock(coding, (int) info[VARIANT]),
	 1 + (mem - 1) / (inMB ? 1048576 : 1024),
	 inMB ? "M" : "k");	       
}



SEXP createSNPmatrix(Long snps, Long individuals,
		     coding_type coding, int variant, Long LDAbitalign,
		     bool OneHaploOnly, SEXP V,
		     option_type *global,  utilsoption_type *utils) {  
  basic_options *opt = &(utils->basic);
  SEXP Code;
  if (coding == FileDot || coding == DotFile) {  
    // NOTE! no attr(Filename) created, no attr(Info) created, as
    // "SNP matrix" is only a virtual one
    Long len = coding == FileDot ? snps : individuals;
    if (V != R_NilValue) {
      if ((Long) LENGTH(V) != (coding != FileDot ? snps :  individuals))
	ERR1("vector must have length equal to number of %.20s.",
	     coding == FileDot ? "individuals" : "snps");
    }
    if (opt->Cprintlevel > 6) {
      PRINTF("Data: %ld SNPs for each of %ld individuals read from file\n",
	     snps, individuals);
    }
    Code = PROTECT(allocVector(REALSXP, len));
    MEMSET(REALX(Code), 0, len * sizeof(double));
  } else {
    if (variant == VARIANT_GPU) check_7_5();      
    Code = PROTECT(CreateEmptyCodeVector(snps, individuals, coding, variant,
					 LDAbitalign, OneHaploOnly,
					 global,utils));
    if (opt->Cprintlevel > 6) start_info(Code);
  }
  
  UNPROTECT(1);
  return Code;
}


SEXP createSNPmatrix(Long snps, Long individuals,
		     coding_type coding, int variant,
		     Long LDAbitalign,
		     option_type *global,  utilsoption_type *utils) {
  return createSNPmatrix(snps, individuals, coding, variant,
			 LDAbitalign, false, R_NilValue, global, utils);
}


void showInfo(Long *info) {
  PRINTF("information variant=%ld\n", info[IMPLEMENTATION]);
  PRINTF("snps=%d\n", SNPS);
  PRINTF("indiv=%d\n", INDIVIDUALS);
  PRINTF("Addr=%d\n", ADDR);
  PRINTF("Align=%d\n", ALIGNADDR);
  PRINTF("coding=%d\n", CODING);
  PRINTF("implementation=%d\n", VARIANT);
  PRINTF("LDAbitalign=%d\n", LDABITALIGN);
  PRINTF("snpXind=%d\n", FILE_SNPxIND_READ_BYROW);
  PRINTF("header=%d\n", FILE_HEADER);
  PRINTF("doubleindiv=%d\n", FILE_DOUBLEINDIV);
  PRINTF("leading=%d\n", FILE_LEADINGCOL);
  PRINTF("mem=%d\n", MEMinUNITS);
  PRINTF("alignedmem=%d\n", ALIGNEDUNITS);
  PRINTF("LDA=%d\n", LDA);
}

void showInfoShort(Long *info) {
  PRINTF("V=%ld, ", info[IMPLEMENTATION]);
  PRINTF("S/I=%d %d, ", SNPS, INDIVIDUALS);
  PRINTF("Ad=%d, ", ADDR);
  PRINTF("Al=%d, ", ALIGNADDR);
  PRINTF("M=%d, ", CODING);
  PRINTF("I=%d\n", VARIANT);
  PRINTF("LDAbitalign=%d, ", LDABITALIGN);
  PRINTF("LDA=%d \n", LDA);
}

#if defined compatibility_to_C_h  
extern utilsoption_type OPTIONS;
#endif
unit_t *Align(SEXP SxI,  basic_options *opt) {
  assert(SxI != R_NilValue);
  int *address = INTEGER(SxI);
  Long *info = GetInfo(SxI);
  coding_type coding = (coding_type) info[CODING];
  if (coding == OrigPlink || coding == OrigPlinkTransposed)
    return (unit_t *) address;
      
   option_type *global;
  utilsoption_type *utils;
#if defined compatibility_to_R_h
  KEY_type *KT = KEYT_M();
  global = &(KT->global);
  utils = &(KT->global_utils);
#else
  global = &OPTIONS_MIRACULIX;
  utils = &OPTIONS;
#endif
  if (opt == NULL) opt = &(utils->basic);

  int cores = GreaterZero(opt->cores);
  Long bitalign = info[LDABITALIGN];
  ASSERT_STORED_ENDIAN;
 
   int
   *algnaddress = (int*) algn_generalL(address, bitalign),
    *infoaddress = (int*) info[ADDR];
    
  bool aligned= algnaddress - address ==
    (int*) algn_generalL(infoaddress, bitalign) - infoaddress;
      
  if (global->messages.warn_address && infoaddress != address) {
    int *recentaddress = (int*) info[RECENTALIGNADDR];
    if (recentaddress != address) {
      PRINTF("Address has changed in a coded object%s (*%u->*%u). %s\n",
	     aligned ?  ", but same modulus" : " by 'gc' or disk saving",
	     (Uint) (uintptr_t) infoaddress, (Uint) (uintptr_t) address,
	     recentaddress == NULL
	     ? "[This messages can be suppressed by 'RFoptions(warn_address=FALSE)']"// OK
	     : "");
      info[RECENTALIGNADDR] = (Long) address;
    }
  }
  
  if (aligned){
    assert((uintptr_t) algnaddress % (bitalign / BitsPerByte) == 0);
    return  (unit_t*) algnaddress;
  }

  Long bytes = info[MEMinUNITS] * BytesPerUnit; // OK
#ifdef DO_PARALLEL
  // prevent two processes to move data at the same time
  if (cores > 1) {
    int mypid, sleep = 1;
    Ext_pid(&mypid);
    if (!info[BLOCKEDINFO]) {
      info[BLOCKEDINFO] = mypid;
      Ext_sleepMicro(&sleep); // wait. Other processes may access simulatenously
    }
    if ((Long) mypid != info[BLOCKEDINFO]) { // other process has been
      // doing the job, maybe simultaneously
      sleep = (int) ((bytes / 1000 + mypid) % MAXINT); // some randomness in sleeping
      for (int i=0; i<10; i++) {
	if (!info[BLOCKEDINFO])
	  return Align(SxI, opt);
 	Ext_sleepMicro(&sleep);
      }
      if (opt->warn_parallel) 
	WARN0("Simultaneous write access: there is some risk of loosing data. You can suppress this message by 'RFoptions(warn_parallel=FALSE)', or avoid any risk by 'RFoptions(cores=1)'."); // ok
      info[BLOCKEDINFO] = 0; 
      return Align(SxI, opt);
    }
  }
#endif
  
  int *OldInNew =
    address + ((Long) info[ALIGNADDR] - (Long) infoaddress);
  MEMMOVE(algnaddress, OldInNew, bytes);
  info[ADDR] = (Long) address;
  info[ALIGNADDR] = (Long) algnaddress;
#ifdef DO_PARALLEL
  info[BLOCKEDINFO] = 0;
#endif

  return (unit_t*) algnaddress;
}


typedef void (*get_uint_t )(unit_t *, Long snps, Long indiv, Long lda,
			    int cores, unit_t *ans, Long ldAns);


void matrix_get_4Byte(unit_t *Code, Long snps, Long indiv, Long lda,
		      coding_type coding,
		      int cores, unit_t *ans, Long ldAns){
  MEMSET(ans, 0, totalMem(snps, indiv, lda, FourByteGeno,
			  false) * BytesPerUnit);// OK
  get_uint_t get = NULL;
  switch (coding) {
  case OneBitGeno : get = get_matrix1_4Byte; break;
  case TwoBitGeno : get = get_matrix2_4Byte ; break;
  case ThreeBit : get = get_matrix3_4Byte; break;
  case Plink : get = get_matrixPlink_4Byte; break;
  case FourByteGeno:  get = get_matrixPlain_4Byte; break;
  case OneByteGeno:  get = get_matrixOneByte_4Byte; break;
  default :
    BUG;
  }
  get(Code, snps, indiv, lda, cores, ans, ldAns);
}


typedef void (*multi_t)(unit_t *x, unit_t *y, Long indiv, Long blocks, Long lda,
			double *ans);
typedef Ulong (*scalar_t)(unit_t *x, unit_t *y, Long blocks, Long Delta);
   
#define ASSERT_NO_TILING(coding, variant)			\
  assert(exists_crossprod(coding) &&				\
	 exists_variant(coding, variant, false, false) &&	\
	 exists_tiling(coding, variant, false));
#define ASSERT_TILING(coding, variant)				\
  assert(exists_crossprod(coding) &&				\
	 exists_variant(coding, variant, false, false) &&	\
	 exists_tiling(coding, variant, true));


void crossprod(unit_t * Code, Long snpsOrig, Long individuals, Long lda,
	       coding_type coding, int variant,
	       Long logSnpGroupSize, bool smoothSnpGroupSize,
	       Long indivGroupSize, basic_options *opt,
	       double *A) {

  Long snps = RoundUpSnps(snpsOrig, coding, variant);
 
  Long snpGroupSize = 1U << logSnpGroupSize;
  int cores =  GreaterZero(opt->cores);
  if (snpGroupSize <= 1) snpGroupSize = 1 - snpGroupSize; // switch 0 and 1
  indivGroupSize = indivGroupSize > 0 ? indivGroupSize :
    cores == 1 ? individuals : 1;

  if (opt->Cprintlevel > 6) 
    PRINTF("cross: %ld/%ld x %ld : %.50s,%d lda=%ld; tiling=%ld x %ld.\n",
	   snpsOrig, snps, individuals,
	   CODING_NAMES[coding], variant, lda,
	   snpGroupSize, indivGroupSize);
  
  assert(individuals > 0);
  double nSq = (double) individuals * (double) individuals;      
  if ((double) snpsOrig * nSq * 4.0 > 4.5035996e+15) // 2^{manissen-bits = 52}
    ERR0("matrix too large to calculate the relationship matrix -- pls contact maintainer of 'miraculix'");
      
  
  Long
    blocks = GetBlocks(snps, coding, variant),
    BpB = GetBytesPerBlock(coding, variant),
    UpB = BpB / BytesPerUnit, // OK
    minigroupsize = 1,
    BlockSizeFactor = 0,
    DeltaPart2 = 0;

  //  printf("crossprod %s %d %d \n", CODING_NAMES[coding], variant, snps);

  scalar_t scalar = NULL;
  multi_t multi = NULL;
  switch (coding) {
  case OneBitGeno : 
    ASSERT_TILING(coding, variant);   
    DeltaPart2 = lda * individuals / UpB;
    BlockSizeFactor = loop_1bit;
    switch(variant) {
    case 512 + VARIANT_B :
      multi  = multi_1v512B;
      scalar = scalar_1v512B;
      break;
    case 512 + VARIANT_A : 
      multi  = multi_1v512A;
      scalar = scalar_1v512A;
      break;
    case 512 :
      multi  = multi_1v512;
      scalar = scalar_1v512;
      break;
    case 256 + VARIANT_B :
      multi  = multi_1v256B;
      scalar = scalar_1v256B;
      break;
    case 256 + VARIANT_A : 
      multi  = multi_1v256A;
      scalar = scalar_1v256A;
      break;
    case 256 :
      multi  = multi_1v256;
      scalar = scalar_1v256;
      break;
    case 128 + VARIANT_B :
      multi  = multi_1v128B;
      scalar = scalar_1v128B;
      break;
    case 128 + VARIANT_A : 
      multi  = multi_1v128A;
      scalar = scalar_1v128A;
      break;
    case 128 :
      multi  = multi_1v128;
      scalar = scalar_1v128;
      break;
    case 64 :
      multi  = multi_1v64;
      scalar = scalar_1v64;
      break;
    case 32 :
      multi  = multi_1v32;
      scalar = scalar_1v32;
      break;
    default: BUG;
    }
    break;
  case TwoBitGeno :
    BlockSizeFactor = loop_2bit;
    switch(variant) {
    case 512  + VARIANT_A:
      ASSERT_TILING(coding, variant);
      minigroupsize = 2;      
      multi =  multi2x2_2v512;
      scalar = scalar_2v512;
      break;
    case 512 :
      ASSERT_TILING(coding, variant);
      multi =  multi_2v512;
      scalar = scalar_2v512;
      break;
    case 256  + VARIANT_A:
      ASSERT_TILING(coding, variant);
      minigroupsize = 2;      
      multi  =  multi2x2_2v256;
      scalar = scalar_2v256;
      break;
    case 256 :
      ASSERT_TILING(coding, variant);
      multi =  multi_2v256;
      scalar = scalar_2v256;
      break;
    case 128  + VARIANT_A:
      ASSERT_TILING(coding, variant);
      minigroupsize = 2;      
      multi =  multi2x2_2v128;
      scalar = scalar_2v128;
      break;
    case 128 :
      ASSERT_TILING(coding, variant);
      multi  =  multi_2v128;
      scalar = scalar_2v128;
      break;
    case 64 :
      ASSERT_NO_TILING(coding, variant);
      crossprod2v64(Code, snps, individuals,lda, cores, A); return;
    case 32 :
      ASSERT_NO_TILING(coding, variant);
      crossprod2v32(Code, snps, individuals,lda, cores, A); return;
    case VARIANT_GPU :
      stopIfNotUInt(snps | individuals | lda);
      crossprod_mmagpu((Uint*) Code, snps,  individuals, (Uint) lda, cores, A);  return;
    default: BUG;
    }
    break;
  case ThreeBit:
    ASSERT_NO_TILING(coding, variant);
    switch(variant) {
    case 64: crossprod3(Code, snps, individuals, lda, cores, A);return;
    default: BUG;
    }
    return;
  case FourByteGeno :
    ASSERT_NO_TILING(coding, variant);
    switch(variant) {
    case 32 : crossprod_Plain (Code, snps, individuals, lda, cores,  A);
      return;
    case VARIANT_GPU :
#ifdef USEGPU
      crossprod_PlainGPU(Code, snps, individuals, lda, cores, A);
      return;
#else
      BUG;
#endif    
    default : BUG;
    }
  case OneByteGeno :
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic) 
#endif
    for (Long j=0; j<individuals; j++) {
      for (Long k=j; k<individuals; k++) {
	char *J = (char*) (Code + lda * j),
	  *K = (char*) (Code + lda * k);
	Long  sum = 0L;
	for (Long i=0; i<snps; i++) sum += J[i] * K[i];
	A[j + individuals * k] = A[k + individuals * j] = (double) sum;
      }
    }
    return;
    
  default:
    ERR1("crossprod does not exist for coding `%.50s'.\n", CODING_NAMES[coding]);
  }

  if (scalar == NULL) BUG;

  //  printf("minigroup = %d %d\n", minigroupsize, BlockSizeFactor);

#define minBlocksPCore 5

  indivGroupSize = (1 + (indivGroupSize -1) / minigroupsize) * minigroupsize;
  const Long
    snpGroupSizePlus = snpGroupSize + BlockSizeFactor / 4,
    blockSize = (snpGroupSize == 0 || snpGroupSize > blocks
		 ? blocks
		 : BlockSizeFactor == 0 || snpGroupSizePlus < BlockSizeFactor
		 || !smoothSnpGroupSize ? snpGroupSize
		 : (snpGroupSizePlus / BlockSizeFactor) * BlockSizeFactor),
    lastRegBlock = blockSize * ((blocks-1) / blockSize),
    blockGsize[2] = {blockSize, blocks - lastRegBlock},
    rest = individuals % indivGroupSize,
    lastRegIndi = individuals - rest,
    grouprest = rest % minigroupsize,
    indivGsize[3] = {indivGroupSize, rest-grouprest, grouprest};
  

 
  //  printf("blocksize = %d %d %d : %d %d; snpGroupSize=%d blocks=%d\n", blockSize, snpGroupSize <= 1 || snpGroupSize > blocks * CpB,	 BlockSizeFactor == 0 || snpGroupSizePlus < BlockSizeFactor || !smoothSnpGroupSize,	 snpGroupSizePlus, BlockSizeFactor, snpGroupSize, blocks	 );
 


  //  printf("\n\nsnps=%d D=%d lda=%d i=%d BitsPBlock=%d/%d\niGrS=%d (%d, %d, %d) bs=%d (%d %d) lastBl=%d %d\n",	 snps, DeltaPart2, lda, individuals , BpB * 8, UpB,	 indivGroupSize, indivGsize[0], indivGsize[1], indivGsize[2], blockSize,	 blockGsize[0],  blockGsize[1],  lastRegBlock, lastRegIndi );
  // BUG;
  
  MEMSET(A, 0, sizeof(*A) * individuals * individuals);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) collapse(2) schedule(static)
#endif
  for(Long iGroup = 0; iGroup < individuals; iGroup+=indivGroupSize) {
    for(Long jGroup = 0; jGroup<individuals; jGroup+=indivGroupSize) {
      if (jGroup < iGroup) continue;
      const bool last_reg = iGroup >= lastRegIndi;
      const Long lEnd_i = iGroup + indivGsize[last_reg];
      const Long lEnd_j = jGroup + indivGsize[jGroup >= lastRegIndi];      
      for(Long b = 0; b < blocks; b += blockSize) {	  
	const Long bsize = blockGsize[b == lastRegBlock];
	unit_t *Code0 = Code + b * UpB;
	for(Long i = iGroup; i < lEnd_i; i+=minigroupsize) {
	  unit_t *Codei = Code0 + i * lda;
	  for(Long j = MAX(i, jGroup); j < lEnd_j; j+=minigroupsize) {
	    //	      printf("iGr=%d %d; i=%d %d b=%d %d blocks=%d lda=%d addr=%ld %ld\n", iGroup, jGroup, i, j, b, bsize, blocks, lda, Code0 - Address, Codei - Address);
	    
	    assert(j * individuals + i < individuals * individuals);
	    // printf("%ld %ld\n", Codei-Code, Code0 + j * lda-Code);
	    
	    multi(Codei, Code0 + j * lda, individuals, bsize, lda,
		  A + j * individuals + i);
	    //	      if (*A > 2000) { printf("%f %d/%d %d/%d %d\n", *A, indivGroupSize, lastRegIndi, blockSize, blocks, minigroupsize); BUG;} BUG;
	  }
	}
      }
    }
  }
    
  if (indivGsize[2] > 0) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
    for(Long iGroup = 0; iGroup < individuals; iGroup+=indivGroupSize) {
      const bool last_reg = iGroup >= lastRegIndi;
      const Long xEnd_i = last_reg ? individuals : iGroup + indivGroupSize,
	lEnd_j = individuals;
      Long iG = iGroup,
	jGroup = individuals - grouprest;
      for(Long b = 0; b < blocks; b += blockSize) {	  
	const Long bsize = blockGsize[b == lastRegBlock];
	unit_t *Code0 = Code + b * UpB;
	for(Long i = iG; i < xEnd_i; ++i) {
	  unit_t *Codei = Code0 + i * lda;
	  for(Long j = MAX(i, jGroup); j < lEnd_j; ++j) {
	    //	    printf("i=%d j=%d\n", i, j);
	    assert(j * individuals + i < individuals * individuals);
	    A[j * individuals + i] += (double)
	      scalar(Codei, Code0 + j * lda, bsize, DeltaPart2);
	  }
	}
      }
    } // for iGroup
  }
     
  for (Long i=0; i<individuals; i++)
    for (Long j=0; j<i; j++)
      A[j * individuals + i] = A[i * individuals + j];
}




Ulong sumGeno(unit_t *Code, Long snps, Long individuals, Long lda,
	      coding_type coding, int variant, int cores,
	      Ulong *sums // [0..individuals-1]
	      ) {
  // total sum of all genos
  variant = main_variant(variant);
  sumGeno_t sG = NULL;
  
  // printf("sumGeno : %sv%d %ld %s\n", CODING_NAMES[coding], variant, (Long) sums, CODING_NAMES[OneByteGeno] );
  
  switch (coding) {
    
  case OneBitGeno : sG = sumGeno1; break;
  case TwoBitGeno :
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 : sG = sumGeno2v256; break;
    case 128: sG = sumGeno2v128; break;
    default :sG = sumGeno2; 
    }
    break;
  case ThreeBit : sG = sumGeno3; break;
  case FourByteGeno: sG = sumGenoPlain; break;
  case OneByteGeno : sG = sumGenoOneByte; break;
  case Plink : sG = sumGenoPlink; break;
  case FiveCodesTransposed : BUG;
  case TwoBitGenoTransposed : BUG;
  case PlinkTransposed : BUG;    
  case OrigPlink : BUG; 
  case OrigPlinkTransposed : BUG; 
  case FiveCodes : 
    // calculation with arguments as information not possible
    BUG;
    break;

  default : BUG;
  }

  return sG(Code, snps, individuals, lda, cores, sums);
  
}


void do_centering(Long Snps, Long Individuals, centering_type centered,
		  normalizing_type normalized, bool squared,
		  LongDouble *Precision,
		  basic_options VARIABLE_IS_NOT_USED *opt,
		  tuning_options VARIABLE_IS_NOT_USED *tuning,
		  double *A, bool transposed) {
  if ((centered == NoCentering && normalized == NoNormalizing )
      || centered == User) return;
  assert(Individuals > 0);

  
  //bool exact = opt->basic.exactness == True || Precision == NULL;
  
  bool exact =  Precision == NULL;
  // printf("exact not exact\n");
  
  const Long cols = transposed ? Snps : Individuals;
  const Long rows = transposed ? Individuals : Snps;
  double totalsum = -1.0;
  const double cols_D = (double) cols;
  const double rows_D = (double) rows;
  LongDouble *FreqSxI=NULL, *pseudoFreq=NULL;
  if (Precision != NULL) {
    Long info[INFO_LAST];
    info[SNPS] = Snps;
    info[INDIVIDUALS] = Individuals;
    totalsum = (double) Precision[transposed ? PseudoSqTotalSumIndex
				  : SqTotalSumIndex];
    FreqSxI = Precision + (transposed ? PseudoFreqSxIIndex : FreqSxIIndex);
    pseudoFreq = Precision + (transposed ? FreqIndex : PseudoFreqIndex);
  }
  if (transposed) centered = centered == RowMeans ? ColMeans
		    : centered == ColMeans ? RowMeans : centered;
  double factor = 1.0,
    *sums = (double *) MALLOC(sizeof(double) * cols);   

  switch(centered) {
  case RowMeans:
  case NoCentering: {
#define maxlong 9223372036854775000.0 // for safty
    if (exact && 2.0 * 4.0 * cols_D * cols_D * rows_D > maxlong) {
      WARN1("Caught overflow %.50s.", CONTACT);
      if (Precision == NULL) ERR0("information missing in do.centering");
      exact = false;
    }
     
    const double  nSq = exact ? cols_D * cols_D : 1.0;
   
    Long nP1 = cols + 1;
    if (exact) {
      totalsum = 0.0;
      for (Long i=0; i<cols; i++) { // ok
	double dummy = 0.0;      
	double *a = A + i * cols;
	for (Long j=0; j<cols; j++) dummy += a[j];// ok
	sums[i] = cols_D * dummy;
	totalsum += dummy;
      }
    } else {
      for (Long i=0; i<cols; i++)
	sums[i] = 2.0 * (double) FreqSxI[i];
    }

    if (centered==RowMeans) {
      for (Long i=0; i<cols; i++) {
	Long idx = i * nP1;
	for (Long j=i; j<cols; j++, idx++) { 
	  A[idx] = nSq * A[idx] - sums[i] - sums[j] + totalsum;
	}
      }
      factor = nSq;
    }
   
  }
    break;
    
  case ColMeans : {
    if (exact) {
      ERR0("colmeans needs frequencies");
    } else {
      for (Long i=0; i<cols; i++) {
	sums[i] = ROUND(2.0 * rows_D * (double) pseudoFreq[i]);//be very precise
      }
    }
    
    for (Long i=0; i<cols; i++) {
      double *a = A + i * cols;
      for (Long j=0; j<cols; j++) {
	a[j] = (a[j] * rows_D - sums[i] * sums[j]) / rows_D;
      }
    }

  }
    break;
    
  default :
    ERR0("user defined centering to be programmed\n");
  }

  
  switch (normalized) {
  case GlobalNormalizing :
    if (exact) ERR0("normalizing needs frequencies")
      else {
	//^= sigma^2 * ind^2
	double sum_P = cols_D * (double) Precision[TotalSumIndex];
	factor *= (sum_P - 0.5 * totalsum) / (cols_D * cols_D);
      }
    break;
  case AccordingCentering :
    switch(centered) {
    case RowMeans : case ColMeans :
      for (Long i=0; i<cols; i++)
	sums[i] = 1 / SQRT(A[i * (cols + 1)]); // 1/ SQRT(diag elements)
      for (Long i=0; i<cols; i++) {
	double *a = A + i * cols;
	for (Long j=0; j<cols; j++) {
	  a[j] *= sums[i] * sums[j];
	}
      }
      break;
    case User :
    case NoCentering : break;
    default: BUG;
    }
    break;
  case NoNormalizing :
    break;
  default : BUG;
  }
  

  if (factor <= 0.0) ERR0("strange input matrix?");
  if (factor != 1.0) {
    Long nsq = cols * cols;
    for (Long i=0; i<nsq; i++) A[i] /= factor;
  }

  if (squared) {
    Long nsq = cols * cols;
    for (Long i=0; i<nsq; i++) A[i] *= A[i];
  }
   
  FREE(sums);
  return;
}


void crossprodI(unit_t * Code, Long snps, Long individuals,  Long lda,
		coding_type coding, int variant,
		centering_type centered, normalizing_type normalized,
		bool squared,
		LongDouble *Precision,
		basic_options *opt, tuning_options *tuning,
		double *A) {
  //  printf("\n%s #%d cores=%d\n", CODING_NAMES[coding], variant, cores);
  crossprod(Code, snps, individuals, lda, coding, variant,
	    tuning->logSnpGroupSize,
	    tuning->smoothSnpGroupSize, 
	    tuning->indivGroupSize,
	    opt, A);
  // printf("crosprod done %s %f\n", CODING_NAMES[coding], *A);
 
  assert(A[0] >= 0.0);

  do_centering(snps, individuals,
	       centered, normalized, squared,
	       Precision,
	       opt, tuning, A, false);
}



double *crossprod(unit_t * Code, Long snps, Long individuals,  Long lda,
		  coding_type coding, int variant,		  
		  centering_type centered, normalizing_type normalized,
		  bool squared,
		  LongDouble *Precision,
		  basic_options *opt, tuning_options *tuning) {
  // never squared here!!!
  double *A = (double *) MALLOC(individuals * individuals * sizeof(double));

  crossprodI(Code, snps, individuals, lda, coding, variant,
	     centered, normalized, squared, Precision,
	     opt, tuning, A);
 
  return A;
}
 


void t_crossprodI(unit_t * Code, Long snps, Long individuals,  Long lda,
		  coding_type coding,
		  int VARIABLE_IS_NOT_USED variant,
		  centering_type centered, // with respect to the
		  // original orientation of the matrix
		  normalizing_type normalized,
		  bool squared,
		  LongDouble *Precision,
		  basic_options *opt, tuning_options *tuning,
		  double *Ans) {
  Long ldAns = snps;
  int cores =  GreaterZero(opt->cores);
  switch(coding) {
  case OneByteGeno :
    if (sizeof(double) != sizeof(Long)) BUG;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(dynamic) 
#endif
    for (Long k=0; k<snps; k++) {
      for (Long i=0; i<individuals; i++) {
	_1Byte* C = (_1Byte*) (Code + i * lda);
 	Long z = C[k];
	Long *a = ((Long*) Ans) + k * ldAns;
	for (Long j=k; j<snps; j++) {
	  a[j] += (Long) C[j] * z;
	}
      }
    }
    for (Long k=0; k<snps; k++) {
      for (Long j=k+1; j<snps; j++)
	Ans[j * ldAns + k] = Ans[k * ldAns + j] =
	  (double) (((Long*) Ans)[k * ldAns + j]);
    }
    break;
  default : BUG
      }
  
  assert(Ans[0] >= 0.0);

  do_centering(snps, individuals,
	       centered, normalized, squared, Precision,
	       opt, tuning, Ans, true);
}
 
     
void allele_sum_I(SEXP SxI,
		  bool pseudo_freq,
		  option_type *global, utilsoption_type *utils,
		  Ulong *ans) {
  basic_options *opt = &(utils->basic);
  Long *info = GetInfo(SxI);		   
  ASSERT_LITTLE_ENDIAN;  
  coding_type coding = (coding_type) info[CODING];
  bool gV = (!pseudo_freq) xor transposedGeno(coding);
  if (!gV) ERR2("pseudo=%d for '%s' not programmed yet.\n",
		pseudo_freq, CODING_NAMES[coding]);
  Long
    lda = info[LDA],
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    rows = pseudo_freq ? individuals : snps,
    cols = pseudo_freq ? snps : individuals;
  unit_t *code = Align(SxI, opt);
  int variant = main_variant(check_variant(coding, (int) info[VARIANT],
					   opt->efficient));
  bool EightByte = coding == EightByteHaplo || coding == EightByteHaplo_Tr;
  Long ldabitalign = info[LDABITALIGN];
  Ulong *sum = ans;
  if (EightByte) {
    sum = (Ulong*) CALLOC(rows << 1, sizeof(*sum));
  }

  switch (coding) {     
  case TwoBitGeno : {
    assert(exists_allelefreq(coding));
 
    switch(variant) {
    case 512:
    case VARIANT_GPU:
    case 256 : if (opt->bigendian)
	allele_sum_2v256(code, snps, individuals, lda, ldabitalign, opt, sum);
      FALLTHROUGH_OK;
    case 128 : if (opt->bigendian)
	allele_sum_2v128(code, snps, individuals, lda, ldabitalign, opt, sum);
      FALLTHROUGH_OK;
    default: allele_sum2(code, snps, individuals, lda, opt, sum);
    }
  }
    break;
    
  
  case ThreeBit :
    assert(exists_allelefreq(coding));
    allele_sum3(code, snps, individuals, lda, opt, sum);
    break;
  case FourByteHaplo:
  case FourByteSingleBit:
    cols *= 2;
    FALLTHROUGH_OK;
  case EightByteHaplo:
  case FourByteGeno:  // ToDo: SIMD version
    assert(exists_allelefreq(coding));
    for (Long i=0; i<cols; i++) {
      unit_t *c = code + lda * i;
      for (Long j=0; j<rows; j++) sum[j] += c[j];
    }
    break;
    
  case EightByteHaplo_Tr:
    cols *= 2;
    FALLTHROUGH_OK;
  case FourByteSingleBit_Tr:
    FALLTHROUGH_OK;
  case FourByte_Tr: // ToDo: SIMD version
    assert(exists_allelefreq(coding));
    for (Long i=0; i<cols; i++) {
      unit_t *c = code + lda * i;
      Ulong S = 0;
      for (Long j=0; j<rows; j++) S += c[j];
      sum[i] = S;
    }
    break;
   
  case OneByteHaplo:
    cols *= 2;
    FALLTHROUGH_OK;
    
  case OneByteGeno: 
    assert(exists_allelefreq(coding));
    allele_sum_OneByte(code, snps, individuals, lda, ldabitalign, opt, sum);
    break;
  
  
  default :
    if (!exists_allelefreq(coding)) BUG;
    Ulong *E = (Ulong *) MALLOC(sizeof(Ulong) * cols);
    for (Long i=0; i<cols; E[i++] = 1.0);
    //   printf("before gv vg Long\n");
    gV_vG_Ulong(SxI,
		E, 1, cols,
		gV, false, 
		global, utils,
		sum, rows);
    FREE(E);
  }

  if (EightByte) {
    for (Long i=0; i<rows; i++) ans[i] = sum[2 * i] + sum[2*i + 1];
    FREE(sum);
  }
}



#define divide_by_nonmissings_fctn(DOUBLE)				\
  void divide_by_nonmissings_##DOUBLE(SEXP SxI,				\
				      /* !pseudo_freq extra needed see below*/ \
					bool pseudo_freq,		\
					Long *mssngs, Long n_miss,	\
				      Ulong *sums,			\
				      DOUBLE *Sums, DOUBLE *ans) {	\
    if (false) PRINTF("divide\n");					\
    Long *info = GetInfo(SxI);						\
    Long lenFreq = info[pseudo_freq ? INDIVIDUALS : SNPS],		\
      summands = info[pseudo_freq ? SNPS : INDIVIDUALS];		\
    for (Long i=0; i<lenFreq; i++) Sums[i] = (DOUBLE) sums[i];		\
    if (mssngs == NULL) {						\
      if (false) {PRINTF("missings = null; %ld\n", n_miss);}			\
      assert(n_miss == 0 || n_miss == NA_LONG);			\
      DOUBLE factor = (DOUBLE) (summands << 1);				\
      for (Long i=0; i<lenFreq; i++) ans[i] = Sums[i] / factor;		\
    } else {								\
      Long missings = n_miss << 1;					\
      if (false) {printf("missings non null; %ld\n", n_miss);printINFO(SxI); } \
      assert(missings > 0);						\
      Long *nonmiss =  (Long*) MALLOC(lenFreq * sizeof(Long));		\
      for (Long i=0; i<lenFreq; nonmiss[i++] = summands);		\
      /* if pseudo missings are summed colwise: */			\
      for (Long i=(Long)pseudo_freq; i<missings; i+=2) nonmiss[mssngs[i]]--; \
      for (Long i=0; i<lenFreq; i++) {					\
 	ans[i] = Sums[i] / (DOUBLE) (2L * nonmiss[i]);			\
      }									\
      FREE(nonmiss);							\
    }									\
    if (false) PRINTF("dividing done \n");				\
  }

#if defined compatibility_to_R_h
divide_by_nonmissings_fctn(double)
#endif
divide_by_nonmissings_fctn(LongDouble)





#define allele_freq_fctn(DOUBLE)					\
  void allele_freq_##DOUBLE(SEXP SxI, Long *mssngs,  Long missings,	\
			    bool pseudo_freq,				\
			    option_type *global, utilsoption_type *utils, \
			    DOUBLE *Sums, /* might be NULL */		\
			    DOUBLE *ans) {				\
    if (false) PRINTF("start allele_freq\n");				\
    Long *info = GetInfo(SxI);						\
    assert(info != NULL);						\
    Long lenFreq = info[pseudo_freq ? INDIVIDUALS : SNPS];		\
    Ulong *sum =  (Ulong*) CALLOC(lenFreq, sizeof(Ulong));		\
    DOUBLE *S = Sums;							\
    if (false) PRINTF("before allele sum \n");				\
    allele_sum_I(SxI, pseudo_freq, global, utils, sum);			\
    if (false) PRINTF("dividing \n");					\
    if (S == NULL) S = (DOUBLE *) MALLOC(lenFreq * sizeof(DOUBLE));	\
    divide_by_nonmissings_##DOUBLE(SxI, pseudo_freq, mssngs, missings, sum, S, ans); \
    if (false) PRINTF("back to allele_freq %ld\n", (Long) sum);		\
    FREE(sum);								\
    if (false) PRINTF(" allele_freq nearly DONE\n");			\
    if (S != Sums) {FREE(S);}						\
    if (false) PRINTF(" allele_freq DONE\n");				\
  }

#if defined compatibility_to_R_h
allele_freq_fctn(double)
#endif
allele_freq_fctn(LongDouble)


LongDouble getTotalSum(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[TotalSumIndex];
}

LongDouble *getFreq(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP)) + FreqIndex;
}


LongDouble *getSum(SEXP SxI){
  Long *info = GetInfo(SxI);
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP)) + SumIndex;
}


LongDouble *getFreqSxI(SEXP SxI){
  Long *info = GetInfo(SxI);
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP)) + FreqSxIIndex;
}


LongDouble getSqTotalSum(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[SqTotalSumIndex];
}


LongDouble getSigmaSq(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[SigmaSqIndex];
}

LongDouble getSumFreq(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[SumFreqIndex];
}

LongDouble *getPseudoFreq(SEXP SxI){
  Long *info = GetInfo(SxI);
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP)) + PseudoFreqIndex;
}


LongDouble *getPseudoSum(SEXP SxI){
  Long *info = GetInfo(SxI);
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP)) + PseudoSumIndex;
}


LongDouble *getPseudoFreqSxI(SEXP SxI){
  Long *info = GetInfo(SxI);
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP)) + PseudoFreqSxIIndex;
}


LongDouble getPseudoSqTotalSum(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[PseudoSqTotalSumIndex];
}


LongDouble getPseudoSigmaSq(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[PseudoSigmaSqIndex];
}

LongDouble getPseudoSumFreq(SEXP SxI){
  SEXP PP = getAttribPointer(SxI, Precise);
  if (PP == R_NilValue) BUG;
  return (LONGREAL(PP))[PseudoSumFreqIndex];
}

  
LongDouble *getFreq(SEXP SxI, double *externalFreq,
		    option_type *global, utilsoption_type *utils){
  const bool debug = false;
  Long *info = GetInfo(SxI);
  //  printf("entering getfreq... %d %d %d longdouble=%d %d\n", (int) LENGTH(SxI), info[30], Precise==R_NilValue, (int) sizeof(LongDouble), length(Precise));

  SEXP PP = getAttribPointer(SxI, Precise);
  if (debug)
    PRINTF("entering getfreq... %ld Nil=%d \n", info[30], PP==R_NilValue);
  if (PP != R_NilValue) return (LONGREAL(PP)) + FreqIndex;
  
  //  printf("entering truely getfreq ... \n");
  // Long *info = GetInfo(SxI);
  SEXP next = getAttribPointer(SxI, Next);
  SEXP Miss = getAttribPointer(SxI, Missings);
  Long *mssngs = Miss == R_NilValue ? NULL : LONG(Miss);
  Long missings = info[MISSINGS];

  basic_options *opt = &(utils->basic);
  Long snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    lenFreq = snps,
    lenFreqSxI = individuals,
    lenPseudoFreq = individuals,
    lenPseudoFreqSxI = snps;

  PP = allocVectorExtra(LONGREALSXP, FreqIndex + 3L * (lenFreqSxI + lenFreq));
  setAttrib(SxI, Precise, PP);
  LongDouble *P = LONGREAL(PP);
  LongDouble *freq = P + FreqIndex,
    *sums = P + SumIndex,
    *freqSxI = P + FreqSxIIndex,
    *pseudofreq =  P + PseudoFreqIndex,
    *pseudosums = P + PseudoSumIndex,
    *pseudofreqSxI = P + PseudoFreqSxIIndex
    ;


  centering_type centered = global->genetics.centered;
  normalizing_type normalized = global->genetics.normalized;
  global->genetics.centered = NoCentering;
  global->genetics.normalized = NoNormalizing;

  allele_freq_LongDouble(SxI, mssngs, missings, false, global, utils,
			 sums, freq);
  if (externalFreq) {
    for (Long j=0; j<snps; j++) freq[j] = (LongDouble) externalFreq[j];
  }

 LongDouble sigmaSq = 0;
  LongDouble sumFreq = 0;
  for (Long i=0; i<lenFreq; i++) {
    sumFreq += freq[i];
    sigmaSq += freq[i] * (1.0 - freq[i]);    
  }
  P[SigmaSqIndex] = 2.0 * sigmaSq;
  P[SumFreqIndex] = sumFreq;


  
  // printf(" ** getFreq nach allelefreq \n");
  vectorGeno_raw_LongDouble(SxI, freq, 1, lenFreq, global, utils,
			    freqSxI, lenFreqSxI);
 //  printf(" ** getFreq nach vectorGeno \n");

   LongDouble freqSxIOne = 0.0;
  for (Long i=0; i<lenFreqSxI; i++) freqSxIOne += freqSxI[i];
  P[SqTotalSumIndex] = (LongDouble) ((Long) (2.0 * freqSxIOne * lenFreqSxI
					       + 0.1));

  if (next == NULL) {
    Ulong *s = (Ulong*) MALLOC(sizeof(Ulong) * lenPseudoFreq);
    Ulong sumgeno = sumGeno(Align(SxI, opt), snps, individuals,
			    info[LDA], (coding_type)info[CODING],
			    (int) info[VARIANT],
			    GreaterZero(opt->cores), s);
    P[TotalSumIndex] = (LongDouble) sumgeno;
    divide_by_nonmissings_LongDouble(SxI, true, mssngs, missings, s,
				     pseudosums, pseudofreq);
    FREE(s);
  } else {
     LongDouble freqOne = 0.0;
    for (Long i=0; i<lenFreq; i++) freqOne += freq[i]; 
    P[TotalSumIndex] = truncl(2.0L * freqOne * lenFreqSxI + 0.1L);

    allele_freq_LongDouble(next, mssngs, missings, true, global, utils,
			   pseudosums, pseudofreq);
  }
 // printf("near end getfreq\n"); 
 
  LongDouble pseudosigmaSq = 0;
  LongDouble pseudosumFreq = 0;
  for (Long i=0; i<lenPseudoFreq; i++) {
    pseudosumFreq += pseudofreq[i];
    pseudosigmaSq += pseudofreq[i] * (1.0 - pseudofreq[i]);
  }
  P[PseudoSigmaSqIndex] = 2.0 * pseudosigmaSq;
  P[PseudoSumFreqIndex] = pseudosumFreq;

  genoVector_raw_LongDouble(SxI, pseudofreq, 1, lenPseudoFreq, global, utils,
			    pseudofreqSxI, lenPseudoFreqSxI);

  LongDouble pseudofreqSxIOne = 0.0;
  for (Long i=0; i<lenPseudoFreqSxI; i++) pseudofreqSxIOne += pseudofreqSxI[i];
  P[PseudoSqTotalSumIndex] = truncl(2.0L * pseudofreqSxIOne * lenPseudoFreqSxI
				    + 0.1L);

  
  global->genetics.centered = centered;
  global->genetics.normalized = normalized;

  //  printf("leaving getfreq\n"); 
 
  return freq;
}

LongDouble *getFreq(SEXP SxI, option_type *global, utilsoption_type *utils){
  return getFreq(SxI, NULL, global, utils);
}



coding_type switch_transposed(coding_type coding) {
 switch(coding) {
 case TwoBitGeno : return(TwoBitGenoTransposed);
 case TwoBitGenoTransposed : return(TwoBitGeno);
 case FiveCodes : return(FiveCodesTransposed);
 case FiveCodesTransposed : return(FiveCodes);
 case Plink : return(PlinkTransposed);
 case PlinkTransposed : return(Plink);
 case OrigPlink : return(OrigPlinkTransposed);
 case OrigPlinkTransposed : return(OrigPlink);
 default : BUG;
 }
}



void sparseTGeno(Uchar *code,
		 Long nrow_code,
		 Long ncol_code,
		 Long lda, // in Byte or unit_t, depending on coding
		 coding_type coding,
		 double *valueB,
		 int nIdx, // 2nd dimension of sparse; == length(rowIdxB)
		 int *rowIdxB,
		 int *colIdxB,
		 bool tSparse,
		 option_type *global, utilsoption_type *utils,
		 double *C,
		 Long Ldc) {
  //   printINFO(SxI);
  MEMSET(C, 0, sizeof(*C) * Ldc * nrow_code);
  if (!has_lda_in_Byte(coding)) lda *= sizeof(unit_t);
  if (tSparse) BUG;
  sparseTGeno_t sTG = NULL;
  switch(coding) {
  case Plink : case OrigPlink : sTG = sparseTGenoPlink; break;    
  case TwoBitGeno : sTG = sparseTGeno2Bit; break;
  case OneByteGeno : sTG = sparseTGenoOneByte; break;
  default : BUG;
  }

  sTG(code, nrow_code, ncol_code, lda, coding,
      valueB, nIdx, rowIdxB, colIdxB, tSparse,
      global, utils, C, Ldc);
  
}


void sparseTGeno(SEXP SxI,
		bool t_SxI,
		double *valueB,
		int nIdx, // 2nd dimension of sparse; == length(rowIdxB)
		int *rowIdxB,
		int *colIdxB,
		int tSparse,
		 option_type *global, utilsoption_type *utils,
		double *C,
		Long Ldc) {
  if (!t_SxI) SxI = getAttribPointer(SxI, Next);
  if (SxI == R_NilValue) BUG;
  Long *info = GetInfo(SxI);
    
  sparseTGeno((Uchar*) Align(SxI, NULL),
	      info[t_SxI ? INDIVIDUALS : SNPS],
	      info[t_SxI ? SNPS : INDIVIDUALS],
	      info[LDA], (coding_type) info[CODING],
	      valueB, nIdx, rowIdxB, colIdxB, tSparse,
	      global, utils,
	      C, Ldc);
}
