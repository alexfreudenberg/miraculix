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
#include "MXinfo.h"
#include "intrinsicsCheck.h"
#include "options.h"
#include "intrinsics_specific.h"





Uint Inti(SEXP X, Uint i) {
  switch(TYPEOF(X)) {
  case INTSXP : return (Uint) INTEGER(X)[i];
  case LGLSXP : return (Uint) LOGICAL(X)[i];
  case REALSXP : return (Uint) REALX(X)[i];
  default : ERR0("not of numerical type");
  }
}




Long *GetInfoUnchecked(SEXP Code) { // has negative integers included!!
  SEXP Info = getAttribPointer(Code, Information);
  if (Info == R_NilValue) return NULL;
  //printf("length Infos=%d %d\n", LENGTH(Infos), TYPEOF(Infos));
#if defined compatibility_to_R_h
  if (TYPEOF(Info) != INTSXP)
#elif defined compatibility_to_C_h
    if (TYPEOF(Info) != LONGSXP)      
#endif
      ERR1("Obsolete storage mode %d", TYPEOF(Info));
  return LONG(Info);
}

Long  snps2entities(Long snps,  Long CodesPerEntity) {
  return 1L + (snps - 1L) / CodesPerEntity;
}			       


Long *GetInfo(SEXP Code, int nullOK) { // has negative integers included!!
  Long *info = GetInfoUnchecked(Code);
  if (info == NULL || info[CODING] == NA_LONG) {//Code=filename if NA_LONG
    if (nullOK) return info == NULL || nullOK == 1 ? NULL : info;
    BUG;
  }
  basic_options opt;
  Ext_get_utils_basic(&opt, false);
  coding_type coding = (coding_type) info[CODING];
  // int variant =
    check_variant(coding, (int) info[VARIANT], opt.efficient);
  if (info[IMPLEMENTATION] !=CURRENT_IMPLEMENTATION)
    ERR2("the stored data format (%ld) does not match the current format (%d).",
	 info[IMPLEMENTATION], CURRENT_IMPLEMENTATION );
  return (Long*) info;
}

Long *GetInfo(SEXP Code) {
   return GetInfo(Code, false);
}

void printINFO(SEXP S) {
  // PRINTF("info:\n");
  const char *eq[]={"!=NULL", "=0"};
  SEXP Info = getAttribPointer(S, Information);
  PRINTF("SEXP: x%s len=%ld, dim=(%d,%d) type=%d dim=%d, Info%s File(%s,%s) prec=%s Class%s Next%s\n",
	 eq[VOIDSXP(S)  == R_NilValue],
	 (Long) (LENGTH(S)),
	 (int) NROWS(S), (int) NCOLS(S),
	 TYPEOF(S),
	 getDimension(S),
	 eq[Info == R_NilValue],
	 eq[getAttribPointer(S, Filecoding)==R_NilValue],
	 eq[getAttribPointer(S, Filename)==R_NilValue],
	 eq[getAttribPointer(S, Precise)==R_NilValue],
	 eq[GET_CLASS(S) == R_NilValue],
	 eq[getAttribPointer(S, Next)==R_NilValue]);
  if (Info != R_NilValue) {
    //    PRINTF("info::\n");
    Long *info = LONG(Info);
    for (int K=0;
	 //	 K<=INFO_GENUINELY_LAST
	 K<INFO_LAST
	   ; K++) {
      //   PRINTF("K=%d\n", K);
      if (K >= 10 && K <= 20 && K!=MEMinUNITS) continue;
      if (K >= 40 && K <= 50) continue;
      if (STRCMP("unused", INFO_NAMES[K])) {
	PRINTF("%2d: %s=\t", K, INFO_NAMES[K]);
        if (STRLEN(INFO_NAMES[K]) <= 10) PRINTF("\t");
	if (K==MEMinUNITS) PRINTF("%lu", info[MEMinUNITS]); else
	  if (info[K]==NA_LONG) PRINTF("NA"); else
	    if (info[K]==0) PRINTF("FALSE"); else // OK
	      if (info[K]==1) PRINTF("TRUE"); else PRINTF("%lu", info[K]); // OK
	switch(K) {
	case 3: PRINTF("\t[%s]\n", CODING_NAMES[info[K]]); break;
	default: PRINTF("\n");
	}
      }
    }
    PRINTF("\n");
  } else PRINTF("Info = R_NilValue!\n");
  if (getAttribPointer(S, Next) != R_NilValue) {
    PRINTF("NEXT:\n");
    printINFO(getAttribPointer(S, Next));
  }
}

  


// LDAbitalign == 0 requests the primitive alignment for Align
 // LDAbitalign > 0 calculates the space needed for a single individual
 //     in snp x individual storage
Long fctnLDAalignBitSnps(Long bits, Long LDAbitalign, Long my_LDABITALIGN,
		    coding_type my_CODING) {
  //  printf("ldabit: %ld %ld; bits=%ld\n", LDAbitalign, my_LDABITALIGN, bits);
  if (LDAbitalign != 0 && LDAbitalign < my_LDABITALIGN)		       
    ERR3("bit alignment for '%s' (%lu) too small (not >= %lu)",	      
	 CODING_NAMES[my_CODING], LDAbitalign, my_LDABITALIGN);
  Long newLDAbitalign = LDAbitalign ? LDAbitalign : my_LDABITALIGN;    
  assert(newLDAbitalign >= BitsPerUnit && newLDAbitalign % BitsPerUnit ==0);// OK
  return(ROUND_GEQ(bits, newLDAbitalign) / BitsPerUnit); // OK
}


