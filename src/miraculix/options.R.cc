
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
#include "MXinfo.h"
#include "options.h"
#include "mmagpu.h"
#include "utils_miraculix.h"
#include "kleinkram.h"

  

SEXP existsVariant(SEXP Coding, SEXP Variant, SEXP Indeed) {
  
  SEXP Ans; 
  coding_type coding = (coding_type) INTEGER(Coding)[0];
  bool indeed = LOGICAL(Indeed)[0];
  if (LENGTH(Variant) == 0) {
    int zaehler = 0,
      curr[NvariantsSHR5];
    for (int i=1; i<NvariantsSHR5; i++) { // nicht die NULL!
      int v = i << 5;
      if (exists_variant(coding, v, indeed, false))
	curr[zaehler++] = v;
    }
    PROTECT(Ans = allocVector(INTSXP, zaehler));
    int *ans = INTEGER(Ans);
    for (int i=0; i<zaehler; i++) ans[i] = curr[i];
  } else {
    PROTECT(Ans = allocVector(LGLSXP, 1));
    LOGICAL(Ans)[0] = exists_variant(coding, INTEGER(Variant)[0], indeed,
				     indeed);
  }
  
  UNPROTECT(1);
  return Ans;
}



SEXP existsTiling(SEXP Coding, SEXP Variant, SEXP tiling) {
  // function assumes that combination coding/variant is valid!
  //  option_type *global = &(KEYT_M()->global);
  SEXP Ans;
  int m = GetName(Coding, (char*) "coding", CODING_NAMES, LastUserCoding + 1);
  if (m<=0 ||  m > LastUserCoding) ERR0("given snp coding not allowed.");
  coding_type coding = (coding_type) m;
  int variant = INTEGER(Variant)[0];
 
  if (LENGTH(tiling) == 0) {
    int zaehler = 0,
      curr[2];
    for (int i=0; i<=1; i++) {
      //printf("i=%d V = %d\n", i,  variants[i]);
      if (exists_tiling(coding, variant, (bool) i)) {
	//printf("     %d %d %d\n", i,    zaehler, variants[i]);
	curr[zaehler++] = i;
      }
    }
    PROTECT(Ans = allocVector(LGLSXP, zaehler));
    int *ans = LOGICAL(Ans);
    for (int i=0; i<zaehler; i++) ans[i] = curr[i];
  } else {
    PROTECT(Ans = allocVector(LGLSXP, 1));
    LOGICAL(Ans)[0] = exists_tiling(coding, variant, LOGICAL(tiling)[0]);
  }
  UNPROTECT(1);
  return Ans;
}



SEXP existsAllelefreq(SEXP Coding) {
  SEXP Ans;
  if (LENGTH(Coding) == 0) {
    int zaehler = 0;
    for (int i=FirstGenoCoding; i<=LastHaploCoding; i++)
      zaehler += exists_allelefreq((coding_type) i);
    Ans = PROTECT(allocVector(INTSXP, zaehler));
    int *ans = INTEGER(Ans);
    zaehler = 0;
    for (int i=FirstGenoCoding; i<=LastHaploCoding; i++)
      if (exists_allelefreq((coding_type) i)) ans[zaehler++] = i;
  } else {
    Ans = PROTECT(allocVector(LGLSXP, 1));
    int m = GetName(Coding, (char*) "coding", CODING_NAMES, nr_coding_type);
    LOGICAL(Ans)[0] = exists_allelefreq((coding_type) m);
  }
  UNPROTECT(1);
  return Ans;
}


SEXP existsCoding(SEXP Coding, SEXP Internal) {
  int internal = (int) LOGICAL(Internal)[0];
  exists_coding_t exists = (internal == NA_INTEGER) ? exists_coding :
    internal ? exists_internal_coding : exists_user_coding;
  SEXP Ans;
  if (LENGTH(Coding) == 0) {
    int zaehler = 0;
    for (int i=FirstGenoCoding; i<=LastHaploCoding; i++)
      zaehler += exists((coding_type) i);
    Ans = PROTECT(allocVector(INTSXP, zaehler));
    int *ans = INTEGER(Ans);
    zaehler = 0;
    for (int i=FirstGenoCoding; i<=LastHaploCoding; i++)
      if (exists((coding_type) i)) ans[zaehler++] = i;
  } else {
    Ans = PROTECT(allocVector(LGLSXP, 1));
    int m = GetName(Coding, (char*) "coding", CODING_NAMES, nr_coding_type);
    LOGICAL(Ans)[0] = exists((coding_type) m);
  }
  UNPROTECT(1);
  return Ans;
}


SEXP existsCrossprod(SEXP Coding) {
  SEXP Ans;
  if (LENGTH(Coding) == 0) {
    int zaehler = 0;
    for (int i=FirstGenoCoding; i<=LastGenoCoding; i++)
      zaehler += exists_crossprod((coding_type) i);
    Ans = PROTECT(allocVector(INTSXP, zaehler));
    int *ans = INTEGER(Ans);
    zaehler = 0;
    for (int i=FirstGenoCoding; i<=LastGenoCoding; i++)
      if (exists_crossprod((coding_type) i)) ans[zaehler++] = i;
  } else {
    Ans = PROTECT(allocVector(LGLSXP, 1));
    int m = GetName(Coding, (char*) "coding", CODING_NAMES, nr_coding_type);
    LOGICAL(Ans)[0] = exists_crossprod((coding_type) m);
  }
  UNPROTECT(1);
  return Ans;
}




void setMoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		bool isList, bool local) {

   
  if (!local && Ext_parallel())
    ERR0("'RFoptions' may not be set from a parallel process.");
      
  option_type *options;
  utilsoption_type *utils;
  WhichOptionList(local, &options, &utils);
  
 //  printf("i=%d %d %d %d %ld\n", i, j, isList, local, options);

  switch(i) {
  case 0: {// genetics
    genetics_options *gp = &(options->genetics);
    switch(j) {
    case 0: gp->digits = NUM; break;
    case 1: {
      Rint m= GetName(el, name, CODING_NAMES, LastUserCoding + 1,
		      gp->coding);
      if ((gp->coding != m && !isList) &&
	  (m < 0 || m > LastUserCoding ||
	   m == FourBit || m == TwoByte ||
	   m == FourByte_Tr ||  m == unused6 || m == unused7 ||
	   m == unused9 || m == unused10
	   )) {
	ERR1("given snp coding (%s) not allowed. (Note that names have changed from version V.0 to V.1 of miraculix.)", CODING_NAMES[m]);
      }
      if (gp->coding != m) {
	options->tuning.variant = 0;
	gp->coding = (coding_type) m;
      }
    }
      break;
    case 2:
      if (TYPEOF(el) == LGLSXP || TYPEOF(el) == STRSXP) {
	gp->centered =  (centering_type) GetName(el, name, CENTERING_NAMES,
						 nr_centering_type,
						 gp->centered);
	if (gp->centered == User) {
	  if (gp->pcentered == NULL) {
	    WARN1("'%.50s' set to 'row means'", name); // OK
	    gp->centered = RowMeans;
	    gp->normalized = AccordingCentering;
	  }
	} else {
	  FREE(gp->pcentered);
	  gp->ncentered = 0;
	  gp->normalized = AccordingCentering;
	}
      } else { // numeric el
	int len = LENGTH(el);
	gp->ncentered = len;
	FREE(gp->pcentered);
	gp->pcentered = (double*) MALLOC(len * sizeof(*(gp->pcentered)));
	Real(el, name, gp->pcentered, len);
	gp->centered = User;
	gp->normalized = AccordingCentering;
      }
      break;      
    case 3:
      gp->normalized =  (normalizing_type) GetName(el, name, NORMALIZING_NAMES,
						   nr_normalizing_type,
						   gp->normalized);
      break;
    case 4 : gp->squared = LOGI; break;
    case 5 : gp->haplo_set1 = USRLOG; break;
    case 6 : gp->interprete_coding_as_is = LOGI; break;
    case 7 : gp->prefer_external_freq = LOGI; break;
   default: BUG; 
    }}
    break;
  case 1: {
    tuning_options *gp = &(options->tuning);
    switch(j) {
    case 0 : gp->savegeno = LOGI; break;
    case 1 : gp->oldversion = LOGI; break;
    case 2 : {
      int m = POS0INT;
#if ! defined USEGPU
      if (m == VARIANT_GPU) {
 	ERR0("'CUDA >= 11' is not available under the current compilation. See the starting message for a remedy.");	
      }
#endif
      //      printf("set options %s %d local=%d %d %ld\n",	     CODING_NAMES[options->genetics.coding],	     m, local, utils->basic.efficient, utils);
      
      if (!exists_variant(options->genetics.coding, m, utils->basic.efficient,
			  true))
	ERR2("required variant '%d' does not match snpcoding='%.50s'.\n",
	     m, CODING_NAMES[options->genetics.coding]);
      gp->variant = m;
    }
      break;
    case 3 : {
      gp->logSnpGroupSize = POS0INT;
      //      printf("option=%d\n", gp->logSnpGroupSize);
    }
      break;
    case 4 : gp->indivGroupSize = POS0INT; break;
    case 5 : gp->smoothSnpGroupSize = LOGI; break;
    case 6 : gp->meanVsubstract = LOGI; break;
    case 7 : gp->meanSxIsubstract = LOGI; break;
    case 8 : gp->floatLoop = POS0INT; break;
    case 9 : gp->missingsFully0 = LOGI; break;
    case 10 : gp->miniRows = POS0INT; break;
    case 11 : gp->miniCols = POS0INT; break;
    case 12 : gp->gpu = LOGI; break;
    case 13 : gp->addtransposed = LOGI; break;
   
   default: BUG; 
      break;
    }
  }
    break;
  case 2: {
    messages_options *gp = &(options->messages);
    switch(j) {
    case 0: gp->warn_address = LOGI; break;
    default: BUG;
    }}
    break;    
  default: BUG; 
  }
  //  printf("end set\n");
}

   

void getMoptions(SEXP sublist, int i, bool local) {
  int k = 0;
  option_type *options;
  utilsoption_type *utils;
   WhichOptionList(local, &options, &utils);
 
  switch(i) {
  case 0 : {
    genetics_options *p = &(options->genetics);
    ADD(ScalarReal(p->digits));
    //    ADD(ScalarString(mkChar(RELSHIP_METH_NAME[p->coding])));
    ADD(ScalarString(mkChar(CODING_NAMES[p->coding])));
    ADD(ScalarString(mkChar(CENTERING_NAMES[p->centered])));
    ADD(ScalarString(mkChar(NORMALIZING_NAMES[p->normalized])));
     ADD(ScalarLogical(p->squared));
    ADD(ExtendedBooleanUsr(p->haplo_set1));
    ADD(ScalarLogical(p->interprete_coding_as_is));
    ADD(ScalarLogical(p->prefer_external_freq));
   }
    break;
  case 1 : {
    tuning_options *p = &(options->tuning);
    ADD(ScalarLogical(p->savegeno));
    ADD(ScalarLogical(p->oldversion));
    ADD(ScalarInteger(p->variant));
    ADD(ScalarInteger(p->logSnpGroupSize));
    ADD(ScalarInteger(p->indivGroupSize));
    ADD(ScalarLogical(p->smoothSnpGroupSize));
    ADD(ScalarLogical(p->meanVsubstract));
    ADD(ScalarLogical(p->meanSxIsubstract));
    ADD(ScalarInteger(p->floatLoop));
    ADD(ScalarLogical(p->missingsFully0));
    ADD(ScalarInteger(p->miniRows));
    ADD(ScalarInteger(p->miniCols));
    ADD(ScalarLogical(p->gpu));
    ADD(ScalarLogical(p->addtransposed));
    //     ADD(ScalarLogical(p->));
   }
    break;
  case 2 : {
    messages_options *p = &(options->messages);
    ADD(ScalarLogical(p->warn_address));
  }
    break;   
  
  default : BUG;
  }
}

#endif

