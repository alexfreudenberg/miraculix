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


#include "Basic_RFUlocal.h"

#if defined compatibility_to_R_h

#include "RandomFieldsUtils.h"
#include "zzz_RFU.h"
#include "xport_import_RFU.h"
#include "kleinkram.h"
#include "extern_RFU.h"
extern const char *basic[basicN];



AVAILABLE_SIMD

//typedef
struct getlist_type {
  int ListNr, i;
}; // getlist_type;


void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);


#define MAXNLIST 7
#define PKGNAMELENGTH 20
int NList = 0; 
int noption_class_list = 0,
  AllprefixN[MAXNLIST] = {0, 0, 0, 0, 0, 0, 0},
  Allversion[MAXNLIST] = {0, 0, 0, 0, 0, 0, 0},
  *AllallN[MAXNLIST] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
const char  *option_class_list[MAXNLIST] =
  {prefixlist[1], NULL, NULL, NULL, NULL, NULL},
  **Allprefix[MAXNLIST] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
  ***Allall[MAXNLIST] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
name_type pkgnames = {"", "", "", "", "", "", ""}; // MAXNLIST !!!
setoptions_fctn setoption_fct_list[MAXNLIST][MAXNLIST] = 
  {{setoptions, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL}};
getoptions_fctn getoption_fct_list[MAXNLIST][MAXNLIST] = 
  {{getoptions, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
   {NULL, NULL, NULL, NULL, NULL, NULL, NULL}};
 finalsetoptions_fctn finaloption_fct_list[MAXNLIST] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL };
deleteoptions_fctn deloption_fct_list[MAXNLIST] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL};
bool installed [MAXNLIST] = { false, false, false, false, false, false, false};
install_modes min_simd_needs[MAXNLIST] = // currently only used if == none or not
  {Inone, Inone, Inone, Inone, Inone, Inone, Inone},
  min_gpu_needs[MAXNLIST] = // currently not used
  {Inone, Inone, Inone, Inone, Inone, Inone, Inone};
Uint simd_infos [MAXNLIST] = {0, 0, 0, 0, 0, 0, 0};
  
/* precide information of min_install is not used yet; just whether it is 
   'none' or not */

void hintVariable(char *name, int warn_unknown_option) {
  static bool printing = true; // da nur Mutterprozess schreiben darf
  if (warn_unknown_option > 0 && OPTIONS.basic.Rprintlevel > 0) {
    PRINTF("'%s' is considered as a variable (not as an option).\n", name);
    if (printing && OPTIONS.basic.helpinfo && !parallel()) {
      PRINTF("[This hint can be turned off by 'RFoptions(%s=-%d)'.]\n",
	     basic[INSTALL_RUN_WARN_OPTION], MAX(1, warn_unknown_option));
      printing = false;
    }
  }
}

void setparameter(SEXP el, char *prefix, char *mainname, bool isList,
		  getlist_type *getlist, int warn_unknown_option, bool local,
		  int calling) {
  int 
    j = NOMATCHING,
    i = NOMATCHING,
    ListNr = NOMATCHING;
  char name[LEN_OPTIONNAME];

  SPRINTF(name, "%.50s%.50s%.50s", prefix, STRLEN(prefix) == 0 ? "" : ".",
	  mainname);

  if (STRCMP(prefix, "")) {
    for (ListNr=0; ListNr<NList; ListNr++) {
      i = Match(prefix, Allprefix[ListNr], AllprefixN[ListNr]);
      if (i != NOMATCHING) break;
    }
    if (i == NOMATCHING) ERR1("option prefix name '%.50s' not found.", prefix); 
    if (i < 0 || STRCMP(prefix, Allprefix[ListNr][i])) {
      for (int k=ListNr + 1; k < NList; k++) {
	int ii = Match(prefix, Allprefix[ListNr], AllprefixN[ListNr]);
	if (ii == NOMATCHING) continue;
	i = MULTIPLEMATCHING;
	if (ii >= 0 && STRCMP(prefix, Allprefix[k][ii]) == 0) {
	  ListNr = k;
	  i = ii;
	  break;
	} // ii >0
      } // for k
      if (i == MULTIPLEMATCHING) 
	ERR1("option prefix name '%.50s' is ambiguous.", prefix);   
    } // prefix == List

    j = Match(mainname, Allall[ListNr][i], AllallN[ListNr][i]);
      
  } else { // (i==0), no prefix given
#define MinNameLength 3
    for (ListNr=0; ListNr<NList; ListNr++) {
      int trueprefixN = AllprefixN[ListNr];
      for (i=0; i < trueprefixN; i++) {	
	j = Match(mainname, Allall[ListNr][i], AllallN[ListNr][i]);
	if (j != NOMATCHING) break;     
      }
      if (j != NOMATCHING) break;
    }
    if (j==NOMATCHING) {

      switch(warn_unknown_option) {
      case WARN_UNKNOWN_OPTION_NONE1 :
	hintVariable(name, warn_unknown_option);
 	FALLTHROUGH_OK;
     case WARN_UNKNOWN_OPTION_NONE : case -WARN_UNKNOWN_OPTION_NONE1 :
	return;
      case WARN_UNKNOWN_OPTION_SINGLE : case -WARN_UNKNOWN_OPTION_SINGLE :
	if (STRLEN(name) != 1) ERR1("Unknown option '%.50s'.", name);
	FALLTHROUGH_OK;
      case WARN_UNKNOWN_OPTION_CAPITAL : case -WARN_UNKNOWN_OPTION_CAPITAL :
	if (name[0] >= 'A' && name[0] <= 'Z') {
	  hintVariable(name, warn_unknown_option);
	  return;
	}
	ERR1("Unknown option: '%.50s'.", name); 
      case WARN_UNKNOWN_OPTION_ALL :
      default :
	ERR1("unknown option '%.50s'.", name);
      }
      return;
    }
    
    // printf("j=%d  %.50s\n", j,  j >=0 ? Allall[ListNr][i][j] : "multi");
    
    if (j < 0  || STRCMP(mainname, Allall[ListNr][i][j])) {
      int starti = i + 1;
      for (int k = ListNr; k<NList; k++, starti=0) {
	int tpN = AllprefixN[k];
	for (int ii=starti; ii < tpN; ii++) {	
	  int jj = Match(mainname, Allall[k][ii], AllallN[k][ii]);
	  if (jj == NOMATCHING) continue;
	  
	  //  printf("listnr=%d %.50s jj=%d %.50s\n", ListNr, Allall[ListNr][i][j], 
	  //	 jj, jj < 0 ? "no" : Allall[k][ii][jj]);
	  j = MULTIPLEMATCHING;
	  if (jj >= 0 && STRCMP(mainname, Allall[k][ii][jj]) == 0) {
	    ListNr = k;
	    i = ii;
	    j = jj;
	    break;
	  } // jj
	} // for ii
	if (j >= 0) break;
      } // for k
    } // if j < 0 || !=
  } // no prefix given

  if (j<0) {
    if (j == NOMATCHING) ERR1("option not recognized: '%.50s'.", name)
      else ERR1("Multiple partial matching for option '%.50s'.", name);
  }
 
  if (getlist != NULL) {
    int k=0;
    while((getlist[k].ListNr != ListNr || getlist[k].i != i)
	  && getlist[k].ListNr >= 0) k++;
    if (getlist[k].ListNr < 0)
      ERR2("Option '%.50s' not allowed for this call.\n   In case you really need this option, use the command 'RFoption(%.50s=...)'", mainname, mainname);
  }

  setoptions_fctn *set = setoption_fct_list[ListNr];
  if (local && set[calling]!= NULL)
    set[calling](i, j, el, name, isList, local);
  else set[ListNr](i, j, el, name, isList, local);
  
  
}




void getListNr(bool save, int t, int actual_nbasic, SEXP which,
	      getlist_type *getlist,
	      int *Nr, int *idx // output
	      ){
  int i, ListNr;
  const char *z;
  if (save && t < noption_class_list) z = option_class_list[t];
  else z = (char*) CHAR(STRING_ELT(which, t - actual_nbasic));
  for (ListNr=0; ListNr<NList; ListNr++) {
    int PreFixN = AllprefixN[ListNr];
    for (i=0; i<PreFixN; i++)
      if (STRCMP(Allprefix[ListNr][i], z) == 0) break;
    if (i < PreFixN) break;
  }
  if (ListNr >= NList) ERR0("unknown value for 'getoptions_'");
  if (getlist != NULL) {
    getlist[t].ListNr = ListNr;
    getlist[t].i = i;
  }
  *Nr = ListNr;
  *idx = i;
}

 

SEXP getRFUoptions(int ListNr, int i, bool local, int calling) {
  SEXP sublist, subnames;
  int elmts = AllallN[ListNr][i];
  PROTECT(sublist = allocVector(VECSXP, elmts));
  PROTECT(subnames = allocVector(STRSXP, elmts));
  for (int k=0; k<elmts; k++) {
    //    printf("getopt %d k=%d elmnts=%d %d %.50s, local=%d\n", i, k, elmts, ListNr, Allall[ListNr][i][k], local);
    SET_STRING_ELT(subnames, k, mkChar(Allall[ListNr][i][k])); 
  }
  getoptions_fctn *get = getoption_fct_list[ListNr];
  if (local && get[calling] != NULL) get[calling](sublist, i, local);
  else get[ListNr](sublist, i, local);
  setAttrib(sublist, R_NamesSymbol, subnames);
  UNPROTECT(2);
  return sublist;
}


SEXP getRFUoptions(SEXP which, getlist_type *getlist, bool save, bool local,
		   int calling) {
  
  int ListNr, idx,
    actual_nbasic = noption_class_list * save,
    totalN = length(which) + actual_nbasic;  

  if (totalN == 0) return R_NilValue;
  if (totalN == 1) {
    getListNr(save, 0, actual_nbasic, which, getlist, &ListNr, &idx);
    return getRFUoptions(ListNr, idx, local, calling);
  }

  SEXP list, names;
  PROTECT(list = allocVector(VECSXP, totalN));
  PROTECT(names = allocVector(STRSXP, totalN));
  for (int t=0; t<totalN; t++) {
    getListNr(save, t, actual_nbasic,  which, getlist, &ListNr, &idx);
    SET_VECTOR_ELT(list, t, getRFUoptions(ListNr, idx, local, calling));
    SET_STRING_ELT(names, t, mkChar(Allprefix[ListNr][idx]));    
  }
  setAttrib(list, R_NamesSymbol, names);     
  UNPROTECT(2);
  return list;
}


SEXP getRFUoptions(bool local, int calling) {
  SEXP list, names;
 
  int  PreFixN, totalN, i, ListNr,
    itot = 0;
  for (totalN=ListNr=0; ListNr<NList; ListNr++) {
    //    printf("totalN=%d\n", totalN);
    PreFixN = AllprefixN[ListNr];
    totalN += PreFixN;
   }

  PROTECT(list = allocVector(VECSXP, totalN));
  PROTECT(names = allocVector(STRSXP, totalN));

  for (ListNr =0; ListNr<NList; ListNr++) {
    //    printf("ListNr %d\n", ListNr);
    PreFixN = AllprefixN[ListNr];    
    for (i=0; i<PreFixN; i++) {
      SET_VECTOR_ELT(list, itot, getRFUoptions(ListNr, i, local, calling));
      SET_STRING_ELT(names, itot, mkChar(Allprefix[ListNr][i]));
      itot ++;
    } 
  }

  setAttrib(list, R_NamesSymbol, names);   
  UNPROTECT(2);

  return list;
}



void setRFUoptions(SEXP el, char *name, bool isList, getlist_type *getlist,
		   int warn_unknown_option, bool local, int calling) {
  int i, len;
  char PreFix[LEN_OPTIONNAME / 2], mainname[LEN_OPTIONNAME / 2];   
  //  printf("splitandset\n");
  len = STRLEN(name);
  for (i=0; i < len && name[i]!='.'; i++);
  if (i == 0) { ERR1("argument '%.50s' not valid\n", name); }
  if (i==len) {
    STRCPY(PreFix, "");
    strcopyN(mainname, name, LEN_OPTIONNAME / 2);
  } else {
    strcopyN(PreFix, name, MIN(i + 1, LEN_OPTIONNAME / 2));
    strcopyN(mainname, name+i+1, MIN(STRLEN(name) - i, LEN_OPTIONNAME / 2) );
  }

  // 
  //  printf("i=%d %d %.50s %.50s local=%d\n", i, len, PreFix, mainname, local);
  setparameter(el, PreFix, mainname, isList && OPTIONS.basic.asList, getlist,
	       warn_unknown_option, local, calling);
  //   printf("ende\n");
}

SEXP RFoptions(SEXP options) {
  return RFUoptions(options, (char*) "RandomFieldsUtils");
}


  
SEXP RFUoptions(SEXP options, char *pkgname) {
  int calling = 0;
  for ( ; calling < NList; calling++)
    if (!STRCMP(pkgname, pkgnames[calling])) break;
  if (calling >= NList) BUG;
  
  //printf("hier\n");
  int i, j, lenlist, lensub;
  SEXP el, list, sublist, names, subnames,
     ans = R_NilValue;
  char *name, *pref;
  bool isList = false;
  int
    warn_unknown_option = OPTIONS.installNrun.warn_unknown_option,
    n_protect = 0;
  bool local = false;

    /* 
     In case of strange values of a parameter, undelete
     the comment for PRINTF
  */

  options = CDR(options); /* skip 'name' */

  // printf(" %d %d\n", OPTIONS.installNrun.warn_unknown_option, options == R_NilValue);
  name = (char*) "";
  if (options != R_NilValue) {
    if (!isNull(TAG(options))) name = (char*) CHAR(PRINTNAME(TAG(options)));
    if (STRCMP(name, "local_") == 0	) {
      el = CAR(options);
      local = (bool) INT;
      options = CDR(options); /* skip 'name' */
     }
  }
  //   printf("hierA\n");
  if (options == R_NilValue || STRCMP(name, "") == 0)
    return getRFUoptions(local, calling);
 
  if (!isNull(TAG(options))) name = (char*) CHAR(PRINTNAME(TAG(options)));
  if (STRCMP(name, "warnUnknown_") == 0) {
    el = CAR(options);
    warn_unknown_option = INT;
    options = CDR(options); /* skip 'name' */
   }
  // printf("hierB\n");
  if (!isNull(TAG(options))) name = (char*) CHAR(PRINTNAME(TAG(options)));
  if ((isList = STRCMP(name, "list_") == 0 )) {
    //     printf("hierC\n");
    if (local) ERR0("'list_' can be used only globally.");
    list = CAR(options);
    if (TYPEOF(list) != VECSXP)
      ERR1("'list_' needs as argument the output of '%.50s'", RFOPTIONS);
    PROTECT(names = getAttrib(list, R_NamesSymbol));
    lenlist = length(list);

    if (lenlist > 0 && !local && parallel())
      ERR0("Global 'RFoptions' such as 'cores', 'seed' and 'printlevel', may be set only outside any parallel code. See '?RandomFieldsUtils::RFoptions' for the complete list of global 'RFoptions'");

   for (i=0; i<lenlist; i++) {
     //     printf("i=%d %d\n", i, lenlist);
      int len;
      pref = (char*) CHAR(STRING_ELT(names, i));  
      sublist = VECTOR_ELT(list, i);
      len = STRLEN(pref);
      for (j=0; j < len && pref[j]!='.'; j++);
      if (TYPEOF(sublist) == VECSXP && j==len) { // no "."
	// so, general parameters may not be lists,
	// others yes
	lensub = length(sublist);
	PROTECT(subnames = getAttrib(sublist, R_NamesSymbol));
	for (j=0; j<lensub; j++) {
	  name = (char*) CHAR(STRING_ELT(subnames, j));
	  // printf("Lxx '%s' %d %d %d local=%d\n", name, isList, i,j, local);
	  setparameter(VECTOR_ELT(sublist, j), pref, name, 
		       isList & OPTIONS.basic.asList, NULL, warn_unknown_option,
		       local, calling);
	}
	UNPROTECT(1);
      } else {   
	setRFUoptions(sublist, pref, isList, NULL, warn_unknown_option, local,
		      calling);
      }
    }
    UNPROTECT(1);
  } else {    // not a list
    // printf("hierD\n");
    getlist_type *getlist = NULL;
    bool save;
    if ((save = STRCMP(name, "saveoptions_") == 0
	 ) ||
	STRCMP(name, "getoptions_") == 0
	) {
      SEXP getoptions = CAR(options);
      options = CDR(options);
      if (options != R_NilValue) {
	//	printf("hier\n");
	int len = length(getoptions) + noption_class_list * save;
	getlist = (getlist_type *) MALLOC(sizeof(getlist_type) * (len + 1));
	getlist[len].ListNr = -1;
      }
      //     printf("l=%d\n", local);
      PROTECT(ans = getRFUoptions(getoptions, getlist, save, local, calling));
      n_protect++;
    }

    if (options != R_NilValue && !local && parallel())
      ERR0("'RFoptions(...)' may be used only outside any parallel code");

    for(i = 0; options != R_NilValue; i++, options = CDR(options)) {
      //      printf("set opt i=%d\n", i);
      el = CAR(options);
      if (!isNull(TAG(options))) {
	name = (char*) CHAR(PRINTNAME(TAG(options)));
	//	printf("AAxx '%s' %d %d local=%d\n", name, isList, i, local);
	setRFUoptions(el, name, isList, getlist, warn_unknown_option, local,
		      calling);
	//printf("done\n");
      }
    }
   FREE(getlist);
   }

  for (i=0; i<NList; i++) {
    if (setoption_fct_list[i][i] == NULL) {
      BUG;
    } else if (finaloption_fct_list[i] != NULL) {
      finaloption_fct_list[i](local);
    }
  }

  if (n_protect > 0) UNPROTECT(n_protect);

  if (!local) OPTIONS.basic.asList = true; // OK
 
  return(ans);
} 


void attachSetNGet(char * callingName, char *pkgname, setoptions_fctn set,
		   getoptions_fctn get) {
  int calling = 0;
  for ( ; calling<NList; calling++)
    if (STRCMP(callingName, pkgnames[calling]) == 0) break;
  if (calling >= NList) {BUG;}  
  for (int ListNr=0; ListNr<calling; ListNr++) {
    if (!STRCMP(pkgname, pkgnames[ListNr])) {
      if (setoption_fct_list[ListNr][calling] == NULL &&
	  getoption_fct_list[ListNr][calling] == NULL) {
	setoption_fct_list[ListNr][calling] = set;
	getoption_fct_list[ListNr][calling] = get;
	return;
      }
     { BUG;}
    }
  }
  { BUG;}
}


/*
Uint get_simd_info(char *pkgname) {
  int calling = 0;
  for ( ; calling<NList; calling++)
    if (STRCMP(callingName, pkgnames[calling]) == 0) break;
  if (calling >= NList) {ERR1("package '%.50s' does not exist\n", pkgname);}
  return simd_infos[NList];
}
*/

SEXP instruction_set(SEXP which,  SEXP pkgs, SEXP Uses) {
#if defined ARM32
#define N_LISTE 3
  const char* liste[N_LISTE] =
    {"CUDA", "SSE2/NEON", "SSSE3/NEON"};
 Uint
    bit[2][N_LISTE] =
    {{gpuMISS, sse2MISS, ssse3MISS},
     {gpuUSE, sse2USE, ssse3USE}};
#elif defined X86_64
#define N_LISTE 6
  const char* liste[N_LISTE] =
    {"CUDA", "SSE2", "SSSE3", "AVX", "AVX2", "AVX512F"};
 Uint
    bit[2][N_LISTE] =
    {{gpuMISS, sse2MISS, ssse3MISS, avxMISS, avx2MISS, avx512fMISS},
     {gpuUSE, sse2USE, ssse3USE, avxUSE, avx2USE, avx512fUSE}};
#else
#define N_LISTE 1
  const char* liste[N_LISTE] =
    {"CUDA"};
 Uint
    bit[2][N_LISTE] =
    {{gpuMISS},
     {gpuUSE}};
#endif
  int
    rowIdx[N_LISTE],
    colIdx[MAXNLIST],
    cols = length(pkgs),
    rows = length(which);
  Uint *simd_bit = bit[LOGICAL(Uses)[0]];
  if (rows == 0) rows = N_LISTE;
  if (cols == 0) cols = NList;
  if (cols > MAXNLIST) ERR0("duplicated package names or request on packages not supported by RandomFieldsUtils");
  if (rows > N_LISTE)
    ERR0("duplicated SIMD names or request on SIMD versions not supported by ");

  SEXP Ans, rownames, colnames;
  PROTECT(rownames = allocVector(STRSXP, rows));
  if (length(which) == 0)
    for (int i=0; i<rows; i++)
      SET_STRING_ELT(rownames, i, mkChar(liste[i]));
  else
    for (int i=0; i<rows; i++)
      SET_STRING_ELT(rownames, i, mkChar((char*) CHAR(STRING_ELT(which, i))));
  for (int i=0; i<rows; i++)
    rowIdx[i] = Match((char*) CHAR(STRING_ELT(rownames, i)), liste, N_LISTE);

  PROTECT(colnames = allocVector(STRSXP, cols));
  if (length(pkgs) == 0)
    for (int i=0; i<cols; i++)
      SET_STRING_ELT(colnames, i, mkChar(pkgnames[i]));
    else 
      for (int i=0; i<cols; i++)
	SET_STRING_ELT(colnames, i, mkChar((char*) CHAR(STRING_ELT(pkgs, i))));
  for (int i=0; i<cols; i++)
    colIdx[i] = Match((char*) CHAR(STRING_ELT(colnames, i)), pkgnames, NList);

  if (cols == 1) {
    PROTECT(Ans = allocVector(LGLSXP, rows));
    setAttrib(Ans, R_NamesSymbol, rownames);    
  } else {
    PROTECT(Ans = allocMatrix(LGLSXP, rows, cols));
    SEXP nameMatrix;
    PROTECT(nameMatrix = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(nameMatrix, 0, rownames);
    SET_VECTOR_ELT(nameMatrix, 1, colnames);
    setAttrib(Ans, R_DimNamesSymbol, nameMatrix);    
    UNPROTECT(1);
  }

  int *ans = (int*) LOGICAL(Ans);
  for (int p=0; p<cols; p++) {   
    if (colIdx[p] >= 0) {
      if (colIdx[p] >= NList) BUG;
      Uint simd_info = simd_infos[colIdx[p]];
      for (int i=0; i<rows; i++, ans++) {
	*ans = rowIdx[i] < 0 ? NA_INTEGER
			   : ((simd_info & 1U << simd_bit[rowIdx[i]]) > 0);
      }
    } else for (int i=0; i<rows; i++, ans++) *ans = NA_INTEGER;
  }
  UNPROTECT(3);
  return Ans;  
}

void attachRFUoptions(char *pkgname,
		      const char **PKGprefixlist, int N, 
		      const char ***PKGall, int *PKGallN,
		      setoptions_fctn set, getoptions_fctn get,
		      finalsetoptions_fctn final,
		      deleteoptions_fctn del,
		      setoptions_fctn setRFU, getoptions_fctn getRFU,
		      int pl_offset, bool basicopt,
		      install_modes gpu_needs, // neu
		      Uint simd_info, // neu
		      int version, // neu -- not used yet
		      int RFUversion,
		      int mem_is_aligned
		     ) {
  if (RFUversion != RFU_VERSION) {
    //printf("RFU=%d rfu=%d\n",  RFU_VERSION, RFUversion);
    if (RFUversion > RFU_VERSION) {
      ERR1("An obsolete version of RandomFieldsUtils has been installed in meanwhile that is incompatible the compiled version of '%.50s'.", pkgname);
    } else {
      ERR2("Package '%.50s' has been compiled against an obsolete version of RandomFieldsUtils. Please recompile '%.50s'.", pkgname, pkgname);
    }
  }
  installNrun_options *ip = &(OPTIONS.installNrun);

  
  if ((usr_bool) mem_is_aligned != ip->mem_is_aligned) {
#if defined SCHLATHERS_MACHINE
    PRINTF("mem alignment: %s=%d  RFU=%d\n",  pkgname, mem_is_aligned, ip->mem_is_aligned);
#else
    if (mem_is_aligned != Nan || ip->mem_is_aligned!=True)
      WARN2("'%.50s' is compiled with an alignment assumption different from the package 'RandomFieldsUtils'. See MEM_IS_ALIGNED and mem_is_aligned in ?RFoptions.\n  Recompile with 'RandomFieldsUtils::RFoptions(install.control=list(MEM_IS_ALIGNED=%.10s))'.", pkgname,
	    ((usr_bool) mem_is_aligned == Nan && ip->mem_is_aligned==True) ||
	    (usr_bool) mem_is_aligned==True ?
	    "TRUE" : "FALSE" // OK
	    );
#endif
  }
  for (int ListNr=0; ListNr<NList; ListNr++) {    
    if (AllprefixN[ListNr] == N && 
	STRCMP(Allprefix[ListNr][0], PKGprefixlist[0]) == 0) {
      if (PL > 0) {
	PRINTF("options starting with prefix '%s' have been already attached (%s %1.1f).",
	       PKGprefixlist[0], pkgname, version / 10.0);
      }
      return;    
    }
  }
  if (basicopt) option_class_list[noption_class_list++] = PKGprefixlist[0];
  if (NList >= MAXNLIST) BUG;
  strcopyN(pkgnames[NList], pkgname, PKGNAMELENGTH);
  Allprefix[NList] = PKGprefixlist;
  AllprefixN[NList] = N;
  Allall[NList] = PKGall;
  AllallN[NList] = PKGallN;
  Allversion[NList] = version;
  
  setoption_fct_list[NList][NList] = set;
  getoption_fct_list[NList][NList] = get;
  finaloption_fct_list[NList] = final;
  deloption_fct_list[NList] = del;

  // printf("simd_needs %d \n", simd_needs);
  if ((simd_info & 1<<avx512fMISS) && avx512fAvail)
    min_simd_needs[NList] = Iavx512f;
  else if ((simd_info & 1<<avx2MISS) && avx2Avail) min_simd_needs[NList]=Iavx2;
  else if ((simd_info & 1<<avxMISS) && avxAvail) min_simd_needs[NList] = Iavx;
  else if ((simd_info & 1<<ssse3MISS) && ssse3Avail)
    min_simd_needs[NList] = Issse3;
  else if ((simd_info & 1<<sse2MISS) && sse2Avail) min_simd_needs[NList]=Isse2;
  else min_simd_needs[NList] = Inone;

  //  printf("min =%d %d %d \n",  min_simd_needs[NList], Iavx512f, Iavx2);
  min_gpu_needs[NList] = gpu_needs;
  simd_infos[NList] = simd_info;
  
  assert(ip->install != Inone || !ip->installPackages);  
  if (ip->install != Inone)
    ip->installPackages |= min_simd_needs[NList] > Inone;
  // if only gpu_needs, the system has been recompiled already,
  // and probably no graphic card has been found

  //  printf("%s %d %d\n", pkgnames[NList],
  //	 ip->installPackages, min_simd_needs[NList]);
  
  NList++;
  PLoffset = pl_offset;
  basic_options *gp =  &(OPTIONS.basic);
  PL = gp->Cprintlevel = gp->Rprintlevel + PLoffset;  
  if (setRFU != NULL) {
    attachSetNGet(pkgname, (char *) "RandomFieldsUtils",
		  setRFU, getRFU);
  }
}



void detachRFUoptions(const char **PKGprefixlist, int N) {
  int ListNr;
  for (ListNr=0; ListNr<NList; ListNr++) {    
    if (AllprefixN[ListNr] == N && 
	STRCMP(Allprefix[ListNr][0], PKGprefixlist[0]) == 0) break;
  }  
  if (ListNr >= NList) {
    ERR1("options starting with prefix '%.50s' have been already detached.",
	PKGprefixlist[0]);
  }

   if (deloption_fct_list[ListNr] != NULL) deloption_fct_list[ListNr](false);
 
  
  int i;
  for (i=0; i<noption_class_list ; i++)
    if (STRCMP(option_class_list[i], PKGprefixlist[0]) == 0) break;
  for (i++ ; i < noption_class_list; i++)
    option_class_list[i - 1] = option_class_list[i];

  for (int call=ListNr+1; call<NList; call++) {    
    for (int Lnr =0; Lnr<ListNr; Lnr++) {
      setoption_fct_list[Lnr][call - 1] = setoption_fct_list[Lnr][call];
      getoption_fct_list[Lnr][call - 1] = getoption_fct_list[Lnr][call];
    }
    for(int Lnr=ListNr+1; Lnr<=call; Lnr++) {
      setoption_fct_list[Lnr -1][call -1] = setoption_fct_list[Lnr][call];
      getoption_fct_list[Lnr -1][call -1] = getoption_fct_list[Lnr][call];
    }
  }

  for (ListNr++; ListNr<NList; ListNr++) {    
    Allprefix[ListNr - 1] = Allprefix[ListNr];
    AllprefixN[ListNr - 1] =  AllprefixN[ListNr];
    Allall[ListNr - 1] = Allall[ListNr];
    AllallN[ListNr - 1] = AllallN[ListNr];
    
    finaloption_fct_list[ListNr - 1] = finaloption_fct_list[ListNr];
    deloption_fct_list[ListNr - 1] = deloption_fct_list[ListNr];
   
    min_simd_needs[ListNr - 1] = min_simd_needs[ListNr];
    min_gpu_needs[ListNr - 1] = min_gpu_needs[ListNr];
    simd_infos[ListNr - 1] = simd_infos[ListNr];
  }

  NList--;
  if (NList <= 1) PLoffset = 0;
}


void resetInstalled() {
  for (int ListNr=0; ListNr<NList; ListNr++)
    installed[ListNr] = min_simd_needs[ListNr] == Inone;
  // printf("installed %d %d\n", installed[0], installed[1]);
}



SEXP getPackagesToBeInstalled(SEXP Force) {
  installNrun_options *ip = &(OPTIONS.installNrun);
  ip->installPackages = false;
  int zaehler =0;
  bool force = ip->mem_is_aligned != MEMisALIGNED ? true : LOGICAL(Force)[0];
  //printf("get %d\n", NList);
  for (int ListNr=0; ListNr<NList; ListNr++)
    // again, do not consider gpu_needs since likely graphic card is missing
    zaehler += force || (!installed[ListNr] && min_simd_needs[ListNr] != Inone);
  //printf("zaehler %d %d\n", zaehler, force);
  if (zaehler == 0) return R_NilValue; // can happen due to user's call
  // and by call from .onLoad
  SEXP str;
  //
  PROTECT(str = allocVector(STRSXP, zaehler));
  zaehler = 0;
  for (int ListNr=0; ListNr < NList; ListNr++) {
    //    printf("%d, %d %s %d\n", zaehler, ListNr, pkgnames[ListNr], min_simd_needs[ListNr]);
    if (force || (!installed[ListNr] && min_simd_needs[ListNr] != Inone)) {
      //      printf("add %d, %s\n", zaehler, pkgnames[zaehler]);
      SET_STRING_ELT(str, zaehler++, mkChar(pkgnames[ListNr]));
      installed[ListNr] = true;
    }
  }
  UNPROTECT(1);
  return str;
}

SEXP isGPUavailable() {
  SEXP ans =  PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] =
#if defined USEGPU
    true;
#else
  GPU_NEEDS==Igpu;
  //  printf("GPU_NEEDS = %d %d\n", GPU_NEEDS,Igpu);
#endif
  UNPROTECT(1);
  return ans;
}


SEXP isNEONavailable() {
  SEXP ans =  PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] =
#if defined NEON
    true;
#elif defined ARM32
  false;
#else
  NA_INTEGER;
#endif
  UNPROTECT(1);
  return ans;
}

SEXP isX86_64() {
    SEXP ans =  PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] =
#if defined X86_64
    true;
#else
  false;
#endif
  UNPROTECT(1);
  return ans;
}

void setCPUs(int *n) {
  //  printf("set cpus = %d %d %s\n", *n, min_simd_needs[0], pkgnames[0]);
    basic_options *gp = &(OPTIONS.basic);    
#ifdef DO_PARALLEL
    if (min_simd_needs[0] == Inone) {
      gp->cores = GreaterZero(*n);
  }
#else
  gp->cores = 1;
#endif
}

//## 5 8, 0 8

void recompilationNeeded(int *n) {
  *n = (int) false;  
  if (INSTALL_DEFAULT == Inone) return;  
  for (int ListNr=0; ListNr<NList; ListNr++) {
    //    printf("%d %d %s\n", ListNr, min_simd_needs[ListNr], pkgnames[ListNr]);
    *n |= min_simd_needs[ListNr] != Inone;
  }
#if defined USERASKED
  if (*n) *n = -1 - USERASKED;
#endif
}


SEXP SIMDmessages(SEXP pkgs) {
#if defined X86_64 || defined NEON
#if defined __APPLE__
  const char OMP[80] = "";
#else
 const char OMP[80] =
    " -Xpreprocessor -fopenmp -pthread' LIB_FLAGS='-lgomp -pthread";
 if (HAS_PARALLEL) STRCPY((char *) OMP, "");
#endif
   
  if (length(pkgs) == 0) {
    PRINTF("If re-compilation since it does not work, consider one of the following options:\n");
    install_modes required = Inone;
    //    printf("AV ok\n");
    for (int ListNr=0; ListNr<NList; ListNr++)
      if (min_simd_needs[ListNr] > required) required = min_simd_needs[ListNr];
    //  printf("AV ok x\n");
#if defined MSDOS_WINDOWS
    if (required > Inone) {
      PRINTF("\n\nBy default the packages are compiled without flag '-mavx' under your OS.\nIf you are sure that AVX2 is available, consider adding the flag '-march=native'\nto 'PKG_CXXFLAGS' in the file src/Makefile.win and then recompile\n'");
      if (required >= Iavx2)
	PRINTF("Or: try adding flag '-mavx' or '-mavx2' to 'PKG_CXXFLAGS'\n");
    }
    if (!HAS_PARALLEL)
      PRINTF("For OMP alone, try adding the flags -Xpreprocessor -fopenmp -pthread to PKG_LIBS and PKG_CXXFLAGS");
#elif defined ARM32
    if (required > Inone) {
      PRINTF("\n\n   install.packages(<package>, configure.args=\"CROSS='FALSE'%s\")\n   install.packages(<package>, configure.args=\"CROSS='FALSE' USE_GPU='yes'%s\")", OMP, OMP); // OK
    }
 #else    
    if (required > Inone) {
      PRINTF("\n\n   install.packages(<package>, configure.args=\"CXX_FLAGS='-march=native%s'\")\n\n   install.packages(<package>, configure.args=\"CXX_FLAGS='-march=native%s' USE_GPU='yes'\")",OMP, OMP);
      if (required > Iavx) {
	PRINTF("\n\n   install.packages(<package>, configure.args=\"CXX_FLAGS='-mavx%s'\")", OMP);
	if (required > Iavx2) 
	  PRINTF("\\n   install.packages(<package>, configure.args=\"CXX_FLAGS='-mavx2%s'\")", OMP);
      }
    }
    /*   HINT && MISS_ANY_SIMD ? "\n\nAlternatively,\n    install.packages(\""#PKG"\", configure.args=\"USE_AVX='yes'\") \nOr, consider installing '"#PKG"'\nfrom https://github.com/schlather/"#PKG", i.e.,\n   install.packages(\"devtools\")\n   library(devtools)\n   devtools::install_github(\"schlather/"#PKG"/pkg\")" : "", */
#endif
 #if ! defined __APPLE__ 
   if (!HAS_PARALLEL)
      PRINTF("\n\nFor OMP alone try\n   install.packages(<package>, configure.args=\"CXX_FLAGS='%s'\")", OMP);
#endif   
  } else {
    if (!STRCMP("OMP", CHAR(STRING_ELT(pkgs, 0))))
      return ScalarString(mkChar(OMP));
    bool all = !STRCMP("all", CHAR(STRING_ELT(pkgs, 0)));
    int len =  all ? NList : length(pkgs);
    for (int i=0; i<len; i++) {
      for (int ListNr=0; ListNr<NList; ListNr++) {
 	//	PRINTF("i=%d %d\n",i, ListNr);
	//	PRINTF("%s\n", pkgnames[ListNr]);
	//	PRINTF("%s\n", CHAR(STRING_ELT(pkgs, i)));
	if ((all && i==ListNr) ||
	    (!all && !STRCMP(pkgnames[ListNr], CHAR(STRING_ELT(pkgs, i))))) {	  
	  Uint simd_info = simd_infos[ListNr];
	  PRINTF("%s ", pkgnames[ListNr]);		 
	  PRINTF("%s ", simd_info & 1<<anyrelevantUSE ?  "sees" : "does not see any of");
	  if (simd_info & 1<<gpuUSE) PRINTF("GPU, ");
	  if (simd_info & 1<<avx512fUSE) PRINTF("AVX512F, ");
	  if (simd_info & 1<<avx2USE) PRINTF("AVX2, ");
	  if (simd_info & 1<<avxUSE) PRINTF("AVX, ");
	  if (simd_info & 1<<ssse3USE) PRINTF("SSSE3, ");
	  if (simd_info & 1<<sse2USE) PRINTF("SSE2, ");
	  if (HAS_PARALLEL) PRINTF("OMP, ");
 	  if (simd_info & 1<<USEnMISS) {
	    PRINTF(" but not ");
	    if (simd_info & 1<<gpuMISS) PRINTF("GPU, ");
	    if ((simd_info & 1<<avx512fMISS) > 0 && avx2Avail)
	      PRINTF("AVX512F, ");
	    if ((simd_info & 1<<avx2MISS) > 0 && avx2Avail) PRINTF("AVX2, ");
	    if ((simd_info & 1<<avxMISS) > 0 && avxAvail) PRINTF("AVX,");
	    if ((simd_info & 1<<ssse3MISS) > 0 && ssse3Avail) PRINTF("SSSE3, ");
	    if ((simd_info & 1<<sse2MISS) > 0 && sse2Avail) PRINTF("SSE2, ");
	    if (!HAS_PARALLEL) PRINTF("OMP.");
	  }
	  PRINTF("\n");
	}
      }
    }
  }
  PRINTF("\n\nOr call 'RFoptions(install=\"no\")' after loading to avoid being asked again.\n");

#else // neither ARM nor X86_64
   const char OMP[80] = "";
#endif
   
  return ScalarString(mkChar(OMP));
}



SEXP Update_utilsoption() {
  update_utilsoption();
  return R_NilValue;
}

#endif
