
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 -- 2017 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#include "RandomFieldsUtils.h"
#include "General_utils.h"
#include "zzz_RandomFieldsUtils.h"
#include "xport_import.h"
#include "kleinkram.h"

#ifdef _WIN32
#include <intrin.h>							
#endif


AVAILABLE


struct getlist_type{
  int ListNr, i;
};


void setpDef(int VARIABLE_IS_NOT_USED  i, 
	     int VARIABLE_IS_NOT_USED  j, 
	     SEXP VARIABLE_IS_NOT_USED  el,
	     char VARIABLE_IS_NOT_USED  name[LEN_OPTIONNAME], 
	     bool VARIABLE_IS_NOT_USED  isList,
	     bool VARIABLE_IS_NOT_USED local) {
  BUG;
}
void getpDef(SEXP VARIABLE_IS_NOT_USED  sublist, int VARIABLE_IS_NOT_USED i,
	     bool VARIABLE_IS_NOT_USED local) {
  BUG;
}

void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);

#define MAXNLIST 7
#define PKGNAMELENGTH 20
int NList = 0; // originally 1
int noption_class_list = 0,
  AllprefixN[MAXNLIST] = {prefixN, 0, 0, 0, 0, 0, 0},
  *AllallN[MAXNLIST] = {allN, NULL, NULL, NULL, NULL, NULL, NULL};
const char  *option_class_list[MAXNLIST] =
  {prefixlist[1], NULL, NULL, NULL, NULL, NULL},
  **Allprefix[MAXNLIST] = {prefixlist, NULL, NULL, NULL, NULL, NULL, NULL},
  ***Allall[MAXNLIST] = { all, NULL, NULL, NULL, NULL, NULL, NULL};
char pkgnames[MAXNLIST][PKGNAMELENGTH+1] = {"", "", "", "", "", "", ""};
setoptions_fctn setoption_fct_list[MAXNLIST] = 
  {setoptions, setpDef, setpDef, setpDef, setpDef, setpDef, setpDef};
getoptions_fctn getoption_fct_list[MAXNLIST] = 
  {getoptions, getpDef, getpDef, getpDef, getpDef, getpDef, getpDef};
finalsetoptions_fctn finaloption_fct_list[MAXNLIST] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL };
deleteoptions_fctn deloption_fct_list[MAXNLIST] =
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL};
bool installed [MAXNLIST] = { false, false, false, false, false, false, false};
install_modes min_avx_needs[MAXNLIST] = // currently only used if == none or not
  {Inone, Inone, Inone, Inone, Inone, Inone, Inone},
  min_gpu_needs[MAXNLIST] = // currently not used
  {Inone, Inone, Inone, Inone, Inone, Inone, Inone};
Uint avx_infos [MAXNLIST] = {0, 0, 0, 0, 0, 0, 0};
  
/* precide information of min_install is not used yet; just whether it is 
   'none' or not */

void hintVariable(char *name, int warn_unknown_option) {
  static bool printing = true; // da nur Mutterprozess schreiben darf
  if (warn_unknown_option > 0 && OPTIONS.basic.Rprintlevel > 0) {
    PRINTF("'%s' is considered as a variable (not as an option).\n", name);
    if (printing && OPTIONS.basic.helpinfo && !parallel()) {
      PRINTF("[These hints can be turned off by 'RFoptions(%s=-%d)'.]\n",
	     basic[BASIC_WARN_OPTION], MAX(1, warn_unknown_option));
      printing = false;
    }
  }
}

void setparameter(SEXP el, char *prefix, char *mainname, bool isList,
		  getlist_type *getlist, int warn_unknown_option, bool local) {
  int 
    j = NOMATCHING,
    i = NOMATCHING,
    ListNr = NOMATCHING;
  char name[LEN_OPTIONNAME];

  SPRINTF(name, "%.50s%.50s%.50s", prefix, STRLEN(prefix)==0 ? "" : ".",
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
	if (ii >= 0 && STRCMP(prefix, Allprefix[k][ii])==0) {
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
      case WARN_UNKNOWN_OPTION_NONE : case -WARN_UNKNOWN_OPTION_NONE1 :
	return;
      case WARN_UNKNOWN_OPTION_CAPITAL : case -WARN_UNKNOWN_OPTION_CAPITAL :
	if (name[0] >= 'A' && name[1] <= 'Z') {
	  hintVariable(name, warn_unknown_option);
	  return;
	}
	ERR1("Unknown option '%.50s'.", name); 
      case WARN_UNKNOWN_OPTION_SINGLE : case -WARN_UNKNOWN_OPTION_SINGLE :
	if (STRLEN(name) == 1 && name[0] >= 'A' && name[1] <= 'Z') {
	  hintVariable(name, warn_unknown_option);
	  return;
	}
	ERR1("Unknown option '%.50s'.", name); 
      case WARN_UNKNOWN_OPTION_ALL :
      default :
	ERR1("Unknown option '%.50s'.", name);
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
	  //	 jj, jj < 0 ? "none" : Allall[k][ii][jj]);
	  j = MULTIPLEMATCHING;
	  if (jj >= 0 && STRCMP(mainname, Allall[k][ii][jj])==0) {
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
  
  setoption_fct_list[ListNr](i, j, el, name, isList, local); 
}


SEXP getRFoptions(int ListNr, int i, bool local) {
  SEXP sublist, subnames;
  int elmts = AllallN[ListNr][i];
  PROTECT(sublist = allocVector(VECSXP, elmts));
  PROTECT(subnames = allocVector(STRSXP, elmts));
  for (int k=0; k<elmts; k++) {
    //   printf("getopt %d k=%d elmnts=%d %d %.50s\n", i, k, elmts, ListNr, Allall[ListNr][i][k]);
    SET_STRING_ELT(subnames, k, mkChar(Allall[ListNr][i][k])); 
  }
  getoption_fct_list[ListNr](sublist, i, local);
  setAttrib(sublist, R_NamesSymbol, subnames);
  UNPROTECT(2);
  return sublist;
}

SEXP getRFoptions(bool local) {
  SEXP list, names;
 
  int  PreFixN, totalN, i, ListNr,
    itot = 0;
  for (totalN=ListNr=0; ListNr<NList; ListNr++) {
    PreFixN = AllprefixN[ListNr];
    for (i=0; i<PreFixN; i++) {
      totalN += STRCMP(Allprefix[ListNr][i], OBSOLETENAME) != 0;
    }
  }

  PROTECT(list = allocVector(VECSXP, totalN));
  PROTECT(names = allocVector(STRSXP, totalN));

  for (ListNr =0; ListNr<NList; ListNr++) {
    //    printf("ListNr %d\n", ListNr);
    PreFixN = AllprefixN[ListNr];    
    for (i=0; i<PreFixN; i++) {
      if (STRCMP(Allprefix[ListNr][i], OBSOLETENAME) == 0) continue;
      SET_VECTOR_ELT(list, itot, getRFoptions(ListNr, i, local));
      SET_STRING_ELT(names, itot, mkChar(Allprefix[ListNr][i]));
      itot ++;
    } 
  }

  setAttrib(list, R_NamesSymbol, names);   
  UNPROTECT(2);

  return list;
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
  if (ListNr >= NList) ERR("unknown value for 'getoptions_'");
  if (getlist != NULL) {
    getlist[t].ListNr = ListNr;
    getlist[t].i = i;
  }
  *Nr = ListNr;
  *idx = i;
}

  
SEXP getRFoptions(SEXP which, getlist_type *getlist, bool save, bool local) {
  
  int ListNr, idx,
    actual_nbasic = noption_class_list * save,
    totalN = length(which) + actual_nbasic;  

  if (totalN == 0) return R_NilValue;
  if (totalN == 1) {
    getListNr(save, 0, actual_nbasic, which, getlist, &ListNr, &idx);
    return getRFoptions(ListNr, idx, local);
  }

  SEXP list, names;
  PROTECT(list = allocVector(VECSXP, totalN));
  PROTECT(names = allocVector(STRSXP, totalN));
  for (int t=0; t<totalN; t++) {
    getListNr(save, t, actual_nbasic,  which, getlist, &ListNr, &idx);
    SET_VECTOR_ELT(list, t, getRFoptions(ListNr, idx, local));
    SET_STRING_ELT(names, t, mkChar(Allprefix[ListNr][idx]));
  }
  setAttrib(list, R_NamesSymbol, names);     
  UNPROTECT(2);
  return list;
}


void splitAndSet(SEXP el, char *name, bool isList, getlist_type *getlist,
		 int warn_unknown_option, bool local) {
  int i, len;
  char PreFix[LEN_OPTIONNAME / 2], mainname[LEN_OPTIONNAME / 2];   
  //  printf("splitandset\n");
  len = STRLEN(name);
  for (i=0; i < len && name[i]!='.'; i++);
  if (i==0) { ERR1("argument '%.50s' not valid\n", name); }
  if (i==len) {
    STRCPY(PreFix, "");
    strcopyN(mainname, name, LEN_OPTIONNAME / 2);
  } else {
    strcopyN(PreFix, name, MIN(i + 1, LEN_OPTIONNAME / 2));
    strcopyN(mainname, name+i+1, MIN(STRLEN(name) - i, LEN_OPTIONNAME / 2) );
  }

  // 
  //  printf("i=%d %d %.50s %.50s\n", i, len, prefix, mainname);
  setparameter(el, PreFix, mainname, isList && OPTIONS.basic.asList, getlist,
	       warn_unknown_option, local);
  //   printf("ende\n");
}


SEXP RFoptions(SEXP options) {
  int i, j, lenlist, lensub;
  SEXP el, list, sublist, names, subnames,
     ans = R_NilValue;
  char *name, *pref;
  bool isList = false;
  int
    warn_unknown_option = OPTIONS.basic.warn_unknown_option,
    n_protect = 0;
  bool local = false;

    /* 
     In case of strange values of a parameter, undelete
     the comment for PRINTF
  */

  //  printf("CORES RFOP = %d %d\n", CORES, OPTIONS.basic.cores);
  
  options = CDR(options); /* skip 'name' */

  // printf(" %d %d\n", OPTIONS.basic.warn_unknown_option, options == R_NilValue);
  name = (char*) "";
  if (options != R_NilValue) {
    if (!isNull(TAG(options))) name = (char*) CHAR(PRINTNAME(TAG(options)));
    if (STRCMP(name, "local_")==0) {
      el = CAR(options);
      local = (bool) INT;
      options = CDR(options); /* skip 'name' */
     }
  }
  
  if (options == R_NilValue || STRCMP(name, "") ==0)
    return getRFoptions(local);
 
  if (!isNull(TAG(options))) name = (char*) CHAR(PRINTNAME(TAG(options)));
  if (STRCMP(name, "warnUnknown_")==0) {
    el = CAR(options);
    warn_unknown_option = INT;
    options = CDR(options); /* skip 'name' */
   }

  if (!isNull(TAG(options))) name = (char*) CHAR(PRINTNAME(TAG(options)));
  if ((isList = STRCMP(name, "list_")==0)) {
    if (local) ERR("'list_' can be used only globally.");
    list = CAR(options);
    if (TYPEOF(list) != VECSXP)
      ERR1("'list_' needs as argument the output of '%.50s'", RFOPTIONS);
    PROTECT(names = getAttrib(list, R_NamesSymbol));
    lenlist = length(list);

    if (lenlist > 0 && !local && parallel())
      ERR("Global 'RFoptions' such as 'cores', 'seed' and 'printlevel', may be set only outside any parallel code. See '?RandomFieldsUtils::RFoptions' for the complete list of global 'RFoptions'");

   for (i=0; i<lenlist; i++) {
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
		       local);
	}
	UNPROTECT(1);
      } else {   
	splitAndSet(sublist, pref, isList, NULL, warn_unknown_option, local);
      }
    }
    UNPROTECT(1);
  } else {    
    getlist_type *getlist = NULL;
    bool save;
    if ((save = STRCMP(name, "saveoptions_")==0 ) ||
	STRCMP(name, "getoptions_")==0) {
      SEXP getoptions = CAR(options);
      options = CDR(options);
      if (options != R_NilValue) {
	//	printf("hier\n");
	int len = length(getoptions) + noption_class_list * save;
	getlist = (getlist_type *) MALLOC(sizeof(getlist_type) * (len + 1));
	getlist[len].ListNr = -1;
      }
      //     printf("l=%d\n", local);
      PROTECT(ans = getRFoptions(getoptions, getlist, save, local));
      n_protect++;
    }

    if (options != R_NilValue && !local && parallel())
      ERR("'RFoptions(...)' may be used only outside any parallel code");

    for(i = 0; options != R_NilValue; i++, options = CDR(options)) {
      //printf("set opt i=%d\n", i);
      el = CAR(options);
      if (!isNull(TAG(options))) {
	name = (char*) CHAR(PRINTNAME(TAG(options)));
	//printf("AAxx '%s' %d %d local=%d\n", name, isList, i, local);
	splitAndSet(el, name, isList, getlist, warn_unknown_option, local);
      }
    }
    FREE(getlist);
  }


  for (i=0; i<NList; i++)
    if (finaloption_fct_list[i] != NULL) {      finaloption_fct_list[i]();
    }

  if (n_protect > 0) UNPROTECT(n_protect);

  if (!local) OPTIONS.basic.asList = true; // OK
 
  return(ans);
} 
 

int PLoffset = 0;
void attachRFoptions(char *pkgname,
		     const char **PKGprefixlist, int N, 
		     const char ***PKGall, int *PKGallN,
		     setoptions_fctn set, finalsetoptions_fctn final,
		     getoptions_fctn get, deleteoptions_fctn del,
		     int pl_offset, bool basicopt,
		     install_modes avx_needs, install_modes gpu_needs,
		     Uint avx_info) {
  for (int ListNr=0; ListNr<NList; ListNr++) {    
    if (AllprefixN[ListNr] == N && 
	STRCMP(Allprefix[ListNr][0], PKGprefixlist[0]) == 0) {
      if (PL > 0) {
	PRINTF("options starting with prefix '%s' have been already attached.",
		PKGprefixlist[0]);
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
  setoption_fct_list[NList] = set;
  finaloption_fct_list[NList] = final;
  getoption_fct_list[NList] = get;
  deloption_fct_list[NList] = del;
  // printf("avx_needs %d \n", avx_needs);
  switch(avx_needs) {
  case Iavx2 : min_avx_needs[NList] = Iavx2; if (avx2) break;
  case Iavx : min_avx_needs[NList] = Iavx; if (avx) break;
  case Issse3 :min_avx_needs[NList] = Issse3; if (ssse3) break;
  case Isse3 :min_avx_needs[NList] = Isse3; if ( sse3) break;
  case Isse2 :min_avx_needs[NList] = Isse2; if (sse2) break;
  case Isse : min_avx_needs[NList] = Isse; if ( sse) break;
  case Inone : min_avx_needs[NList] = Inone; break;
  default: BUG;
  }
  min_gpu_needs[NList] = gpu_needs;
  avx_infos[NList] = avx_info;
  
  basic_options *gp = &(OPTIONS.basic);
  assert((OPTIONS.basic.install != Inone) || !gp->installPackages);  
  if (OPTIONS.basic.install != Inone)
    gp->installPackages |= min_avx_needs[NList] > Inone;
  // if only gpu_needs, the system has been recompiled already,
  // and probably no graphic card has been found

  //  printf("%s %d %d\n", pkgnames[NList],
  //	 gp->installPackages, min_avx_needs[NList]);
  
  NList++;
  PLoffset = pl_offset;
  PL = gp->Cprintlevel = gp->Rprintlevel + PLoffset;  
  CORES = gp->cores;
}

void detachRFoptions(const char **PKGprefixlist, int N) {
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
  
  for (ListNr++; ListNr<NList; ListNr++) {    
    Allprefix[ListNr - 1] = Allprefix[ListNr];
    AllprefixN[ListNr - 1] =  AllprefixN[ListNr];
    Allall[ListNr - 1] = Allall[ListNr];
    AllallN[ListNr - 1] = AllallN[ListNr];
    setoption_fct_list[ListNr - 1] = setoption_fct_list[ListNr];
    finaloption_fct_list[ListNr - 1] = finaloption_fct_list[ListNr];
    getoption_fct_list[ListNr - 1] = getoption_fct_list[ListNr];
    min_avx_needs[ListNr - 1] = min_avx_needs[ListNr];
    min_gpu_needs[ListNr - 1] = min_gpu_needs[ListNr];
    avx_infos[ListNr - 1] = avx_infos[ListNr];
  }

  NList--;
  if (NList <= 1) PLoffset = 0;
}

void getUtilsParam(utilsoption_type **global) { 
  // printf("GLO %ld\n", &OPTIONS);
  *global = &OPTIONS; 
}


void resetInstalled() {
  for (int ListNr=0; ListNr<NList; ListNr++)
    installed[ListNr] = min_avx_needs[ListNr] == Inone;
  // printf("installed %d %d\n", installed[0], installed[1]);
}



SEXP getPackagesToBeInstalled(SEXP Force) {
  basic_options *gp = &(OPTIONS.basic);
  gp->installPackages = false;
  int zaehler =0;
  bool force = LOGICAL(Force)[0];
  //printf("get %d\n", NList);
  for (int ListNr=0; ListNr<NList; ListNr++)
    // again, do not consider gpu_needs since likely graphic card is missing
    zaehler += force || (!installed[ListNr] && min_avx_needs[ListNr] != Inone);
  //printf("zaehler %d %d\n", zaehler, force);
  if (zaehler == 0) return R_NilValue; // can happen due to user's call
  // and by call from .onLoad
  SEXP str;
  //
  PROTECT(str = allocVector(STRSXP, zaehler));
  zaehler = 0;
  for (int ListNr=0; ListNr < NList; ListNr++) {
    //    printf("%d, %d %s %d\n", zaehler, ListNr, pkgnames[ListNr], min_avx_needs[ListNr]);
    if (force || (!installed[ListNr] && min_avx_needs[ListNr] != Inone)) {
      //      printf("add %d, %s\n", zaehler, pkgnames[zaehler]);
      SET_STRING_ELT(str, zaehler++, mkChar(pkgnames[ListNr]));
      installed[ListNr] = true;
    }
  }
  
  UNPROTECT(1);
  return str;
}

void setCPUs(int *n) {
  //  printf("set cpus = %d %d %s\n", *n, min_avx_needs[0], pkgnames[0]);
  basic_options *gp = &(OPTIONS.basic);
  if (min_avx_needs[0] == Inone) {
    CORES = gp->cores = *n;
    //   printf("CORES %d\n", CORES);
  }
}

//## 5 8, 0 8

void recompilationNeeded(int *n) {
  //  printf("set cpus = %d %d %s\n", *n, min_avx_needs[0], pkgnames[0]);
  *n = (int) false;
  for (int ListNr=0; ListNr<NList; ListNr++) {
    //    printf("%d %d %s\n", ListNr, min_avx_needs[ListNr], pkgnames[ListNr]);
    *n |= min_avx_needs[ListNr] != Inone;
  }
}


SEXP AVXmessages(SEXP pkgs) {
  const char OMP[80] =
    " -Xpreprocessor -fopenmp -pthread' LIB_FLAGS='-lgomp -pthread";
  if (HAS_PARALLEL) STRCPY((char *) OMP, "");
   
  if (length(pkgs) == 0) {
    install_modes required = Inone;
    //    printf("AV ok\n");
    for (int ListNr=0; ListNr<NList; ListNr++)
      if (min_avx_needs[ListNr] > required) required = min_avx_needs[ListNr];
    //  printf("AV ok x\n");
 #ifdef WIN32
    if (required > Inone) {
      PRINTF("\n\nBy default the packages are compiled without flag '-mavx' under your OS.\nIf you are sure that AVX2 is available, consider adding the flag '-march=native'\nto 'PKG_CXXFLAGS' in the file src/Makefile.win and then recompile\n'");
      if (required >= Iavx2)
	PRINTF("Or: try adding flag '-mavx' or '-mavx2' to 'PKG_CXXFLAGS'\n");
    }
    if (!HAS_PARALLEL)
      PRINTF("For OMP alone, try adding the flags -Xpreprocessor -fopenmp -pthread to PKG_LIBS and PKG_CXXFLAGS");
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
	  Uint avx_info = avx_infos[ListNr];
	  PRINTF("%s ", pkgnames[ListNr]);		 
	  PRINTF("%s ", avx_info & 1 ?  "sees" : "does not see any of");
	  if (avx_info & 1<<1) PRINTF("GPU, ");
	  if (avx_info & 1<<2) PRINTF("AVX2, ");
	  if (avx_info & 1<<3) PRINTF("AVX, ");
	  if (avx_info & 1<<4) PRINTF("SSSE3, ");
	  if (avx_info & 1<<5) PRINTF("SSE2, ");
	  if (avx_info & 1<<6) PRINTF("SSE, ");
	  if (HAS_PARALLEL) PRINTF("OMP, ");
 	  if (avx_info & 1<<10) {
	    PRINTF(" but not ");
	    if (avx_info & 1<<11) PRINTF("GPU, ");
	    if ((avx_info & 1<<12) > 0 && avx2) PRINTF("AVX2, ");
	    if ((avx_info & 1<<13) > 0 && avx) PRINTF("AVX,");
	    if ((avx_info & 1<<14) > 0 && ssse3) PRINTF("SSSE3, ");
	    if ((avx_info & 1<<15) > 0 && sse2) PRINTF("SSE2, ");
	    if ((avx_info & 1<<16) > 0 && sse) PRINTF("SSE, ");
	    if (!HAS_PARALLEL) PRINTF("OMP.");
	  }
	  PRINTF("\n");
	}
      }
    }
  }
  PRINTF("\n"); 
  
  return ScalarString(mkChar(OMP));
}


