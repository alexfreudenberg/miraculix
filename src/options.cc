/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2001 -- 2016 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Linpack.h>
#include <stdio.h>  
//#include <stdlib.h>
#include <unistd.h>
#include <string.h>

// ACHTUNG: Reihenfolge nicht aendern!
#include "Miraculix_aux.h"
#include "miraculix.h"
#include "options.h"
//#include "def.h"
#include "xport_import.h"
#include "AutoMiraculix.h"
#include "intrinsics.h"

 
CALL1(void, getErrorString, errorstring_type, errorstring)
CALL1(void, setErrorLoc, errorloc_type, errorloc)



const char * prefixlist[prefixN] = 
  {"relationshipmatrix"};


// IMPORTANT: all names of general must be at least 3 letters long !!!
const char *relationship[relationshipN] = 
  {"digits", "relship_method", "centered", "normalized", "per_snp",
   "returnsigma"};

  

int PL=C_PRINTLEVEL,
  CORES = INITCORES; // err
globalparam GLOBAL = {
relationship_START
};
utilsparam *GLOBAL_UTILS;


const char **all[prefixN] = {relationship};
int allN[prefixN] = {relationshipN};
 

void setparameter(int i, int j, SEXP el, char name[200], 
		  bool VARIABLE_IS_NOT_USED isList, int local) {
#ifdef DO_PARALLEL
  if (local != isGLOBAL) ERR1("Options specific to RandomFieldsUtils, here '%.50s', can be set only via 'RFoptions' outside any parallel code.", name);
#endif  
  globalparam *options = &GLOBAL; 
  switch(i) {
  case 0: {// relationship
    relationship_param *gp;
    gp = &(options->relationship);
    switch(j) {
    case 0: gp->digits = NUM; break;
    case 1: {
      int m = TYPEOF(el) != STRSXP ? POS0NUM
	: GetName(el, name, RELSHIP_METH_NAME, last_usr_meth + 1, gp->method);
#if not defined SSE2 and not defined AVX
      if (m == Shuffle || m == Hamming2)
	ERR("'Shuffle' and 'Hamming2' are not available on that machine. Consider recompilation of the package with -avx or -msse2 or similar.\n");
#endif
#if not defined SSSE3 
      if (m == Hamming3)
	ERR("'hamming3' is not available on that machine. Consider recompilation of the package with -msse3 or similar.\n");
#endif
      gp->method = m;
      break;
    }
    case 2:
      if (TYPEOF(el) == LGLSXP) {
	gp->normalized = (gp->centered = USRLOG);
	if (gp->centered == Nan) {
	  if (gp->pcentered == NULL) {
	    WARN1("'%.50s' set to TRUE", name); // OK
	    gp->normalized = (gp->centered = True);
	  }
	} else {
	  FREE(gp->pcentered);
	  gp->ncentered = 0;
	}
      } else {
	int len = length(el);
	gp->ncentered = len;
	FREE(gp->pcentered);
	gp->pcentered = (double*) MALLOC(len * sizeof(double));
	Real(el, name, gp->pcentered, len);
	gp->centered = Nan;
	gp->normalized = False;
      }
      break;      
    case 3:
      gp->normalized = LOGI;
      if (gp->normalized && gp->centered != True) { 
	warn("'normalized=TRUE' only allowed with 'centered=TRUE'.\n'normalized=FALSE' is kept"); // OK
	gp->normalized = false;
      }
      break;
    case 4: gp->per_snp = LOGI; break;
    case 5 : gp->returnsigma = LOGI; break;
    default: BUG; 
    }
  }
    break;
  default: BUG;
  }
}


#define PLoffset -10
  void finalparameter(int VARIABLE_IS_NOT_USED local) {
  PL = GLOBAL_UTILS->basic.Cprintlevel - PLoffset;
  CORES = GLOBAL_UTILS->basic.cores;
}


void getparameter(SEXP sublist, int i, int VARIABLE_IS_NOT_USED local) {
  int k;
#ifdef DO_PARALLEL
  //  if (local != isGLOBAL) ERR("Options specific to RandomFieldsUtils can be obtained only on a global level and outside any parallel code.");
#endif  
  globalparam *options = &GLOBAL; 
 switch(i) {
  case 0 : {
    k = 0;
    relationship_param *p = &(options->relationship);
    ADD(ScalarReal(p->digits));
    //    ADD(ScalarString(mkChar(RELSHIP_METH_NAME[p->method])));
    ADD(ScalarInteger(p->method));
    ADD(ExtendedBooleanUsr(p->centered));    
    ADD(ScalarLogical(p->normalized));
    ADD(ScalarLogical(p->per_snp));
    ADD(ScalarLogical(p->returnsigma));
  }
    break;
  default : BUG;
  }
}


void attachmiraculix() {
  includeXport();
  Ext_getUtilsParam(&GLOBAL_UTILS);
  GLOBAL_UTILS->solve.max_chol = 8192;
  GLOBAL_UTILS->solve.max_svd = 6555;  
/*
  spam.min.n = as.integer(400) # 400
  spam.tol = 1e-8
  spam.min.p = 0.8 # 0.8
  spam.n = as.integer(1:500)
  spam.factor = as.integer(4294967)# prime, so that n * factor is still integer
  silent = T RUE
*/

  finalparameter(isGLOBAL);
  Ext_attachRFoptions(prefixlist, prefixN, all, allN,
		      setparameter, finalparameter, getparameter,
		      NULL, -10, false);

  finalparameter(isGLOBAL);
}



void detachmiraculix() {
  Ext_detachRFoptions(prefixlist, prefixN);
}

void RelaxUnknownRFoption(int *RELAX) { 
  Ext_relaxUnknownRFoption((bool) *RELAX); 
}
