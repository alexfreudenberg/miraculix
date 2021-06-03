/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2018 -- 2019 Martin Schlather, 

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


//#include <R.h>
//#include <Rinternals.h>
//#include <Rdefines.h>
//#include <R_ext/Linpack.h>
//#include <stdio.h>  
//#include <stdlib.h>
//#include <unistd.h>
//#include <string.h>

// ACHTUNG: Reihenfolge nicht aendern!
#include "IntrinsicsBase.h"
#include "intrinsics.h"
#include "error.h"
#include "miraculix.h"
#include "options.h"
#include "AutoMiraculix.h"
#include "xport_import.h" // important! must be very last!
//                           (HAS_XXX in zzz_RandomFieldsUtils.h)
#include "mmagpu.h"
#include "utils.h"
#include "kleinkram.h"


 const char * prefixlist[prefixN] = {"genetics", "genetics_messgaes"};

// IMPORTANT: all names of general must have at least 3 letters !!!
const char *genetics[geneticsN] = 
  {"digits", "snpcoding", "centered", "normalized", "returnsigma", "efficient"};

const char *messages[messagesN] = 
  {"warn_address"}; 
  

option_type OPTIONS = { // OK
  genetics_START,
  messages_START
};


const char **all[prefixN] = {genetics,messages};
int allN[prefixN] = {geneticsN,messagesN};


void setoptions(Rint i, Rint j, SEXP el, char name[LEN_OPTIONNAME], 
		  bool VARIABLE_IS_NOT_USED isList, bool local) {
  if (!local && parallel())
    ERR("'RFoptions' may not be set from a parallel process.");
      
  option_type *options = WhichOptionList(local);
  
  switch(i) {
  case 0: {// genetics
    genetics_options *gp;
    gp = &(options->genetics);
    switch(j) {
    case 0: gp->digits = NUM; break;
    case 1: {
      Rint m;
      if (TYPEOF(el) != STRSXP) m = POS0NUM;
      else m= GetName(el, name, SNPCODING_NAMES, LastGenuineMethod + 1, gp->method);
#if !defined USEGPU
      if (m == MMAGPU && !HAS_CUDA) {
 	ERR1("'%.20s', which needs 'CUDA >= 11', is not available under the current compilation. See the starting message for a remedy.", SNPCODING_NAMES[m]);	
      }
#endif      
#if !defined SSE2
      if (m == Hamming2) {
 	ERR1("'%.20s', which needs 'SSE2', is not available under the current compilation. See the starting message for a remedy.", SNPCODING_NAMES[m]);
      }
      if (m == Packed ||  m == Multiply) {
	ERR1("'%.20s', which needs 'SSSE3' is not available under the current compilation. Set 'RFoptions(efficient=TRUE)' or see the starting message for a remedy.", SNPCODING_NAMES[m]);
      }
#endif
#if !defined SSSE3 
      if (m == Hamming3) {
	ERR1("'%.20s', which needs 'SSSE3' is not available under the current compilation. See the starting message for a remedy.", SNPCODING_NAMES[m]);
      }
      if (m == Shuffle) {
	ERR1("'%.20s', which needs 'SSSE3' is not available under the current compilation. Set 'RFoptions(efficient=TRUE)' or see the starting message for a remedy.", SNPCODING_NAMES[m]);
      }
#endif
#if !defined AVX2
      if (m == Packed256 || m == Multiply256) {
	ERR1("'%.20s', which needs 'AVX2', is not available under the current compilation. See the starting message for a remedy.",
	    SNPCODING_NAMES[m]);
      }
      if (m == Shuffle256) {
	ERR1("'%.20s', which needs 'AVX2', is not available under the current compilation. Set 'RFoptions(efficient=TRUE)' or see the starting message for a remedy.", SNPCODING_NAMES[m]);
      }
#endif
      if (m > LastGenuineMethod) ERR("given snp coding not allowed");
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
	Uint len = length(el);
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
	WARN2("'%.20s=TRUE' only allowed with 'centered=TRUE'.\n'%.20s=FALSE' is kept", genetics[i], genetics[i]); // OK
	gp->normalized = false;
      }
      break;
    case 4 : gp->returnsigma = LOGI; break;
    case 5 : gp->efficient = LOGI; break;
    default: BUG; 
    }}
    break;
  case 1: {
    messages_options *gp = &(options->messages);
    switch(j) {
    case 0: gp->warn_address = LOGI; break;
    default: BUG;
    }}
    break;    
  default: BUG; 
  }
}

   

void getoptions(SEXP sublist, Rint i, bool local) {
  Uint k = 0;
  option_type *options = WhichOptionList(local);
  switch(i) {
  case 0 : {
    genetics_options *p = &(options->genetics);
    ADD(ScalarReal(p->digits));
    //    ADD(ScalarString(mkChar(RELSHIP_METH_NAME[p->method])));
    ADD(ScalarInteger(p->method));
    ADD(ExtendedBooleanUsr(p->centered));    
    ADD(ScalarLogical(p->normalized));
    ADD(ScalarLogical(p->returnsigma));
    ADD(ScalarLogical(p->efficient));
  }
    break;
  case 1 : {
    messages_options *p = &(options->messages);
    ADD(ScalarLogical(p->warn_address));
  }
    break;   
  
  default : BUG;
  }
}


snpcoding getAutoCodingIntern() { // used also in options.h
  return
#if defined USEGPU
    HAS_CUDA ? MMAGPU :
#endif
#if defined AVX2
    Shuffle256
#elif defined SSSE3
    Shuffle
#elif defined SSE2
    Packed
#else
     TwoBit
#endif
    ;  
}


SEXP getAutoCoding() {
  SEXP dummy;
  PROTECT(dummy=allocVector(INTSXP, 1));
  INTEGER(dummy)[0] = getAutoCodingIntern();
  UNPROTECT(1);
  return dummy;
}


SEXP hasSSE2() {
  SEXP Ans;
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = HAS_SSE2;
  UNPROTECT(1);
  return Ans;
}


SEXP hasSSSE3() {
  SEXP Ans;
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = HAS_SSSE3;
  UNPROTECT(1);
  return Ans;
}


SEXP hasAVX() {
  SEXP Ans;
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = HAS_AVX;
  UNPROTECT(1);
  return Ans;
}


SEXP hasAVX2() {
  SEXP Ans;
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = HAS_AVX2;
  UNPROTECT(1);
  return Ans;
}

SEXP hasCUDA() {
  SEXP Ans;
  PROTECT(Ans = allocVector(LGLSXP, 1));
  LOGICAL(Ans)[0] = HAS_CUDA;
  UNPROTECT(1);
  return Ans;
}
