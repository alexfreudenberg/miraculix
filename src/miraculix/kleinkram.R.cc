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


#include "def.h" // nur wg STAND_ALONE
#include "intrinsics.h"
#include "Basic_RandomFieldsUtils.h"
#include "kleinkram.h"



#if defined compatibility_to_R_h

#include "kleinkram.h"
#include "General_utils.h"
// R FU INCLUDE #include  "xport_import_RFU.h"


const char // constants cannot be exported; 
*KKR_TYPE_NAMES[LAST_R_TYPE_NAME + 1] = { // never change ! see AutoRFU.cc
  "NILSXP" /* 0 */,
  "SYMSXP", "LISTSXP", "CLOSXP", "ENVSXP", "PROMSXP",
  "LANGSXP", "SPECIALSXP", "BUILTINSXP", "CHARSXP", "LGLSXP" /* 10 */,
  "??", "??", "INTSXP", "REALSXP", "CPLXSXP",
  "STRSXP", "DOTSXP", "ANYSXP", "ECSXP", "EXPRSXP" /*20 */,
  "BCODESXP", "EXTPTRSXP", "WEAKREFSXP", "RAWSXP", "S4SXP" /* 25 */,
  "", "", "", "", "NEWSXP" /* 30 */,
  "FREESXP", "??SXP"};


SEXP TooLarge(int *n, Long l){
#define nTooLarge 2 // mit op
  const char *tooLarge[nTooLarge] = {"size", "msg"};
  SEXP namevec, msg;
  PROTECT(msg=allocVector(VECSXP, nTooLarge));
  PROTECT(namevec = allocVector(STRSXP, nTooLarge));
  for (Long i=0; i<nTooLarge; i++)
    SET_STRING_ELT(namevec, i, mkChar(tooLarge[i]));
  setAttrib(msg, R_NamesSymbol, namevec);
  Long i=0;
  SET_VECTOR_ELT(msg, i++, Int(n, l));
  SET_VECTOR_ELT(msg, i,
		 mkString("too many elements - increase max.elements"));
  UNPROTECT(2);
  return msg;
}
SEXP TooLarge(Long n){int nn=(int) n;  return TooLarge(&nn, 1); }
SEXP TooLarge(Long row, Long col){
   int nn[2] = {(int) row, (int) col};
   return TooLarge(nn, 2);
}
 
SEXP TooSmall(){
  SEXP namevec;
  const char *msg = "value has not been initialized";
  PROTECT(namevec = allocVector(STRSXP, 1));
  SET_STRING_ELT(namevec, 0, mkChar(msg));
  UNPROTECT(1);
  return namevec;
}


SEXP Int(int *V, Long n, Long max) {
  SEXP Ans;
  if (V==NULL) return allocVector(INTSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();   
  PROTECT(Ans=allocVector(INTSXP, (int) n));
  MEMCOPY(INTEGER(Ans), V, n * sizeof(*V));
  UNPROTECT(1);
  return Ans;
}

SEXP Int(int* V, Long n) { return Int(V, n, n); }


SEXP Logic(bool* V, Long n, Long max) {
  SEXP Ans;
  if (V==NULL) return allocVector(VECSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();
  PROTECT(Ans=allocVector(LGLSXP, (int) n));
  int *ans = LOGICAL(Ans);
  for (Long i=0; i<n; i++) ans[i] = V[i];
  UNPROTECT(1);
  return Ans;
}
SEXP Logic(bool* V, Long n) { return Logic(V, n, n); }


SEXP Num(double* V, Long n, Long max) {
  SEXP Ans;
  if (V==NULL) return allocVector(REALSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();
  PROTECT(Ans=allocVector(REALSXP, (int) n));
  MEMCOPY(REAL(Ans), V, n * sizeof(*V));
  UNPROTECT(1);
  return Ans;
}
SEXP Num(double* V, Long n) {  return Num(V, n, n); }


SEXP Char(const char **V, Long n, Long max) {
  SEXP Ans;
  if (V==NULL) return allocVector(STRSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();
  PROTECT(Ans=allocVector(STRSXP, (int) n));
  for (Long i=0; i<n; i++) SET_STRING_ELT(Ans, i, mkChar(V[i]));  
  UNPROTECT(1);
  return Ans;
}

SEXP Char(const char **V, Long n) { return Char(V, n, n); }

SEXP Mat(double* V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP Ans;
  PROTECT(Ans=allocMatrix(REALSXP, (int) row, (int) col));
  MEMCOPY(REAL(Ans), V, n * sizeof(*V));
  UNPROTECT(1);
  return Ans;
}

SEXP Mat(double* V, Long row, Long col) {
  return Mat(V, row, col, MAXINT);
}


SEXP Mat_t(double* V, Long row, Long col, Long max) {
  Ulong ldAns = col;
  if (V==NULL) return allocMatrix(REALSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP Ans;
  PROTECT(Ans=allocMatrix(REALSXP, (int) row, (int) col));
  double *ans = REAL(Ans);
  for (Long k=0, j=0; j<col; j++)
     for (Long i=0; i<row; i++) ans[k++] = V[j + ldAns * i];
   UNPROTECT(1);
  return Ans;
}

SEXP Mat_t(double* V, Long row, Long col) {
  return Mat_t(V, row, col, MAXINT);
}


SEXP MatString(char **V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(STRSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP Ans;
  PROTECT(Ans=allocMatrix(STRSXP, (int) row, (int) col));  
  for (Long k=0; k<n; k++) SET_STRING_ELT(Ans, k, mkChar(V[k]));
  UNPROTECT(1);
  return Ans;
}

SEXP MatString(char** V, Long row, Long col) {
  return MatString(V, row, col, MAXINT);
}

SEXP MatInt(int* V, Long row, Long col, Long max) {
  if (V==NULL) return allocMatrix(INTSXP, 0, 0);
  Long n = row * col;
  if (n>max) return TooLarge(row, col);
  SEXP Ans;
  PROTECT(Ans=allocMatrix(INTSXP, (int) row, (int) col));
  MEMCOPY(INTEGER(Ans), V, n * sizeof(*V));
  UNPROTECT(1);
  return Ans;
}

SEXP MatInt(int* V, Long row, Long col) {
  return MatInt(V, row, col, MAXINT);
}

SEXP Array3D(double** V, Long depth, Long row, Long col, Long max) {
  if (V==NULL) return alloc3DArray(REALSXP, 0, 0, 0);
  Long
    m = row * col,
    n = row * col * depth;
  if (n>max) {
    int nn[3] = { (int) row, (int) col, (int) depth };
    return TooLarge(nn, 3);
  }
  SEXP Ans;
  PROTECT(Ans=alloc3DArray(REALSXP, (int) depth, (int) row, (int) col));
  double *ans = REAL(Ans);
  for (Long j=0; j<depth; j++) 
    for (Long i=0; i<m; i++) 
      ans[j*m+i] = V[j][i];
  UNPROTECT(1);
  return Ans;
}

SEXP Array3D(double** V, Long depth, Long row, Long col) {
  return Array3D(V, depth, row, col, MAXINT);
}




usr_bool UsrBoolRelaxed(SEXP p, char *name, Long idx) {
  double tmp = Real(p, name, idx);
  if (!R_finite(tmp)) return Nan;
  return tmp==0.0 ? False : True ;
}


usr_bool UsrBool(SEXP p, char *name, Long idx) {
  double tmp = Real(p, name, idx);
  usr_bool ans;
  if (tmp == 0.0) ans= False;
  else if (tmp == 1.0) ans= True;
  else if (ISNAN(tmp)) ans= Nan;
  else  RFERROR2("invalid value (%d) for boolean variable '%.50s'.",
		 (int) tmp, name);
  return ans;
}




SEXP String(char *V) {
  SEXP str;
  PROTECT(str = allocVector(STRSXP, 1)); 
  SET_STRING_ELT(str, 1, mkChar(V));
  UNPROTECT(1);
  return str;
}

SEXP String(char V[][MAXCHAR], Long n, Long max) {
  SEXP str;
  if (V==NULL) return allocVector(STRSXP, 0);
  if (n>max) return TooLarge(n);
  if (n<0) return TooSmall();
  PROTECT(str = allocVector(STRSXP, (int) n)); 
  for (Long i=0; i<n; i++) SET_STRING_ELT(str, i, mkChar(V[i]));
  UNPROTECT(1);
  return str;
}

SEXP String(char V[][MAXCHAR], Long n) {return String(V, n, n);}
 

SEXP String(int *V, const char * List[], Long n, Long endvalue) {
  SEXP str;
  if (V==NULL || n <= 0) return allocVector(STRSXP, 0);
  Long k;
  for (k=0; k<n; k++) {
    //    printf("k=%d %d: %d %d\n", k,n, V[k], NoInversionMethod);
    assert(V[k] <= endvalue);
    if (V[k] == endvalue) break;
  }
  PROTECT(str = allocVector(STRSXP, (int) k)); 
  for (Long i=0; i<k; i++) SET_STRING_ELT(str, i, mkChar(List[V[i]]));
  UNPROTECT(1);
  return str;
}

double Real(SEXP p, char *name, Long idx) {
  if (p != R_NilValue) {
    assert(idx < length(p));
    switch (TYPEOF(p)) {
    case REALSXP :  return REAL(p)[idx];
    case INTSXP :
      if (INTEGER(p)[idx]==NA_INTEGER) return RF_NA;
      else return((double) INTEGER(p)[idx]);
    case LGLSXP :
      if (LOGICAL(p)[idx]==NA_LOGICAL) return(RF_NA);
      else return((double) LOGICAL(p)[idx]);
    default : {}
    }
  }

  RFERROR2("'%.50s' can not be transformed to double! (got'%.50s')\n",
	   name,
	   TYPEOF(p) <= LAST_R_TYPE_NAME ? KKR_TYPE_NAMES[TYPEOF(p)] :
	   "something else"
	   );  
  return RF_NA;  // to avoid warning from compiler
}



void Real(SEXP el,  char *name, double *vec, Long maxn) {
  if (el == R_NilValue) {
    RFERROR1("'%.50s' cannot be transformed to double.\n", name);
  }
  Long n = length(el);
  for (Long j=0, i=0; i<maxn; i++) {
    vec[i] = Real(el, name, j);
    if (++j >= n) j=0;
  }
  return;
}


int Integer(SEXP p, char *name, Long idx, bool nulltoNA) {
  //printf("integer %s %d %d len=%d\n",  name, idx, nulltoNA, length(p));
  if (p != R_NilValue) {
    assert(idx < length(p));
    switch(TYPEOF(p)) {
    case INTSXP : 
      return INTEGER(p)[idx]; 
    case REALSXP : 
      double value;
      value = REAL(p)[idx];
      if (ISNAN(value)) {
	return NA_INTEGER;
      }
      int intvalue;
      intvalue = (int) value;
      if (value == intvalue) return intvalue;      
      else {
	RFERROR2("%.50s: integer value expected. Got %10e.", name, value);
      }
    case LGLSXP :
      if (LOGICAL(p)[idx]==NA_LOGICAL) return(NA_INTEGER);
      else return((int) LOGICAL(p)[idx]);
    default : {}
    }
  } else if (nulltoNA) return NA_INTEGER;

  RFERROR2("%.50s: incorrect type. Got '%.50s'.",
	   name,
	   TYPEOF(p) <= LAST_R_TYPE_NAME ? KKR_TYPE_NAMES[TYPEOF(p)]
	   : "something else");

 return NA_INTEGER; // compiler warning vermeiden
}

int Integer(SEXP p, char *name, Long idx) {
  return Integer(p, name, idx, false);
}


void Integer(SEXP el, char *name, int *vec, Long maxn) {
  if (el == R_NilValue) {
    RFERROR1("'%.50s' cannot be transformed to integer.\n",name);
  }
  Long n = length(el);
  for (Long j=0, i=0; i<maxn; i++) {
    vec[i] = Integer(el, name, j);
    if (++j >= n) j=0;
  }
}




void Integer2(SEXP el, char *name, int *vec) {
  Long n = length(el);
  if (n == 0)  RFERROR1("'%.50s' cannot be transformed to integer.\n",name);
 
  vec[0] = Integer(el, name, 0);
  if (vec[0] != NA_INTEGER && vec[0] < 1)
    RFERROR1("first component of '%.50s' must be at least 1", name);
  if (n == 1) vec[1] = vec[0];
  else {
    vec[1] = Integer(el, name, n-1);    
    if ( vec[1] != NA_INTEGER && vec[1] < vec[0])
      RFERROR1("'%.50s' must be increasing", name);
    if (n > 2) {
      vec[1] = vec[0];
      for (Long i = 1; i<n; i++)
	if (Integer(el, name, i) != ++(vec[1]))
	  RFERROR1("'%.50s' is not a sequence of numbers",name);
    }
  }
}





bool Logical(SEXP p, char *name, Long idx) {
   if (p != R_NilValue)
    assert(idx < length(p));
   switch (TYPEOF(p)) {
    case REALSXP:
      if (ISNAN(REAL(p)[idx])) return NA_LOGICAL ;
      else return (bool) REAL(p)[idx];
    case INTSXP :
      if (INTEGER(p)[idx]==NA_INTEGER) return NA_LOGICAL;
      else return (bool) INTEGER(p)[idx];
    case LGLSXP : return LOGICAL(p)[idx];
    default : {}
    }
  RFERROR1("'%.50s' cannot be transformed to logical.\n", name);  
  return NA_LOGICAL;  // to avoid warning from compiler
}


char Char(SEXP el, char *name) {
  SEXPTYPE type;
  if (el == R_NilValue) goto ErrorHandling;
  type = TYPEOF(el);
  if (type == CHARSXP) return CHAR(el)[0];
  if (type == STRSXP) {
    if (length(el) == 1) {
      if (STRLEN(CHAR(STRING_ELT(el,0))) == 1)
	return (CHAR(STRING_ELT(el,0)))[0];
      else if (STRLEN(CHAR(STRING_ELT(el,0))) == 0)
	return '\0';
    }
  }
 
 ErrorHandling:
  RFERROR1("'%.50s' cannot be transformed to character.\n",  name);  
  return 0; // to avoid warning from compiler
}


void String(SEXP el, char *name, char names[][MAXCHAR], Long maxlen) {
  Long l = length(el);
  SEXPTYPE type;  
  if (el == R_NilValue) goto ErrorHandling;
  if (l > maxlen)  {
    RFERROR1("number of variable names exceeds %d. Take abbreviations?",
	     (int) maxlen);
  }
  type = TYPEOF(el);
  if (type == CHARSXP) {
    for (Long i=0; i<l; i++) {
      names[i][0] = CHAR(el)[i];
      names[i][1] = '\0';
    }
  } else if (type == STRSXP) {
    for (Long i=0; i<l; i++) {
      strcopyN(names[i], CHAR(STRING_ELT(el, i)), MAXCHAR);
    }
  } else goto ErrorHandling;
  return;
 
 ErrorHandling:
  RFERROR1("'%.50s' cannot be transformed to character.\n",  name);  
}


int NonNegInteger(SEXP el, char *name) {
  int num = INT;
  if (num<0) {
    num=0; 
    WARN1("'%.50s', which has been negative, is set 0.\n",name);
  }
  return num;
}

double NonNegReal(SEXP el, char *name) {
  double num = NUM;
  if (num<0.0) {
    num=0.0; 
    WARN1("%.50s, which has been negative, is set 0.\n",name);
   }
  return num;
}

double NonPosReal(SEXP el, char *name) {
  double num = NUM;
  if (num>0.0) {
    num=0.0; 
    WARN1("%.50s, which has been positive, is set 0.\n",name);
  }
  return num;
}

int PositiveInteger(SEXP el, char *name) {
  int num = INT;
  if (num <= 0) {
    WARN2("'%.50s', which has been %.50s, is set 1.\n",
	  name, num ? "negative" : "0");
    num=1;
  }
   return num;
}

double PositiveReal(SEXP el, char *name) {
  double num = NUM;
  if (num<=0.0) {
     WARN2("'%.50s', which has been %.50s, is set 1.\n",
	   name, num==0.0 ? "0" : "negative");
     num=1.0; 
   }
  return num;
}

SEXP ExtendedInteger(double x) {
  return ScalarInteger(R_FINITE(x) ? x : NA_INTEGER);
}

SEXP ExtendedBooleanUsr(usr_bool x) {
  return ScalarLogical((int) x);
}




void GetName(SEXP el, char *name, const char * List[], Ulong n,
	     int defaultvalue, int endvalue, int *ans, Ulong maxlen_ans) {
   int
    k = 0, // globale variable !
    len_el = length(el);
   bool
     relax = defaultvalue >= 0;

  if (len_el > maxlen_ans) 
    RFERROR2("option '%.50s' is too lengthy. Maximum length is %lu.",
	     name, maxlen_ans);
  
  if (TYPEOF(el) == STRSXP) {    
    for (; k<len_el; k++) {
      ans[k] = Match(name, (char*) CHAR(STRING_ELT(el, k)), List, n,
		     '\0', NULL, relax);
      if (ans[k] < 0) {
	assert(relax); // otherwise Match would have failed
	ans[0] = defaultvalue;
	for (int i=1; k<maxlen_ans; i++) ans[i] = endvalue;
	return;
      }
    }
  } else {
    Integer(el, name, ans, maxlen_ans);
    for (; k<len_el; k++) {
      if (ans[k] < 0 || ans[k] >= n) matchError(name, NULL, ans[k], List, n, 0);
    }
  }
  
  for (k=len_el; k<maxlen_ans; k++) ans[k] = endvalue;
}


int GetName(SEXP el, char *name, const char * List[], int n,
	    int defaultvalue) {
  int i;
  GetName(el, name, List, n, defaultvalue, defaultvalue, &i, 1);
  return i;
}


int GetName(SEXP el, char *name, const char * List[], int n) {
 return GetName(el, name, List, n, -1);
}


#endif



