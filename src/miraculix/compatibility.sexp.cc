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




#include "def.h"
#include "Basic_RandomFieldsUtils.h"
#include "compatibility.SEXP.h"
#include "errors_messages.h"



SEXP
  Class = R_NilValue,
  Information = R_NilValue,
  Filename = R_NilValue,
  Filecoding = R_NilValue,
  Next = R_NilValue,
  Missings = R_NilValue,
  Doubled = R_NilValue,
  Precise = R_NilValue;


void *typeError(int type, int required) {
  char msg[100];
  SPRINTF(msg, "type should be %d, but is %d\n", required, type);
  RFERROR(msg);
  return NULL;
}

#define copyAttr(WHICH)  tmp=getAttribValue(From, WHICH); setAttrib(To, WHICH, tmp)


#define condCopy(WHICH) {				       \
    oldAttr = getAttribPointer(To, WHICH);		       \
    if (unconditional || oldAttr == R_NilValue) {	       \
      if (oldAttr != R_NilValue) { FREE_SEXP(&oldAttr); }      \
   copyAttr(WHICH);					       \
  }							       \
}


Long Memory(int type) {
  switch(type) {
  case INTSXP : return sizeof(int); 
  case LGLSXP : return sizeof(int); 
  case REALSXP : return sizeof(double);
  default:
#if defined compatibility_to_R_h
    BUG; // may never appear!
#else    
    switch(type) {
    case LONGREALSXP : return sizeof(LongDouble);    
    case LONGSXP :     return sizeof(Long);
    case POINTERSXP :  return sizeof(void*);
    default : RFERROR("R type unknown or not programmed yet");
    }
#endif
    return 0;
  }
}



#if defined compatibility_to_R_h


void install_default() {
  Information = install("information");
  Filecoding = install("filecoding");
  Filename = install("filename");
  Next = install("next");
  Precise = install(PreciseName);
  Doubled = install("longDouble");
  Missings = install("missings");
 }




SEXP allocVectorPointer( int type, Long len, void *x) {
  SEXP Ans = PROTECT(allocVector(type, len));
  MEMCOPY(VOIDSXP(Ans), x, Memory(type) * len);
  UNPROTECT(1);
  return Ans;
}

/*
SEXP allocVectorGeneralIntern(int type, Long len) {
    switch(type) {
    LONGSXP : ?????
    }
  SEXP Ans = allocObject(type, len); // PROTECT() 
  Ans->object = SXPvector;
 return Ans;
}
*/

#else



#include "compatibility.SEXP.h"


SEXP allocObject(int type, Long len, void *x) {
  SEXP Ans = (SEXP) CALLOC(1, sizeof(SEXPtype));
  Ans->type = type;
  
  if (type == STRSXP) {
    if (len != 1) BUG; // not programmed!
    return Ans;
  }
  Ans->len = len;
  //printf("****** type = %d %ld %d %ld\n", type, len, Memory(type), len * Memory(type));  
  Ans->pointer = x != NULL;
  if (Ans->pointer) {
    Ans->x = CALLOC(1, sizeof(void*));
    *((void**) (Ans->x)) = x;
  } else Ans->x = len == 0 ? R_NilValue : CALLOC(Memory(type), len);
  // printf("don %ld\n", (Long) (Ans->x));
  return Ans;
}

SEXP allocObject(int type, Long len) {
  return allocObject(type, len, NULL);
}

SEXP allocVectorPointer( int type, Long len, void *x) {
  SEXP Ans = allocObject(type, len, x); /* PROTECT() */
  Ans->object = SXPvector;
  return Ans;
}

SEXP allocVectorIntern(int type, Long len) {
  return allocVectorPointer(type, len, NULL);
}

SEXP allocMatrixPointer(int type, Long row, Long col, void *x) {
    SEXP Ans = allocObject(type, row * col, x ); /* PROTECT() */
  Ans->dim[0] = row;
  Ans->dim[1] = col;
  Ans->object = SXPmatrix;
  return Ans;
}

SEXP allocMatrixIntern(int type, Long row, Long col) {
  return allocMatrixPointer(type, row, col, NULL);
 }


void FREE_SEXP(SEXP *S){
  if (*S != NULL) {
    if ((*S)->len > 0 && !(*S)->pointer) { FREE((*S)->x); }
    FREE_SEXP(&((*S)->Class));
    FREE_SEXP(&((*S)->Information));
    FREE_SEXP(&((*S)->Filecoding));
    FREE_SEXP(&((*S)->Filename));
    FREE_SEXP(&((*S)->Next));
    FREE_SEXP(&((*S)->Missings));
    FREE_SEXP(&((*S)->Precise));
    FREE(*S);
  }
}


void SET_STRING_ELT(SEXP S, int elmnt, const char *c) {
  if (S->type != STRSXP) RFERROR("SEXP not a string element");
  if (elmnt != 0) RFERROR("only first element can be set in standalone mode");
  FREE(S->x);
  Long n = STRLEN(c) + 1;
  S->x = MALLOC(sizeof(char) * n);
  MEMCOPY(S->x, c, n * sizeof(char));
}

SEXP COPY_SEXP_STRING(SEXP X) {
  SEXP Ans = allocVector(STRSXP, 1);
  SET_STRING_ELT(Ans, 0, mkChar(CHAR(STRING_ELT(X, 0))));
  return Ans;
}


#endif



#define showSEXP(WHICH)				\
  if (getAttribPointer(S, WHICH) != R_NilValue) {	\
  for (int i=0; i<=level; i++) PRINTF("  ");		\
  PRINTF("%s:", ""#WHICH);				\
  printSEXP(getAttribPointer(S, WHICH), technical, level+1, n);	\
  }

void printSEXP(SEXP S, bool technical, int level, int n) {
  if (S == R_NilValue) { PRINTF("NilValue\n"); return; }
  const char *eq[]={"!=NULL", "=0"};
  // for (int i=0; i<level; i++) PRINTF("  ");
  PRINTF(" ");
  if (technical) PRINTF(" x%s len=%ld, dim=(%ld,%ld) type=%d dim=%d\n",
			eq[VOIDSXP(S)  == NULL],
			LENGTH(S),
			NROWS(S), NCOLS(S),
			TYPEOF(S),
			getDimension(S));
  else {
    Long dim = getDimension(S);
    int type = TYPEOF(S);
    Long cols = NCOLS(S);
    Long rows = NROWS(S);
    
    switch(type) {
    case INTSXP : PRINTF("int"); break;
    case REALSXP : PRINTF("real"); break;
    default:
      if (type == LONGREALSXP) PRINTF("longreal");
      else PRINTF("type=%d", type);   
    }
    
    if (dim == 0) {PRINTF("[] type = %d\n", type); return;}
    if (dim == 1) {Long t = rows; rows = cols; cols = t; PRINTF("[%ld] ", t);}
    else if (dim == 2) PRINTF("[%ld,%ld]\n", rows, cols);
    else PRINTF("[%ld,...]\n", rows);
    Long c = cols;
    Long r = rows;
    if (r > n + 1) r = n;
    if (c > n + 1) c = n;
    if (type != INTSXP && type != REALSXP && type != LONGREALSXP) {
      PRINTF("\n");
      return;
    }

    PRINTF(" ");
    for (Long j=0; j<r; j++) {
      if (dim > 1) for (int i=0; i<=level; i++) PRINTF("  ");
      for (Long i=0; i<c; i++) {
	switch(type) {
	case INTSXP :PRINTF("%d ", INTEGER(S)[i * r  + j]); break;
	case REALSXP : PRINTF("%f ", REAL(S)[i * r  + j]); break;
	default:
	  if (type == LONGREALSXP) PRINTF("%e ", (double)LONGREAL(S)[i*r+j]); 
	  else BUG;
	}
      }
      if (cols > c) PRINTF("...[%ld]", cols -c);
      PRINTF("\n");
    }
  }
  showSEXP(Information);
  showSEXP(Filecoding);
  showSEXP(Filename);
  showSEXP(Class);
  showSEXP(Next);
  showSEXP(Missings);
  showSEXP(Precise);
}


void printSEXP(SEXP S, int n) { PRINTF("\nSEXP:"); printSEXP(S, false, 0, n); }
void printSEXP(SEXP S) {printSEXP(S, 5);}


int objectError(int object, int required_type) {
  const char *objNames[4] = {"null", "vector", "matrix", "array"};
  char msg[100];
  if (required_type >= SXPempty && required_type <= SXParray)  {
    SPRINTF(msg, "object should be %s, but ", objNames[required_type-SXPempty]);
    if (object >= SXPempty && object <= SXParray)
      SPRINTF(msg, "is %s\n", objNames[object - SXPempty]);
    else 
      SPRINTF(msg, "is of type %d\n", object);
  } else BUG;
  RFERROR(msg);
  return 0;
}



SEXP copySEXP(SEXP From) {
   // printf("hier %ld\n", From==R_NilValue);
  SEXP To;
  if (From == R_NilValue) return R_NilValue;
  if (isVector(From)) To = allocVector(TYPEOF(From), LENGTH(From));
  else if (isMatrix(From))
    To = PROTECT(allocMatrix(TYPEOF(From), NROWS(From), NCOLS(From)));
  else BUG;
  MEMCOPY(VOIDSXP(To), VOIDSXP(From), LENGTH(From)*Memory(TYPEOF(From)));
  SEXP tmp;
  copyAttr(Information);
  copyAttr(Filecoding);
  copyAttr(Filename);
  copyAttr(Class);
  copyAttr(Next);
  copyAttr(Missings);
  copyAttr(Precise);
  UNPROTECT(1);
  return To;
}


void copySEXPattrib(SEXP To, SEXP From, bool unconditional) {
   //  printf("condSEXPattrib **** hier %d\n", From==R_NilValue);
  if (From == R_NilValue) return;
  SEXP tmp = R_NilValue;
  SEXP oldAttr = R_NilValue;

  condCopy(Information);
  condCopy(Filecoding);
  condCopy(Filename);
  condCopy(Class);
  condCopy(Next);
  condCopy(Missings);
  condCopy(Precise);
}

