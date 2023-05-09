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

#ifndef compatibility_SEXP_h 
#define compatibility_SEXP_h 1

#include "compatibility.general.h"

void *typeError(int type, int required);
int objectError(int object, int required);


#define SXPempty    60
#define SXPvector   61
#define SXPmatrix   62
#define SXParray    63
#define PreciseName "longDoublePrecision"

#if defined compatibility_to_R_h


#include <Rdefines.h>
#define FREE_SEXP(X)
#define getAttribPointer getAttrib
#define getAttribValue getAttrib
#define getDimension(X) (X == R_NilValue ? 0 :	       \
			 isVector(X) ? 1 :	       \
			 isMatrix(X) ? 2 :	       \
			 isArray(X) ? -1 : \
			 objectError(TYPEOF(X), SXPmatrix))
#define COPY_SEXP_STRING(X) ScalarString(mkChar(CHAR(STRING_ELT(X, 0))))
#define NROWS(X) (X==R_NilValue ? 0 : isVector(X) ? LENGTH(X) :	\
		  isMatrix(X) ? nrows(X) : objectError(TYPEOF(X), SXPmatrix))
#define NCOLS(X) (X==R_NilValue ? 0 : isVector(X) ? 1 : \
		  isMatrix(X) ? ncols(X) : objectError(TYPEOF(X), SXPmatrix))
#define VOIDSXP(X) (TYPEOF(X)==INTSXP ? (void*) INTEGER(X) :\
		    TYPEOF(X)==REALSXP ? (void*) REAL(X) :		\
		    TYPEOF(X)==LGLSXP ? (void*) LOGICAL(X) :		\
		    typeError(TYPEOF(X), BUGSXP))

#define BUGSXP      50
#define LONGREALSXP REALSXP
#define LONGSXP     INTSXP


#define LONG(X) ((Long *) INTEGER(X))
#define LONGREAL(X) ((LongDouble *) REAL(X))

#define allocVectorLONGREALSXP(Y) \
  allocVector(INTSXP, sizeof(LongDouble) / sizeof(double) * (Y))
#define allocVectorLONGSXP(Y) \
  allocVector(INTSXP, sizeof(Long) / sizeof(int) * (Y))
#define allocVectorExtra(X,Y) allocVector##X(Y)
    

#define allocMatrixLONGSXP(Y,Z)					\
  allocMatrix(INTSXP, sizeof(Long) / sizeof(int) * (Y), Z)
#define allocMatrixLONGREALSXP(Y,Z) \
  allocMatrix(INTSXP, sizeof(LongDouble) / sizeof(double) * (Y), Z)
#define allocMatrixExtra(X,Y,Z) allocMatrix##X(Y,Z)

void install_default();


#else


#define NILSXP	     0	 
#define SYMSXP	     1	  
#define LISTSXP	     2	  
#define CLOSXP	     3	  
#define ENVSXP	     4	  
#define PROMSXP	     5	  
#define LANGSXP	     6	  
#define SPECIALSXP   7	  
#define BUILTINSXP   8	  
#define CHARSXP	     9	  
#define LGLSXP	    10	  

#define INTSXP	    13	  
#define REALSXP	    14	  
#define CPLXSXP	    15	  
#define STRSXP	    16	  
#define DOTSXP	    17	  
#define ANYSXP	    18	  

#define VECSXP	    19	  
#define EXPRSXP	    20	  
#define BCODESXP    21    
#define EXTPTRSXP   22    
#define WEAKREFSXP  23    
#define RAWSXP      24    
#define S4SXP       25    

#define NEWSXP      30    
#define FREESXP     31

#define BUGSXP 50
#define LONGREALSXP 51
#define LONGSXP 52
#define POINTERSXP 53

//#define FLOATSXP    55
//#define BYTESXP     56


#define FUNSXP      99


typedef struct SEXPtype SEXPtype; // OK
typedef SEXPtype *SEXP; 
struct SEXPtype {
  void *x;
  bool pointer;
  Long len; // nr of integer values
  Long dim[2];
  int type;
  int object;
  SEXP Class, Information, Filecoding, Filename, Next, Missings, Precise;
};



#define PROTECT(X) X
#define UNPROTECT(X)
#define LENGTH(X) ((X)->len)
#define TYPEOF(X) ((X)->type)
#define POINTER_OF(X) ((X)->pointer ? *((void**)((X)->x)) : (X)->x) 
#define INTEGER(X) (TYPEOF(X) == INTSXP	 ? (int*) POINTER_OF(X) :	\
		    (int*) typeError(TYPEOF(X), INTSXP))
#define LONG(X) (TYPEOF(X) == LONGSXP ? (Long*) POINTER_OF(X) :\
		 (Long*) typeError(TYPEOF(X), LONGSXP))
#define LONGREAL(X) (TYPEOF(X) == LONGREALSXP ? (LongDouble*) POINTER_OF(X) : \
		    (LongDouble*) typeError(TYPEOF(X), LONGREALSXP))
#define LOGICAL(X) (TYPEOF(X) == LGLSXP ? (int*) POINTER_OF(X) : \
		    (int*) typeError(TYPEOF(X), LGLSXP))
#define STRING_ELT(X, N) (TYPEOF(X) == STRSXP && (N) > 0 ? POINTER_OF(X) : \
			  typeError(TYPEOF(X), LGLSXP))
#define VOIDSXP(X) POINTER_OF(X)
#define CHAR
// #define allocVector(T,L) (printf("%s %d\n", __FILE__, __LINE__),  allocVectorIntern(T, L))
#define allocVector allocVectorIntern
#define allocVectorExtra allocVectorIntern
#define allocMatrix allocMatrixIntern
#define allocMatrixExtra allocMatrixIntern
SEXP allocVectorIntern(int type, Long len);
SEXP allocMatrixIntern(int type, Long row, Long col);
#define NROWS(X) ((X)->dim[0])
#define NCOLS(X) ((X)->dim[1])
#define nrows(X) ((X)->object==SXPmatrix ? (int) NROWS(X)	\
		  : (int) objectError((X)->object, SXPmatrix))
#define ncols(X) ((X)->object==SXPmatrix ? (int) NCOLS(X) \
		  : (int) objectError((X)->object, SXPmatrix))
#define isMatrix(X) ((X)->object == SXPmatrix)
#define isVector(X) ((X)->object == SXPvector)
void FREE_SEXP(SEXP *S);
#define getAttribPointer(S,A) ((S)->A)
#define getAttribValue(S, A) copySEXP(getAttribPointer(S,A))
#define GET_CLASS(S) getAttribPointer(S,Class)
#define getDimension(S) ((S)->object-SXPempty)
#define mkChar (char*)
void SET_STRING_ELT(SEXP S, int elmnt, const char *);
#define setAttrib(S,A,T) {      	\
  if ((S)->A != T) {			\
    FREE_SEXP(&((S)->A));	       	\
    (S)->A = T;				\
  }					\
}
#define SET_CLASS(S,C) setAttrib(S,Class,C)

// in the miraculix package, the call of REAL() in a non-R-compatibility
// environment should be an exception. To control this, REALX is introduced.
#define REALX(X) (TYPEOF(X) == REALSXP ? (double*) POINTER_OF(X) : \
		  (double*) typeError(TYPEOF(X), REALSXP))
#define REAL REALX

SEXP COPY_SEXP_STRING(SEXP X);
SEXP copySEXP(SEXP X);

   
#endif


void printSEXP(SEXP S);
void printSEXP(SEXP S, int n);
void copySEXPattrib(SEXP To, SEXP From, bool unconditional);
SEXP allocVectorPointer(int type, Long len, void *x);
SEXP allocMatrixPointer(int type, Long row, Long col, void *x);
 

#endif
