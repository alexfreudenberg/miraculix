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



#ifndef compatibility_to_C_h 
#define compatibility_to_C_h 1

#include "errors_messages.h"


#if defined compatibility_to_R_h
#error melange of C and R in header file calls of compatibility_to_*.h
#endif

// #include "compatibility.C.h"
//#include "options_RFU.h"

extern double ownNA;
extern double ownNaN;
extern int errno ;

#if !defined INT_MIN
#define INT_MIN (-2147483647-1)
//#define INT_MIN (-2147483648)
#endif
#define NA_INTEGER INT_MIN  // NA_INTEGER OK
#define R_FINITE(X) (FABS(X) <= 1e300)
#define ISNA(X) ((X) == RF_NA)
#define isNaN(X) ((X) == RF_NAN)
#define ISNAN(X) (ISNA(X) || isNaN(X))
#define warning(X) PRINTF("WARNING: %s", X);
// void warning(const char* format, ...);
#define R_NilValue NULL

#define NA_REAL ownNA
#define RF_NA NA_REAL
#define RF_NAN ownNaN
#define RF_NEGINF (-INFINITY)
#define RF_INF INFINITY
#define ACOS acos
#define ASIN asin
#define ATAN atan
#define ATANH atanh
#define ACOSH acosh
#define ASINH asinh
#define EXPM1 expm1
#define LOG1P log1p
#define FROUND fround
#define COS cos
#define EXP exp
#define FABS(X) fabs((double) X) // OK; keine Klammern um X!
#if ! defined MALLOCX
#define MALLOCX malloc
#define FLOOR floor
#define SQRT(X) sqrt((double) X) // OK
#define CEIL(X) ceil((double) X) // OK; keine Klammern um X!
#define FREEX free
#endif
#define LOG log
#define POW(X, Y) pow((double) X, (double) Y) // OK; keine Klammern um X!
#define SIN sin
#define STRCMP(A, B) strcmp(A, B) // OK
#define STRCPY(A, B) strcpy(A, B) // OK
#define STRLEN strlen
#define STRNCMP(A, B, C) strncmp(A, B, C) // OK
#define TAN tan
#define MEMCOPYX memcpy
#define MEMMOVE memmove
#define MEMSET memset  
#define MEMCMP memcmp
#define AALLOC aligned_alloc
#define CALLOCX calloc
#define SPRINTF sprintf // Rprintf
#define QSORT qsort
#define RFERROR(X) {fprintf(stderr, "%s", X); exit(EXIT_FAILURE); }
#define LOGGAMMAFN lgamma
#define GAMMAFN tgamma
#define RPRINTF printf // Rprintf
#define ATAN2 atan2
#define COSH cosh
#define SINH sinh
#define TANH tanh
#define SIGN(X) (((X) > 0) - ((X) < 0))
#define TRUNC(X) trunc((double) X) // OK; keine Klammern um X!


double gauss_random(double mu, double sigma);
double uniform_random(double a, double b); 
double poisson_random(double lambda);
#define GAUSS_RANDOM(SIGMA) gauss_random(0.0, SIGMA)
#define UNIFORM_RANDOM uniform_random(0.0, 1.0) // OK
#define POISSON_RANDOM(x) poisson_random(x)

#if defined USE_FC_LEN_T
#undef USE_FC_LEN_T
#warning USE_FC_LEN_T was already defined
#endif


#define F77_NAME(F) F##_  // Basic_utils
#define F77_CALL(F) F##_  // Basic_utils


F77name(dgeqrf)(int*, int*, double* , int *, double*, double*, int*, int*);
F77name(dsyevr)(const char*,  const char *, const char *, int*, double* ,
		int*, double* , double* ,
		int *, int*, double*,
		int*, double* , double* ,
		int *, int*, double*,
		int*, int*, int*, int*);
F77name(dgetrf)(int*, int*, double* , int *, int*, int*);
F77name(dgetrs)(const char *, int*, int*, double* , int *, int*,double*, int*,
		int*);
F77name(dgetri)(int*,  double* , int *, int*,double*, int*, int*);
F77name(dgesv)(int*, int*, double* , int *, int*,double*, int*, int*);
F77name(dpotrf)(const char*, int*, double* , int *, int*);
F77name(dtrmv)(const char*,  const char *, const char *,
	       int*, double* , int*,double*, int*);

F77name(dgesdd)(const char*, int*, int*, double* ,
		int*, double* , double* ,
		int*, double* , int*, double* ,
		int *, int*, int*);
F77name(dgemv)(const char*,  int*, int*, double* , double*,
	       int*, double* ,
	       int*, double* , double*,
	       int *);
F77name(ddot)(int*, double* , int*, double*, int*);
F77name(dsyrk)(const char*,  const char *,
	       int*, int*, double* ,  double*,
	       int *, double*,double*,
	       int*);

#endif
