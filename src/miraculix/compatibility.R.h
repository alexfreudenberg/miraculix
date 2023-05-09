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

#ifndef compatibility_to_R_h 
#define compatibility_to_R_h 1

#if defined compatibility_to_C_h
#error melange of R and C in header file calls of compatibility_to_*.h
#endif

// #include "compatibility.R.h"


#define RF_NA NA_REAL 
#define RF_NAN R_NaN
#define RF_NEGINF R_NegInf
#define RF_INF R_PosInf
#define ACOS std::acos
#define ASIN std::asin
#define ATAN std::atan
#define ATANH std::atanh
#define ACOSH std::acosh
#define ASINH std::asinh
#define EXPM1 std::expm1
#define LOG1P std::log1p
#define FROUND fround
#define COS std::cos
#define EXP std::exp
#define FABS(X) std::fabs((double) X) // OK; keine Klammern um X!
#if ! defined MALLOCX
#define MALLOCX std::malloc
#define FLOOR std::floor
#define SQRT(X) std::sqrt((double) X) // OK
#define CEIL(X) std::ceil((double) X) // OK; keine Klammern um X!
#define FREEX std::free
#endif
#define LOG std::log
#define POW(X, Y) R_pow((double) X, (double) Y) // OK; keine Klammern um X!
#define SIN std::sin
#define STRCMP(A, B) std::strcmp(A, B) // OK
#define STRCPY(A, B) std::strcpy(A, B) // OK
#define STRLEN std::strlen
#define STRNCMP(A, B, C) std::strncmp(A, B, C) // OK
#define TAN std::tan
#define MEMCOPYX std::memcpy
#define MEMMOVE std::memmove
#define MEMSET std::memset  
#define MEMCMP std::memcmp
#define AALLOC std::aligned_alloc
#define CALLOCX std::calloc
#define SPRINTF std::sprintf // Rprint 
#define QSORT std::qsort
#define RFERROR error
#define LOGGAMMAFN lgammafn
#define GAMMAFN gammafn
#define RPRINTF Rprintf
#define ATAN2 std::atan2
#define COSH std::cosh
#define SINH std::sinh
#define TANH std::tanh
#define SIGN(X) sign((double) X) // OK
#define TRUNC(X) truncf((double) X) // OK; keine Klammern um X!
#define REALX REAL

#define complex Rcomplex

#define LENGTH length // to avoid the unvoluntiered use of LENGTH defined by R
#define GAUSS_RANDOM(SIGMA) rnorm(0.0, SIGMA)
#define UNIFORM_RANDOM unif_rand()
#define POISSON_RANDOM(x) rpois(x)

#if ! defined USE_FC_LEN_T
#define USE_FC_LEN_T
#endif


#endif
