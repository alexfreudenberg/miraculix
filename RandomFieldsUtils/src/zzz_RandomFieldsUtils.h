


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017 Martin Schlather

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

#ifndef rfutils_init_H
#define rfutils_init_H 1


 /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
  */

#include "General_utils.h"
#include "Utils.h"
#include "options.h"


#ifdef HAVE_VISIBILITY_ATTRIBUTE
  # define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
  # define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MY_PACKAGE "RandomFieldsUtils"
#define MY_ACRONYM XX
#include "zzz_calls.h"

  /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
!!!!! auch kein mit MALLOC kreiertes Objekt  !!!!
  */

  
  DECLARE1(void, utilsoption_DELETE, utilsoption_type *, S)
  DECLARE1(void, utilsoption_NULL, utilsoption_type *, S)
  DECLARE1(void, solve_DELETE, solve_storage**, S)
  DECLARE1(void, solve_NULL, solve_storage*, x)
  DECLARE7(int, solvePosDef, double*, M, int, size, bool, posdef, 
	   double *, rhs, int, rhs_cols, double *, logdet, solve_storage *, PT)
  DECLARE8(int, solvePosDefSp, double *, M, int, size, bool, posdef,
	   double *, rhs, int, rhs_cols, double *,logdet,
	   solve_storage *, Pt, solve_options *,sp)
  //  DECLARE8(int, solvePosDefResult, double*, M, int, size, bool, posdef, 
  //	   double *, rhs, int, rhs_cols, double *, result, double*, logdet, 
  //	   solve_storage*, PT)
  DECLARE4(int, sqrtPosDefFree, double *, M, int, size, solve_storage *, pt,
	   solve_options *, sp)
  DECLARE3(int, sqrtRHS, solve_storage *, pt, double*, RHS, double *, res)
  DECLARE2(int, invertMatrix, double *, M, int, size)
  DECLARE2(double, StruveH, double, x, double, nu)
  DECLARE3(double, StruveL, double, x, double, nu, bool, expScale1d)
  DECLARE1(double, I0mL0, double, x)
  DECLARE3(double, WM, double, x, double, nu, double, factor)
  DECLARE3(double, DWM, double, x, double, nu, double, factor)
  DECLARE3(double, DDWM, double, x, double, nu, double, factor)
  DECLARE3(double, D3WM, double, x, double, nu, double, factor)
  DECLARE3(double, D4WM, double, x, double, nu, double, factor)
  DECLARE4(double, logWM, double, x, double, nu1, double, nu2, double, factor)
  DECLARE1(double, Gauss, double, x)
  DECLARE1(double, DGauss, double, x)
  DECLARE1(double, DDGauss, double, x)
  DECLARE1(double, D3Gauss, double, x)
  DECLARE1(double, D4Gauss, double, x)
  DECLARE1(double, logGauss, double, x)
  DECLARE0(int, cores1)
  DECLARE0(int, cpus)
  
  DECLARE1(void, getUtilsParam, utilsoption_type **, up)
  DECLARE14(void, attachRFoptions, char *, name,
	    const char **, prefixlist, int, N, 
	   const char ***, all, int *, allN, setoptions_fctn, set, 
	    finalsetoptions_fctn, final, getoptions_fctn, get,
	    deleteoptions_fctn, del, int, PLoffset, bool, basicopt,
	    install_modes, avx_needs, install_modes, gpu_needs, Uint, avx_info)
  DECLARE2(void, detachRFoptions, const char **, prefixlist, int, N)

  DECLARE3(void, sorting, double*, data, int, len, usr_bool, NAlast)
  DECLARE3(void, sortingInt, int*, data, int, len, usr_bool, NAlast)
  DECLARE4(void, ordering, double*, data, int, len, int, dim, int *, pos)
  DECLARE4(void, orderingInt, int*, data, int, len, int, dim, int *, pos)
  DECLARE4(double, scalarX, double *, x, double *, y, int, len, int, n)
  //  DECLARE4(int, scalarInt, int *, x, int *, y, int, len, int, n)
  DECLARE2(double, detPosDef, double *, M, int, size) // destroys M!
  DECLARE3(double, detPosDefsp, double *, M, int, size, solve_options *, sp)
  DECLARE8(int, XCinvXdet,double*, M, int, size, double *,X, int, X_cols,
	  double *, XCinvX, double *, det, bool, log, solve_storage, *PT)
  DECLARE10(int, XCinvYdet,double*, M, int, size, bool, posdef,
	    double *, X, double *, Y, int, cols,
	    double *, XCinvY, double *, det, bool, log, solve_storage, *PT)
  //  DECLARE5(double, XCinvXlogdet, double *, M, int, size, double *, X,
  //	   int, X_cols, solve_storage *, PT)
  DECLARE2(bool, is_positive_definite, double *, C, int, dim)
  DECLARE2(void, chol2inv, double *, MPT, int, size)
  DECLARE2(int, chol, double *, MPT, int, size)
  DECLARE1(void, pid, int *, i)
  DECLARE1(void, sleepMicro, int *, i)
  // DECLARE7(int, cholGPU, bool, copy, double*, M, int, size, double*, rhs, int, rhs_cols, double *, LogDet, double *, RESULT); // entkommentieren
  

 /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
  */
 
  /*

    See in R package RandomFields, /src/userinterfaces.cc 
          CALL#(...)
    at the beginning for how to make the functions available
    in a calling package

  */
#ifdef __cplusplus
}
#endif


/*

install.packages("RandomFieldsUtils_0.5.21.tar.gz", configure.args="CXX_FLAGS=-march=native", repo=NULL); library(RandomFieldsUtils, verbose=100); q()

*/

#ifdef DO_PARALLEL
#define HAS_PARALLEL true
#else
#define HAS_PARALLEL false
#endif
#if defined USEGPU
#define HAS_GPU true
#else
#define HAS_GPU false
#endif
#if defined AVX2
#define HAS_AVX2 true
#else
#define HAS_AVX2 false
#endif
#if defined AVX
#define HAS_AVX true
#else
#define HAS_AVX false
#endif
#if defined SSSE3
#define HAS_SSSE3 true
#else
#define HAS_SSSE3 false
#endif
#if defined SSE2
#define HAS_SSE2 true
#else
#define HAS_SSE2 false
#endif
#if defined SSE
#define HAS_SSE true
#else
#define HAS_SSE false
#endif

#define USES_GPU (HAS_GPU && NEED_GPU)
#define USES_AVX2 (HAS_AVX2 && NEED_AVX2)
#define USES_AVX (HAS_AVX && NEED_AVX)
#define USES_SSSE3 (HAS_SSSE3 && NEED_SSSE3)
#define USES_SSE2 (HAS_SSE2 && NEED_SSE2)
#define USES_SSE (HAS_SSE && NEED_SSE)

#define MISS_GPU (!HAS_GPU && NEED_GPU)
#define MISS_AVX2 (!HAS_AVX2 && NEED_AVX2)
#define MISS_AVX (!HAS_AVX && NEED_AVX)
#define MISS_ANY_SIMD (MISS_AVX2 || MISS_AVX || !HAS_SSE2)
#define MISS_SSSE3 (!HAS_SSSE3 && NEED_SSSE3)
#define MISS_SSE2 (!HAS_SSE2 && NEED_SSE2)
#define MISS_SSE (!HAS_SSE && NEED_SSE)


#define AttachMessageN 2000
#define HAS_ONE_RELEVANT (HAS_PARALLEL || USES_AVX2 || USES_AVX || USES_SSSE3 ||  USES_SSE2 || USES_SSE)
#define HAS_ALL_RELEVANT (HAS_PARALLEL && !MISS_AVX2 && !MISS_AVX && !MISS_SSSE3 &&  !MISS_SSE2 && !MISS_SSE)


#define AVX_INFO	  \
  (HAS_ONE_RELEVANT * 1 + \
   USES_GPU  *  (1<<1) + \
   USES_AVX2 * (1<<2) +	  \
   USES_AVX  * (1<<3) +	  \
   USES_SSSE3 * (1<<4) +	  \
   USES_SSE2  * (1<<5) +	  \
   USES_SSE  * (1<<6) +				     \
   (HAS_ONE_RELEVANT && !HAS_ALL_RELEVANT) * (1<<10) +  \
   MISS_GPU  * (1<<11) +				     \
   MISS_AVX2 * (1<<12)+				     \
   MISS_AVX  * (1<<13) +				     \
   MISS_SSSE3 * (1<<14) +				     \
   MISS_SSE2  * (1<<15) +				     \
   MISS_SSE * (1<<16))



#if defined WIN32
  #define AVX_NEEDS Inone
  #define GPU_NEEDS Inone
#else
#define AVX_NEEDS					      \
 (MISS_AVX2 ? Iavx2 : MISS_AVX ? Iavx : MISS_SSSE3 ? Issse3 : \
   MISS_SSE2 ? Isse2 :  MISS_SSE ? Isse : Inone)
  #define GPU_NEEDS (MISS_GPU ? Igpu : Inone) 
#endif


#define EAX 0
#define EBX 1
#define ECX 2
#define EDX 3
#define sse Available(1, EDX,25)
#define sse2 Available(1, EDX,26)
#define sse3 Available(1, ECX,0)
#define ssse3 Available(1, ECX,9)
#define sse4a Available(1, ECX,6)
#define fma3 Available(1, ECX,12)
#define avx Available(1, ECX,28)
#define avx2 Available(7, EBX,5)
#define avx512F Available(7, EBX,16)
#define avx512PF Available(7, EBX,26)
#define avx512AR Available(7, EBX,27)
#define avx512CD Available(7, EBX,28)

#ifdef _WIN32
#define AVAILABLE							\
  static inline bool Available(unsigned Blatt, int Register, int Bit) {	\
    uint32_t s[4];							\
    __cpuid((int *)s, (int) Blatt);						\
    return s[Register] & (1 << (Bit));					\
  }
#else 
#define AVAILABLE							\
static inline bool Available(unsigned Blatt, int Register, int Bit) {	\
  uint32_t s[4];							\
  asm volatile								\
    ("cpuid": "=a"(s[0]), "=b"(s[1]),"=c"(s[2]),"=d"(s[3]):"a"(Blatt),"c"(0)); \
  return s[Register] & (1 << (Bit));					\
}
#endif

#define WARN_PARALLELXX							\
  if (OPTIONS_UTILS->basic.warn_parallel && mypid == parentpid) 	\
    PRINTF("Do not forget to run 'RFoptions(storing=FALSE)' after each call of a parallel command (e.g. from packages 'parallel') that calls a function in 'RandomFields'. (OMP within RandomFields is not affected.) This message can be suppressed by 'RFoptions(warn_parallel=FALSE)'.") /*// ok */ \
    
#define WARN_PARALLEL 


#endif
