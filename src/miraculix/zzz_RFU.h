
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


#ifndef rfutils_init_H
#define rfutils_init_H 1


 /* 
!!!!! HIER NIE EIN S E X P OBJEKT ZURUECKGEBEN  !!!!  
  */

#include "options_RFU.h"
#include "rfu.h"

#define AttachMessageN 2000
#define LEN_OPTIONNAME 201 // zwingend ungerade

 
typedef void (*finalsetoptions_fctn) (int);
typedef void (*deleteoptions_fctn) (bool);


#define MATERN_NU_THRES 100
#define BESSEL_NU_THRES 100
#define LOW_MATERN 1e-20
#define LOW_BESSEL 1e-20


typedef struct whittle_work_type { // OK
  double loggamma1old, nu1old,
    loggamma2old, nu2old,
    loggamma_old,nuOld,
    gamma, nuAlt;
} whittle_work_type;


#ifdef HAVE_VISIBILITY_ATTRIBUTE
  # define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
  # define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MY_PACKAGE "RandomFieldsUtils"
  //#define MY_ACRONYM XX
#include "zzz_calls.h"

  /* 
!!!!! HIER NIE EIN S E X P OBJEKT ZURUECKGEBEN  !!!!  
!!!!! auch kein mit MALLOC kreiertes Objekt  !!!!
  */


#if defined compatibility_to_R_h
  typedef void (*setoptions_fctn) (int, int, SEXP, char[LEN_OPTIONNAME],
				 bool, bool);
  typedef void (*getoptions_fctn) (SEXP, int, bool);
  DECLARE2(SEXP, RFUoptions, SEXP options, char * calling)
  DECLARE3(void, getoptionsRFU, SEXP sublist, int i, utilsoption_type *options)
  DECLARE6(void, setoptionsRFU, int i, int j, SEXP el,
	   char name[LEN_OPTIONNAME], bool isList, utilsoption_type *options) 
  DECLARE2(void, detachRFUoptions, const char ** prefixlist, int N)
  DECLARE18(void, attachRFUoptions, char * name,
	    const char ** prefixlist, int N, 
	    const char *** all, int * allN,
	    setoptions_fctn set, getoptions_fctn get,
	    finalsetoptions_fctn final,
	    deleteoptions_fctn del,
	    setoptions_fctn setRFU, getoptions_fctn getRFU,
	    int pl_offset, bool basicopt,
	    install_modes gpu_needs, Uint avx_info, int version,
	    int RFUversion, int mem_is_aligned)
  DECLARE4(void, attachSetNGet, char*calling,
	   char *pkgname, setoptions_fctn set,
	   getoptions_fctn get)
#endif

  
  DECLARE1(void, del_utilsoption, utilsoption_type * S)
#define PIVOT_IDX_N 0
#define  N_UTILS_PARAM (PIVOT_IDX_N + 1)
  DECLARE2(void, params_utilsoption,  int local, int * params)
  DECLARE2(void, get_utilsoption, utilsoption_type * S, int local)
  DECLARE2(void, get_utils_basic, basic_options * S, int local)
  DECLARE2(void, push_utilsoption, utilsoption_type * S, int local)
  DECLARE0(void, update_utilsoption)
  DECLARE0(void, startRFU)
  
   
  DECLARE1(void, solve_DELETE, solve_storage** S)
  DECLARE1(void, solve_NULL, solve_storage* x)
  DECLARE4(int, sqrtRHS, solve_storage * pt, double* RHS, double * res,
	   int cores)
 		   
  
  DECLARE2(double, StruveH, double x, double nu)
  DECLARE3(double, StruveL, double x, double nu, bool expScale1d)
  DECLARE1(double, I0mL0, double x)
  DECLARE4(double, WM, double x, double nu, double factor,
	   whittle_work_type *work)
  DECLARE4(double, DWM, double x, double nu, double factor,
	   whittle_work_type *work)
  DECLARE4(double, DDWM, double x, double nu, double factor,
	   whittle_work_type *work)
  DECLARE4(double, D3WM, double x, double nu, double factor,
	   whittle_work_type *work)
  DECLARE4(double, D4WM, double x, double nu, double factor,
	   whittle_work_type *work)
  DECLARE5(double, logWM2, double x, double nu1, double nu2, double factor,
	   whittle_work_type *work)
  DECLARE1(double, Gauss, double x)
  DECLARE1(double, DGauss, double x)
  DECLARE1(double, DDGauss, double x)
  DECLARE1(double, D3Gauss, double x)
  DECLARE1(double, D4Gauss, double x)
  DECLARE1(double, logGauss, double x)
  DECLARE0(int, cores1)
  DECLARE0(int, cpus)
  
 

  DECLARE3(void, sorting, double* data, int len, usr_bool NAlast)
  DECLARE3(void, sortingL, double* data, Long len, usr_bool NAlast)
  DECLARE3(void, sortingInt, int* data, int len, usr_bool NAlast)
  DECLARE3(void, sortingLong, Long* data, Long len, usr_bool NAlast)
  DECLARE4(void, ordering, double* data, int len, int dim, int * pos)
  DECLARE4(void, orderingL, double* data, Long len, int dim, Long * pos)
  DECLARE4(void, orderingInt, int* data, int len, int dim, int * pos)
  DECLARE4(void, orderingLong, Long* data, Long len, int dim, Long * pos)

  DECLARE4(double, scalarX, double * x, double * y, Long len, Long n)
  DECLARE4(Long, scalarXint, int * x, int * y, Long len, Long n)
   //  DECLARE4(int, scalarInt, int * x, int * y, int len, int n)
  DECLARE3(void, chol2inv, double * MPT, int size, int cores)
  DECLARE1(void, pid, int * i)
  DECLARE0(bool, parallel)
  DECLARE1(void, sleepMicro, int * i)
  // DECLARE7(int, cholGPU, bool copy, double* M, int size, double* rhs, int rhs_cols, double * LogDet, double * RESULT); // entkommentieren


  DECLARE9(int, xCinvXdet,double* M, int size, double *X, Long X_cols,
	   double * XCinvX, double * det, bool log, solve_storage *PT,
	   int cores)
  DECLARE11(int, xCinvYdet,double* M, int size, bool posdef,
	    double * X, double * Y, Long cols,
	    double * XCinvY, double * det, bool log, solve_storage *PT,
	    int cores)
  DECLARE3(int, cholesky, double * MPT, int size, int cores)
  DECLARE8(int, SolvePosDef, double* M, int size, bool posdef, 
	   double * rhs, Long rhs_cols, double * logdet, solve_storage * PT,
	   int cores)
  DECLARE9(int, SolvePosDefSp, double * M, int size, bool posdef,
	   double * rhs, Long rhs_cols, double *logdet,
	   solve_storage * Pt, solve_options *sp, int cores)
  DECLARE5(int, SqrtPosDefFree, double * M, int size, solve_storage * pt,
	   solve_options * sp, int cores)
  DECLARE4(double, DetPosDefsp, double * M, int size, solve_options * sp,
	   int cores)
  DECLARE3(int, InvertMatrix, double * M, int size, int cores)
  DECLARE3(double, DetPosDef, double * M,  int size, int cores) // destroys M!
  DECLARE3(bool, Is_positive_definite, double * C, int  dim, int cores)


  


#define cRFUMethods				\
  CDEF(sleepMilli,  1, int_arg),		\
    CDEF(sleepMicro, 1, int_arg),		\
    CDEF(pid, 1, int_arg),			\
    CDEF(hostname, 2, host_arg),		\
    CDEF(setCPUs, 1, int_arg),			\
    CDEF(recompilationNeeded, 1, int_arg),	\
    CDEF(loadoptionsRFU, 0, none),		\
    CDEF(detachoptionsRFU, 0, none),
  
#define callRFUMethods				\
  CALLDEF(SIMDmessages, 1),			\
    CALLDEF(DebugCall, 0),			\
    CALLDEF(Chol, 1),				\
    CALLDEF(debuggingLevel, 0),			\
    CALLDEF(scalarR, 3),			\
    CALLDEF(SolvePosDefR, 3),			\
    CALLDEF(struve, 4),				\
    CALLDEF(besselk_simd, 2),			\
    CALLDEF(I0ML0, 1),				\
    CALLDEF(gaussr, 2),				\
    CALLDEF(WMr, 4),				\
    CALLDEF(logWM2r, 4),			\
    CALLDEF(sortX, 4),				\
    CALLDEF(orderX, 4),				\
    CALLDEF(DivByRow, 2),			\
    CALLDEF(colMaxs, 1),			\
    CALLDEF(quadratic, 2),			\
    CALLDEF(dotXV, 2),				\
    CALLDEF(rowMeansX, 2),			\
    CALLDEF(rowProd, 1),			\
    CALLDEF(dbinorm, 2),			\
    CALLDEF(chol2mv, 2),			\
    CALLDEF(tcholRHS, 2),			\
    CALLDEF(crossprodX, 3),			\
    CALLDEF(getPackagesToBeInstalled, 1),	\
    CALLDEF(isGPUavailable,0),			\
    CALLDEF(isNEONavailable,0),			\
    CALLDEF(isX86_64,0),			\
    CALLDEF(gpu_info,1),			\
    CALLDEF(instruction_set, 3),		\
    CALLDEF(testStrassen, 4),			\
    CALLDEF(Update_utilsoption, 0),		\


#define  extRFUMethods				\
  EXTDEF(RFoptions, -1), 


#define RFU_CALLABLE				\
  CALLABLE(startRFU);				\
    CALLABLE(del_utilsoption);			\
    CALLABLE(get_utilsoption);			\
    CALLABLE(get_utils_basic);			\
    CALLABLE(push_utilsoption);			\
    CALLABLE(params_utilsoption);		\
    CALLABLE(update_utilsoption);		\
						\
    CALLABLE(solve_DELETE);			\
    CALLABLE(solve_NULL);			\
						\
    CALLABLE(SolvePosDef);			\
    CALLABLE(SolvePosDefSp);			\
    CALLABLE(SqrtPosDefFree);			\
    CALLABLE(xCinvXdet);			\
    CALLABLE(xCinvYdet);			\
    CALLABLE(DetPosDefsp);			\
    CALLABLE(InvertMatrix);			\
    CALLABLE(cholesky);				\
    CALLABLE(DetPosDef);			\
    CALLABLE(Is_positive_definite);		\
    						\
    CALLABLE(sqrtRHS);				\
    CALLABLE(chol2inv);				\
						\
    CALLABLE(StruveH);				\
    CALLABLE(StruveL);				\
    CALLABLE(I0mL0);				\
						\
    CALLABLE(WM);				\
    CALLABLE(DWM);				\
    CALLABLE(DDWM);				\
    CALLABLE(D3WM);				\
    CALLABLE(D4WM);				\
    CALLABLE(logWM2);				\
						\
    CALLABLE(Gauss);				\
    CALLABLE(DGauss);				\
    CALLABLE(DDGauss);				\
    CALLABLE(D3Gauss);				\
    CALLABLE(D4Gauss);				\
    CALLABLE(logGauss);				\
						\
    CALLABLE(attachRFUoptions);			\
    CALLABLE(detachRFUoptions);			\
    /*  CALLABLE(linkRFUoptions); */		\
    CALLABLE(RFUoptions);			\
    CALLABLE(attachSetNGet);			\
    CALLABLE(getoptionsRFU);			\
    CALLABLE(setoptionsRFU);			\
						\
    CALLABLE(ordering);				\
    CALLABLE(orderingL);			\
    CALLABLE(orderingInt);			\
    CALLABLE(orderingLong);			\
    CALLABLE(sorting);				\
    CALLABLE(sortingL);				\
    CALLABLE(sortingInt);			\
    CALLABLE(sortingLong);			\
    CALLABLE(scalarX);				\
    CALLABLE(scalarXint);			\
    /*  CALLABLE(scalarInt); */			\
						\
    CALLABLE(pid);				\
    CALLABLE(parallel);				\
    CALLABLE(sleepMicro); /* problem? */	\

  
 /* 
!!!!! HIER NIE EIN S E X P OBJEKT ZURUECKGEBEN  !!!!  
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


#endif


