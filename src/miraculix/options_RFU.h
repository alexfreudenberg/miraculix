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


  
#ifndef rfutils_options_H
#define rfutils_options_H 1

#include "AutoRandomFieldsUtilsLocal.h"
#include "Basic_RandomFieldsUtils.h"


#define SOLVE 0
#define MATRIXSQRT 1
#define DETERMINANT 2

#define SOLVE_METHODS 3


#if defined SCHLATHERS_MACHINE
#define DETERM_LAMODE false
#else
#define DETERM_LAMODE true
#endif

#define basicN 12
#define GENERAL_EXACTNESS 11
// IMPORTANT: all names of basic must be have least 3 letters !!!
typedef // benoetigt
struct basic_options {
  int  
  Rprintlevel, Cprintlevel, seed, cores,
    efficient,//allow for different level later on
    dummy0[4];
  usr_bool exactness;
  bool skipchecks, helpinfo, asList /* hidden:verbose */,
    bigendian, warn_parallel, 
    dummy6, dummy7;
  int dummy8[8];
  double NA, NaN;
} basic_options;
#define basic_START \
  { R_PRINTLEVEL, C_PRINTLEVEL,		\
      NA_INTEGER, INITCORES,				\
      true, /* different levels later on */		\
    {0, 0, 0, 0},					\
      Nan,						\
      false, true, true, false,				\
      false, false, false,			\
	{0,0,0,0, 0,0,0,0},						\
  1.79769313486231570E+308, 1.79769313486231570E+308 /* reset when started */ \
  }

#define installNrunN 9
#define MAX_GPU_DEVICES 16
typedef // benoetigt
#define INSTALL_RUN_WARN_OPTION 1
struct installNrun_options {
  int  
   warn_unknown_option, LaMaxTakeIntern,
    gpu_devices[MAX_GPU_DEVICES], Ngpu_devices, maxStreams,
    dummy0[4];
  install_modes install, dummy1;
  la_modes la_usr, la_mode, dummy2;
  usr_bool mem_is_aligned;
  bool installPackages, determineLAmode,
    kahanCorrection, dummy3,
    dummy4, dummy5, dummy6, dummy7;
  int dummy8[8];
} installNrun_options;
#define installNrun_START \
  { WARN_UNKNOWN_OPTION_ALL, MAXINT,				\
      {0}, 0, 0,						\
      {0, 0, 0, 0},						\
      INSTALL_DEFAULT, Inone,					\
      LA_AUTO, LA_R, LA_AUTO, /*LA_R  */			\
      MEMisALIGNED,						\
      true, false, DETERM_LAMODE,				\
      false,							\
      false, false, false, false,				\
      {0,0,0,0, 0,0,0,0}					\
  }


#define SOLVE_SVD_TOL 3
#define solveN 26
typedef // benoetigt
struct solve_options {
  usr_bool sparse, pivot_check, dummy0, dummy1;
  bool det_as_log, pivot_partialdet, pseudoinverse, dummy2, dummy3;
  double spam_tol, spam_min_p[2], svd_tol, eigen2zero, pivot_relerror,
    max_deviation, max_reldeviation, StrassenFactor,
    dummy4[4];
  InversionMethod Methods[SOLVE_METHODS], dummy5;
  int spam_min_n[2], spam_sample_n, spam_factor, pivotsparse, max_chol,
    max_svd,
    pivot, // obsolete
     actual_size,
    *pivot_idx, n_pivot_idx,//permutation; phys+logi laenge
    tinysize, AtAmode, AtAnrow, AtAncol, StrassenMin,
    dummy6[6];
  //  bool tmp_delete;
  pivot_modes actual_pivot,pivot_mode, dummy7;
  int dummy8[10];
 } solve_options;
#ifdef SCHLATHERS_MACHINE
#define svd_tol_start 1e-08
#else
#define svd_tol_start 0
#endif
#define solve_START							\
  False, False, False, False,						\
    true, false, false,	false, false,					\
    2.220446e-16, {0.8, 0.9}, svd_tol_start, 1e-12, 1e-11,		\
    1e-10, 1e-10, 1.5,							\
    {0.0, 0.0, 0.0, 0.0},						\
 {NoInversionMethod,  NoFurtherInversionMethod},NoInversionMethod,	\
    {400, 10000}, 500, 4294967, PIVOTSPARSE_MMD, 16384,			\
    10000,  /* never change -- see RFoptions.Rd */			\
    PIVOT_NONE, /* obsolete */						\
    0, NULL, 0, 3, 1, 0, 0, 512,					\
    {0,0,0,0,0, 0,},						\
   PIVOT_UNDEFINED, PIVOT_AUTO, PIVOT_UNDEFINED, /* PIVOT_NONE */	\
    {0,0,0,0,0, 0,0,0,0,0}

typedef // benoetigt
struct dummy_options {
  int dummy[30];
} dummy_options;

#define dummy_START \
  { {0,0,0,0, 0,0,0,0, 0,0,			\
	0,0,0,0, 0,0,0,0, 0,0,			\
	0,0,0,0, 0,0,0,0, 0,0}			\
  }


typedef // benoetigt
struct utilsoption_type{
  basic_options basic;
  installNrun_options installNrun;
  solve_options solve;
  dummy_options dummy;
} utilsoption_type;



#define ADD(ELT) SET_VECTOR_ELT(sublist, k++, ELT)
#define ADDCHAR(ELT) x[0] = ELT; ADD(ScalarString(mkChar(x)))

//int own_chol_up_to(int size, int maxtime);
//int own_chol_up_to();
void SetLaMode();
void SetLaMode(la_modes, int cores);
void resetInstalled();


#endif
