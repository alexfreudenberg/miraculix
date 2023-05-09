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


#include <unistd.h>
#include "parallel_simd.h"
#ifdef TIME_AVAILABLE
#include <time.h>
#endif




#include "Basic_RFUlocal.h" // must be before anything else
#include "kleinkram.h"
#include "zzz_RFU.h"
#include "options_RFU.h"
#include "extern_RFU.h"
#if defined compatibility_to_R_h     
#include "xport_import_RFU.h"
#endif

// IMPORTANT: all names of general must have at least 3 letters !!!
const char *basic[basicN] =
  { "printlevel","cPrintlevel", "seed", "cores",
    "skipchecks", "asList", "verbose", "helpinfo", "efficient",
    "bigendian","warn_parallel", "exactness"
  };

const char *installNrun[installNrunN] =
  { "kahanCorrection", "warn_unknown_option", "la_mode", 
    "install","installPackages", "determineLAmode", "mem_is_aligned",
    "gpuDevices", "maxStreams"
  };

const char * solve[solveN] = 
  { "use_spam", "spam_tol", "spam_min_p", "svdtol", "eigen2zero",
    "solve_method", "spam_min_n", "spam_sample_n", "spam_factor", "spam_pivot",
    "max_chol", "max_svd", "pivot",
    "pivot_idx", // dynamic parameter
    "pivot_relerror", "pivot_max_deviation", "pivot_max_reldeviation",
    "det_as_log", "pivot_actual_size", "pivot_check", "pseudoinverse",
    "AtAmode", "AtAnrow", "AtAncol", "Strassen_min", "Strassen_factor"
    //, "tmp_delete"
  };


const char * prefixlist[prefixN] = {"basic", "installNrun", "solve"};
const char **allOptions[prefixN] = {basic, installNrun, solve};
int allOptionsN[prefixN] = {basicN, installNrunN, solveN};


utilsoption_type OPTIONS = { // OK
  basic_START,
  installNrun_START,
  { solve_START },
  dummy_START
};




int doPosDefIntern(double *M0, int size, bool posdef,
		   double *rhs0, Long rhs_cols, double *result,
		   double *logdet, int calculate, solve_storage *Pt,
		   solve_options *sp, int VARIABLE_IS_NOT_USED cores);

void params_utilsoption(int local, int *params) {
  utilsoption_type *from = WhichOptionList(local);
  params[PIVOT_IDX_N] = from->solve.n_pivot_idx;
}

void get_utils_basic(basic_options *opt, int local) {
  utilsoption_type *from = WhichOptionList(local);
  MEMCOPY(opt, from, sizeof(basic_options));
}

void get_utilsoption(utilsoption_type *S, int local) {
  assert(solveN == 21 && basicN == 9 && installNrunN == 10 && prefixN==3);
  utilsoption_type *from = WhichOptionList(local);
  assert(from->solve.n_pivot_idx!=0 xor from->solve.pivot_idx == NULL);
  assert(S->solve.n_pivot_idx!=0 xor S->solve.pivot_idx == NULL);
  if (S->solve.n_pivot_idx != from->solve.n_pivot_idx) BUG;
  int *save_idx = S->solve.pivot_idx,
    n = S->solve.n_pivot_idx < from->solve.n_pivot_idx ?
    S->solve.n_pivot_idx : from->solve.n_pivot_idx;
  MEMCOPY(S, from, sizeof(utilsoption_type)); // OK
  S->solve.pivot_idx = save_idx;
  if (n > 0) {
    MEMCOPY(S->solve.pivot_idx, from->solve.pivot_idx,
	    sizeof(*(S->solve.pivot_idx)) * n);
  }
}

void push_utilsoption(utilsoption_type *S, int local) {
  utilsoption_type *to = WhichOptionList(local);
  assert(to->solve.n_pivot_idx!=0 xor to->solve.pivot_idx == NULL);
  assert(S->solve.n_pivot_idx!=0 xor S->solve.pivot_idx == NULL);
  int *save_idx = to->solve.pivot_idx;
  if (to->solve.n_pivot_idx != S->solve.n_pivot_idx) {
    FREE(to->solve.pivot_idx);
    to->solve.pivot_idx = (int*) MALLOC(S->solve.n_pivot_idx * sizeof(int));
    save_idx = to->solve.pivot_idx;
  }
  MEMCOPY(to, S, sizeof(utilsoption_type)); // OK
  to->solve.pivot_idx = save_idx;
  if (S->solve.n_pivot_idx > 0) {
    MEMCOPY(to->solve.pivot_idx, S->solve.pivot_idx,
	    sizeof(*(to->solve.pivot_idx)) * S->solve.n_pivot_idx);
  }
}

void del_utilsoption(utilsoption_type *S) {
  FREE(S->solve.pivot_idx);
  S->solve.n_pivot_idx = 0;
}

void update_utilsoption() {
  push_utilsoption(&OPTIONS, true);
}

#define FASTER 1.3 // 1.3 as for typical application of likelihood,
// the determinant calculation in RandomFieldsUtils is for free. Somehow a
// balance
int own_chol_up_to(int size, int maxtime, int VARIABLE_IS_NOT_USED cores) {
#ifdef TIME_AVAILABLE
  if (size <= 0) return true;
  Long delta[2];
  solve_options sp;
  solve_storage pt;
  solve_NULL(&pt);
  MEMCOPY(&sp, &(OPTIONS.solve), sizeof(solve_options)); 
  sp.Methods[0] = Cholesky;					     
  sp.Methods[1] = NoFurtherInversionMethod;
  sp.pivot_mode = PIVOT_NONE;	  
  sp.sparse = False;	  
  double old_quotient = RF_NAN;
  // basic assumption is that R implementation's getting faster
  // for larger matrices
  //  printf("**** start\n");
  while (true) {
    //    printf("x \n");
    int sizeP1 = size + 1,
      sizesq = size * size,
      loops = size > 64 ? 1 : 16384 / ((size + 8) * (size+8)) / 4;
    double *M = (double*) MALLOC(sizesq * sizeof(double));
    for (int j=0; j<=1; j++) {
      //      printf("%d,%d\n", j, loops);
      SetLaMode(j == 0 ? LA_INTERN : LA_R, cores);
      clock_t start = clock();
      for (int k=0; k<loops; k++) {      
	MEMSET(M, 0, sizesq * sizeof(*M));
	for (int i=0; i<sizesq; i+=sizeP1) M[i] = 1.0;
	if (size > 1) M[1] = M[size] = 1e-5;
	//printf("size=%d\n", size);
	doPosDefIntern(M, size, true, NULL, 0, NULL, NULL, MATRIXSQRT, &pt, &sp,
		       cores);
	//printf("doen\n");
      }
      delta[j] = (Long) clock() - start; 
      if (delta[j] < 0) delta[j] += MAXINT; // manual: 32bit repres. of clock
    }
    FREE(M);
    if (PL > 2)
      PRINTF("Cholesky decomposition for a %d x %d matrix needs %ld and %ld [mu s] on R and facile code on %d cores (#%d), respectively.\n", size, size, delta[1], delta[0], cores, loops);

    //  printf("delta %d %d %d\n", delta[0], delta[1], maxtime);
    
    if (delta[0] > maxtime || delta[1] > maxtime ||
	(double) delta[0] >= FASTER * (double) delta[1]){
      solve_DELETE0(&pt);
      if ((maxtime > 0 &&
	   (delta[0] > 10 * maxtime || delta[1] > 10 * maxtime)) ||
	  delta[0] > 2 * delta[1] || delta[1] > 2 * delta[0]) {
	// seems to be time consuming. So stop.
	return (double) delta[0] < FASTER * (double) delta[1]
	  ? MAXINT : (size <= 0 ? 0 : size / 2);
      }
      break;
    }
    old_quotient = (double) delta[0] / (double) delta[1];
    size *= 2;
  }
  double new_quotient = (double) delta[0] / (double) delta[1];
  if (new_quotient < FASTER) return MAXINT;
  if (size <= 0) return(0);
  if (!R_FINITE(old_quotient)) {
    // printf("halfsize\n");
    int compare = own_chol_up_to(size / 2, 0, cores);
   return compare == MAXINT ? size : compare;
  }
  double x0 = 0.5 * size * (1.0 + (FASTER - old_quotient) /
			    (new_quotient - old_quotient)); //lin interpol
  assert(x0 >= 0.5 * size && x0 <= size);
  int compare = own_chol_up_to((int) x0, 0, cores);
  //  printf("%f %f %f %f %d %d\n", x0,FASTER,  old_quotient, new_quotient, size, compare);
  return (int) (compare == MAXINT ?  x0 : 0.5 * size);
#else 
  ERR0("option 'LA_AUTO' is available only on linux systems");
  return 0;
#endif
}
int own_chol_up_to(int VARIABLE_IS_NOT_USED cores) {
  own_chol_up_to(256, 0, cores); //warm up for some BLAS implementatioan
  //  own_chol_up_to(256, 50000);
  //  own_chol_up_to(8, 50000);
  return own_chol_up_to(256, 50000, cores);
}


  
void SetLaMode(la_modes usr_mode, int VARIABLE_IS_NOT_USED cores) {
  utilsoption_type *utils = &OPTIONS;
  la_modes la_mode = usr_mode;
  utils->solve.tinysize =
    utils->installNrun.LaMaxTakeIntern = -1;
#define TINY_SIZE_MAX 3  
  if (la_mode == LA_INTERN) {
    utils->solve.tinysize = TINY_SIZE_MAX;
    utils->installNrun.LaMaxTakeIntern = MAXINT;
  } else if (la_mode == LA_AUTO) {  
    la_mode = HAS_GPU ? LA_GPU : LA_R ;
#if defined TIME_AVAILABLE
#  ifdef SCHLATHERS_MACHINE
#else    
    int PLalt = PL;
    PL = 0;
#  endif  
    utils->installNrun.LaMaxTakeIntern = own_chol_up_to(cores);
    utils->solve.tinysize = MIN(TINY_SIZE_MAX, utils->installNrun.LaMaxTakeIntern);
    if (PL > 0)
      PRINTF("Limit size for facile Cholesky algorithm  = %d\n",
	     utils->installNrun.LaMaxTakeIntern);
#  ifdef SCHLATHERS_MACHINE
#else   
    PL = PLalt;
#  endif
#endif
  }
  
  if ((la_mode == LA_GPU || la_mode == LA_R) &&
      utils->solve.pivot_mode > PIVOT_AUTO)
    ERR0("Pivotized Cholesky decomposition has not been implemented yet for GPU and the LAPACK library");
  
  utils->installNrun.la_mode = la_mode;
}
void SetLaMode() {
  utilsoption_type *utils = &OPTIONS;
  SetLaMode(utils->installNrun.la_usr,
	    GreaterZero(utils->basic.cores));
}





utilsoption_type *WhichOptionList(bool local) {  
  if (local) {
 #if defined compatibility_to_R_h     
    KEY_type *KT = KEYT();
    if (KT == NULL) BUG;
    return &(KT->global_utils);
#else
    ERR0("'local' must be 'false' in standalone packages.")
#endif
  }
  return &OPTIONS;
}
