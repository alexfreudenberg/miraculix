/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 -- 2017  Martin Schlather

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

#include "Basic_utils.h" // must be before anything else
#ifdef DO_PARALLEL
#include <omp.h>
#endif
#include <unistd.h>
#include "options.h"
#include "General_utils.h"
#include "kleinkram.h"
#include "zzz_RandomFieldsUtils.h"
#include "xport_import.h"

#ifdef TIME_AVAILABLE
#include <time.h>
#endif

#define PLverbose 2

// IMPORTANT: all names of general must have at least 3 letters !!!
const char *basic[basicN] =
  { "printlevel","cPrintlevel", "seed", "cores", "skipchecks",
    "asList", "verbose", "kahanCorrection", "helpinfo", "warn_unknown_option",
    "la_mode", "warn_parallel", "install", "installPackages"};

const char * solve[solveN] = 
  { "use_spam", "spam_tol", "spam_min_p", "svdtol", "eigen2zero",
    "solve_method", "spam_min_n", "spam_sample_n", "spam_factor", "spam_pivot",
    "max_chol", "max_svd", "pivot",
    "pivot_idx", // dynamic parameter
    "pivot_relerror", "pivot_max_deviation", "pivot_max_reldeviation",
    "det_as_log", "pivot_actual_size", "pivot_check", "pseudoinverse"
    //, "tmp_delete"
  };


const char * prefixlist[prefixN] = {"basic", "solve"};
const char **all[prefixN] = {basic, solve};
int allN[prefixN] = {basicN, solveN};


utilsoption_type OPTIONS = { // OK
  basic_START,
  { solve_START }
};


//#if defined(unix) || defined(__unix__) || defined(__unix)
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
#else
int numCPU = MAXINT;
#endif




int doPosDef(double *M0, int size, bool posdef,
	     double *rhs0, int rhs_cols, double *result,
	     double *logdet, int calculate, solve_storage *Pt,
	     solve_options *sp);


void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		       bool isList, bool local) {
  //  printf("i=%d j=%d\n", i, j);
  if (local || parallel())
    ERR1("Option '%.25s' can be set only through 'RFoptions' at global level",
	 all[i][j])
  utilsoption_type *options = &OPTIONS; 
  switch(i) {
  case 0: {// general
    basic_options *gp = &(options->basic);
    switch(j) {
    case 0: {  // general options
      int threshold = 1000; //PL_ERRORS;
      gp->Rprintlevel = INT;
      if (gp->Rprintlevel > threshold) gp->Rprintlevel = threshold;
      PL = gp->Cprintlevel = gp->Rprintlevel + PLoffset;
    }
      break;
    case 1: PL = gp->Cprintlevel = INT + PLoffset ;
      break;
    case 2: gp->seed = Integer(el, name, 0, true); break;
    case 3: gp->cores = POSINT;
      if (gp->cores > numCPU) {
	WARN1("number of 'cores' is set to %d", numCPU);
	gp->cores = numCPU;
      }
#ifdef DO_PARALLEL
#else
      if (gp->cores != 1) {
	gp->cores = 1;
	PRINTF("The system does not allow for OpenMP: the value 1 is kept for 'cores'.");
      }      
#endif
      SetLaMode();
      CORES = gp->cores;
      break;
    case 4: gp->skipchecks = LOGI;    break;
    case 5: gp->asList = LOGI; break;
    case 6 : if (!isList) {
	PL = gp->Cprintlevel = gp->Rprintlevel = 1 + (LOGI * (PLverbose - 1));
      }
      break;
    case 7: gp->kahanCorrection = LOGI; break;
    case 8: gp->helpinfo = LOGI; break;
    case BASIC_WARN_OPTION: gp->warn_unknown_option = INT;
	break;
    case 10: {
      int neu;
      if (TYPEOF(el) == STRSXP)
	neu = GetName(el, name, LA_NAMES, LA_LAST + 1, gp->la_usr);
      else {
	neu = POS0INT;
	if (neu > LA_LAST) ERR("wrong value for 'la_mode'");
      }				  
#ifdef USEGPU
#else
      if (neu > LA_GPU)
	ERR1("In order to use '.20s', install the package with apropriate compiler options.",
	     LA_NAMES[LA_GPU]);
#endif
      SetLaMode((la_modes) neu);  // koennen noch fehler auftreten
      //      printf("simu ende\n");
     gp->la_usr = (la_modes) neu;
    }
      break;
     
    case 11 : gp->warn_parallel = LOGI;  break;
    case 12 : {
      install_modes old =  gp->install;
      gp->install = (install_modes )
	GetName(el, name, INSTALL_NAMES, INSTALL_LAST + 1, Iask);
      if (gp->install == Inone) gp->installPackages = false;
      else if (old != gp->install) {
	gp->installPackages = true;
	resetInstalled();
      }
    }
      break;
    case 13 :
      // gp->installPackages != LOGI;  
      break;
    default: BUG;
    }}
    break;
 
 case 1: {
   //   printf("name = %.50s %d\n", name, j);
   
    solve_options *so = &(options->solve);
    switch(j) {
    case 0: so->sparse = USRLOG;
      if (so->sparse != False) {
	so->sparse = False;
	ERR("'spam' is currently disabled.")
      }      
      break; // USRLOGRELAXED??
    case 1: so->spam_tol = POS0NUM; break;      
    case 2: Real(el, name, so->spam_min_p, 2);
      for (int u=0; u<=1; u++)
	so->spam_min_p[u] = so->spam_min_p[u] < 0 ? 0 :
			    so->spam_min_p[u] > 1.0 ? 1.0 : so->spam_min_p[u];
      break;      
    case SOLVE_SVD_TOL: so->svd_tol = POS0NUM; break;        
    case 4: so->eigen2zero = POS0NUM; break;        
    case 5:
      GetName(el, name, InversionNames, nr_user_InversionMethods,
		    (int) NoInversionMethod, (int) NoFurtherInversionMethod, 
		    (int *)so->Methods, SOLVE_METHODS);
      break;
    case 6: Integer(el, name, so->spam_min_n, 2);      
      break; 
    case 7: so->spam_sample_n = POSINT; break;      
    case 8: so->spam_factor = POSINT; break;      
    case 9: so->pivotsparse = POSINT; 
      if (so->pivotsparse > 2) so->pivotsparse = PIVOT_NONE;
      break;    
    case 10:
      //      printf("max chol = %d\n",  so->max_chol);
      so->max_chol = POSINT;
      //      printf("X max chol = %d\n",  so->max_chol);
       break;      
    case 11: so->max_svd = POS0INT; break;    
      //    case 11: so->tmp_delete = LOGI; break;    
    case 12: {
      int neu;
      if (TYPEOF(el) == STRSXP)
	neu = GetName(el, name, PIVOT_NAMES, PIVOT_LAST + 1, so->pivot_mode);
      else {
	neu = POS0INT;
	if (neu > PIVOT_LAST) ERR("wrong value for 'privot'");
      }
      so->pivot_mode = (pivot_modes) neu;
    }
      break;    
    case 13: if (!isList) {
      int len = length(el);
      if (len == 0) {
	if (so->pivot_idx_n > 0) FREE(so->pivot_idx);
      } else {
	if (so->pivot_idx_n != len) {
	  FREE(so->pivot_idx);
	  so->pivot_idx = (int*) MALLOC(len * sizeof(int));
	}
	for (int L=0; L<len; L++) so->pivot_idx[L] = Integer(el, name, L);
      }
      so->pivot_idx_n = len;     
    }
      break; 
    case 14: so->pivot_relerror = POS0NUM; break;    
    case 15: so->max_deviation = POSNUM; break;    
    case 16: so->max_reldeviation = POS0NUM; break;    
    case 17: so->det_as_log = LOGI; break;    
    case 18: so->actual_size = POS0NUM; break;    
    case 19: so->pivot_check = USRLOG; break;    
    case 20: so->pseudoinverse = LOGI; break;    
   default: BUG;
    }}
    break;
    
  default: BUG;
  }
}


void getoptions(SEXP sublist, int i,
		       bool VARIABLE_IS_NOT_USED local) {
  int  k = 0;
  utilsoption_type *options = &OPTIONS; 
  switch(i) {
  case 0 : {
    //    printf("OK %d\n", i);
    basic_options *p = &(options->basic);
    ADD(ScalarInteger(p->Rprintlevel));    
    ADD(ScalarInteger(p->Cprintlevel - PLoffset));
    ADD(ScalarInteger(p->seed));    
    ADD(ScalarInteger(p->cores));    
    ADD(ScalarLogical(p->skipchecks));
    ADD(ScalarLogical(p->asList));   
    ADD(ScalarLogical(p->Rprintlevel >= PLverbose));
    ADD(ScalarLogical(p->kahanCorrection));   
    ADD(ScalarLogical(p->helpinfo));    
    ADD(ScalarInteger(p->warn_unknown_option));    
    ADD(ScalarString(mkChar(LA_NAMES[p->la_usr])));
    ADD(ScalarLogical(p->warn_parallel));
    ADD(ScalarString(mkChar(INSTALL_NAMES[p->install])));
     ADD(ScalarLogical(p->installPackages));
  }
    break;
  
  case 1 : {
    solve_options *p = &(options->solve);
    //    printf("sparse user interface %d %d %d\n", p->sparse, NA_LOGICAL, NA_INTEGER);
    ADD(ExtendedBooleanUsr(p->sparse));
    //    printf("A\n");
    ADD(ScalarReal(p->spam_tol));    
    SET_VECTOR_ELT(sublist, k++, Num(p->spam_min_p, 2));
    ADD(ScalarReal(p->svd_tol)); 
    ADD(ScalarReal(p->eigen2zero));    
    SET_VECTOR_ELT(sublist, k++,
		   String((int*) p->Methods, InversionNames, SOLVE_METHODS,
			  (int) NoFurtherInversionMethod));	
    SET_VECTOR_ELT(sublist, k++, Int(p->spam_min_n, 2));
    ADD(ScalarInteger(p->spam_sample_n));    
    ADD(ScalarInteger(p->spam_factor));    
    ADD(ScalarInteger(p->pivotsparse));    
    ADD(ScalarInteger(p->max_chol));    
    ADD(ScalarInteger(p->max_svd)); 
    ADD(ScalarString(mkChar(PIVOT_NAMES[p->pivot_mode])));
    //if (true)
    SET_VECTOR_ELT(sublist, k++, Int(p->pivot_idx, p->pivot_idx_n));
    //  else ADD(ScalarInteger(NA_INTEGER));
    //    ADD(ScalarLogical(p->tmp_delete));
    ADD(ScalarReal(p->pivot_relerror));    
    ADD(ScalarReal(p->max_deviation));    
    ADD(ScalarReal(p->max_reldeviation));    
    ADD(ScalarLogical(p->det_as_log));
    ADD(ScalarInteger(p->actual_size));
    ADD(ExtendedBooleanUsr(p->pivot_check));
    ADD(ScalarLogical(p->pseudoinverse));
  }
    break;
  default : BUG;
  }
  //printf("EE A\n");
}

void utilsoption_NULL (utilsoption_type *S) {
  assert(solveN == 21 && basicN == 12 && prefixN==2);
  MEMCOPY(S, &OPTIONS, sizeof(utilsoption_type));
  S->solve.pivot_idx = NULL;
  S->solve.pivot_idx_n = 0;
}

void utilsoption_DELETE(utilsoption_type *S) {
  FREE(S->solve.pivot_idx);
  S->solve.pivot_idx_n = 0;
}



void deloptions(bool VARIABLE_IS_NOT_USED local) {
#ifdef DO_PARALLEL
  if (local) RFERROR("'pivot_idx' cannot be freed on a local level");
#endif  
  utilsoption_type *options = &OPTIONS;
  FREE(options->solve.pivot_idx);
}



#define FASTER 1.3 // 1.3 as for typical application of likelihood,
// the determinant calculation in RandomFieldsUtils is for free. Somehow a
// balance
int own_chol_up_to(int size, int maxtime) {
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
  while (true) {
    int sizeP1 = size + 1,
      sizesq = size * size,
      loops = size > 64 ? 1 : 16384 / ((size + 8) * (size+8));
    double *M = (double*) MALLOC(sizesq * sizeof(double));
    for (int j=0; j<=1; j++) {
      SetLaMode(j == 0 ? LA_INTERN : LA_R);
      clock_t start = clock();
      for (int k=0; k<loops; k++) {      
	MEMSET(M, 0, sizesq * sizeof(double));
	for (int i=0; i<sizesq; i+=sizeP1) M[i] = 1.0;
	if (size > 1) M[1] = M[size] = 1e-5;
	doPosDef(M, size, true, NULL, 0, NULL, NULL, MATRIXSQRT, &pt, &sp);   
      }
      delta[j] = (Long) clock() - start; 
      if (delta[j] < 0) delta[j] += MAXINT; // manual: 32bit repres. of clock
    }
    FREE(M);
    if (PL > 2)
      PRINTF("Cholesky decomposition for a %d x %d matrix needs %ld and %ld [mu s] on R and facile code on %d cores (#%d), respectively.\n", size, size, delta[1], delta[0], CORES, loops);
    if (delta[0] > maxtime || delta[1] > maxtime ||
	delta[0] >= FASTER * delta[1]) break;
    old_quotient = (double) delta[0] / delta[1];
    size *= 2;
  }
  solve_DELETE0(&pt);
  double new_quotient = (double) delta[0] / delta[1];
  if (new_quotient < FASTER) return MAXINT;
  if (size <= 0) return(0);
  if (!R_FINITE(old_quotient)) {
    int compare = own_chol_up_to(size / 2, 0);
   return compare == MAXINT ? size : compare;
  }
  double x0 = 0.5 * size * (1.0 + (FASTER - old_quotient) /
			    (new_quotient - old_quotient)); //lin interpol
  assert(x0 >= 0.5 * size && x0 <= size);
   int compare = own_chol_up_to((int) x0, 0);
  //  printf("%f %f %f %f %d %d\n", x0,FASTER,  old_quotient, new_quotient, size, compare);
  return (int) (compare == MAXINT ?  x0 : 0.5 * size);
#else 
  ERR("option 'LA_AUTO' is available only on linux systems");
  return 0;
#endif
}
int own_chol_up_to() {
  own_chol_up_to(256, 0); //warm up for some BLAS implementatioan
  //  CORES = GL OBAL.basic.cores = 4;
  //  own_chol_up_to(256, 50000);
  //  own_chol_up_to(8, 50000);
  return own_chol_up_to(256, 50000);
}


  
void SetLaMode(la_modes usr_mode) {
  //  printf("stlamode=%d\n", usr_mode);
  utilsoption_type *global = &OPTIONS; 
  la_modes la_mode = usr_mode;
  global->solve.tinysize = global->basic.LaMaxTakeIntern = -1;
#define TINY_SIZE_MAX 3  
  if (la_mode == LA_INTERN) {
    global->solve.tinysize = TINY_SIZE_MAX;
    global->basic.LaMaxTakeIntern = MAXINT;
  } else if (la_mode == LA_AUTO) {  
    la_mode = HAS_GPU ? LA_GPU : LA_R;
#if defined TIME_AVAILABLE
#  ifdef SCHLATHERS_MACHINE
#else    
    int PLalt = PL;
    PL = 0;
#  endif  
    global->basic.LaMaxTakeIntern = own_chol_up_to();
    global->solve.tinysize = MIN(TINY_SIZE_MAX, global->basic.LaMaxTakeIntern);
    if (PL > 0)
      PRINTF("Limit size for facile Cholesky algorithm  = %d\n",
	     global->basic.LaMaxTakeIntern);
#  ifdef SCHLATHERS_MACHINE
#else   
    PL = PLalt;
#  endif
#endif
  }
  
  if ((la_mode == LA_GPU || la_mode == LA_R) &&
      global->solve.pivot_mode > PIVOT_AUTO)
    ERR("Pivotized Cholesky decomposition has not been implemented yet for GPU and the LAPACK library");
  
  global->basic.la_mode = la_mode;
}
void SetLaMode() {
  utilsoption_type *global = &OPTIONS; 
   SetLaMode(global->basic.la_usr);
}


