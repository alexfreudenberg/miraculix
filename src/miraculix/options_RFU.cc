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

#if defined compatibility_to_R_h
#include "kleinkram.h"
#include "zzz_RFU.h"
#include "xport_import_RFU.h"
#include "options_RFU.h"
#include "extern_RFU.h"


#define PLverbose 2



//#if defined(unix) || defined(__unix__) || defined(__unix)
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
#else
int numCPU = MAXINT;
#endif

void setoptionsRFU(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		   bool isList, utilsoption_type *options) {
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
    case 3: {
      int cores = POSINT;
      //      if (cores > 60)
      //	ERR("performance limit of omp is about 60 cores.");
      if (cores > numCPU) {
	WARN1("number of 'cores' is set to %d", numCPU);
	gp->cores = numCPU;
      } else gp->cores = cores;
#ifdef DO_PARALLEL
#else
      if (gp->cores != 1) {
	gp->cores = 1;
	PRINTF("The system does not allow for OpenMP: the value 1 is kept for 'cores'.");
      }    
#endif
      if (options->installNrun.determineLAmode) SetLaMode();
    }
      break;
    case 4: gp->skipchecks = LOGI;    break;
    case 5: gp->asList = LOGI; break;
    case 6 : if (!isList) {
	PL = gp->Cprintlevel = gp->Rprintlevel = 1 + (LOGI * (PLverbose - 1));
      }
      break;
    case 7: gp->helpinfo = LOGI; break;
    case 8 : gp->efficient = LOGI; break;
    case 9: break; // bigendian; read only
    case 10 : gp->warn_parallel = LOGI;  break;
    case GENERAL_EXACTNESS: gp->exactness = USRLOG; break; 
 
    default: BUG;
    }}
    break;
    
  case 1: {
    installNrun_options *gp = &(options->installNrun);
    switch(j) {
    case 0: gp->kahanCorrection = LOGI; break;
    case INSTALL_RUN_WARN_OPTION: gp->warn_unknown_option = INT;
      break;
    case 2: {
      int neu;
      if (TYPEOF(el) == STRSXP) {
	neu = GetName(el, name, LA_NAMES, LA_LAST + 1, gp->la_usr);
	//if (neu == LA_QUERY) {
	// if (!local) // OK
	//   PRINTF("internal mode = '%.10s'\nmax size for internal Cholesky = %d\nmax size for tiny size algorithm = %d\n", LA_NAMES[gp->la_mode], gp->LaMaxTakeIntern,  options->solve.tinysize);
	//}
      } else {
	neu = POS0INT;
	if (neu > LA_LAST) ERR0("wrong value for 'la_mode'");
      }				  
      if (neu != LA_QUERY) {
#ifdef USEGPU
#else
	if (neu > LA_GPU)
	  ERR1("In order to use '%.20s', install the package with apropriate compiler options.",
	       LA_NAMES[LA_GPU]);
#endif
	SetLaMode((la_modes) neu,
		  GreaterZero(options->basic.cores));  // koennen noch fehler auftreten
	gp->la_usr = (la_modes) neu;
      }
    }
      break;     
    case 3 : {
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
    case 4 :
      // gp->installPackages != LOGI;  
      break;
    case 5 : gp->determineLAmode = LOGI;  break;
    case 6 :
      // gp->mem_is_aligned = LOGI;
      break;
    case 7 :  Integer(el, name, gp->gpu_devices, MAX_GPU_DEVICES) ;
      gp-> Ngpu_devices = MIN(length(el), MAX_GPU_DEVICES);
      break;
     case 8 :
       gp->maxStreams = POS0INT;
      break;
   default: BUG;
    }}
    break;
    
  case 2: {
    //   printf("name = %.50s %d\n", name, j);
    
    solve_options *so = &(options->solve);
    switch(j) {
    case 0: so->sparse = USRLOG;
      if (so->sparse != False) {
	so->sparse = False;
	ERR0("'spam' is currently disabled.")
	  }      
      break; // USRLOGRELAXED??
    case 1: so->spam_tol = POS0NUM; break;
    case 2: Real(el, name, so->spam_min_p, 2);
      for (int u=0; u<=1; u++)
	so->spam_min_p[u] = so->spam_min_p[u] < 0 ? 0
  	: so->spam_min_p[u] > 1.0 ? 1.0 : so->spam_min_p[u];
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
    case 12: 
      so->pivot_mode = (pivot_modes) GetName(el, name, PIVOT_NAMES,
					     PIVOT_LAST + 1, so->pivot_mode);
      break;    
    case 13: if (!isList) {
      int len = length(el);
      if (len == 0) {
	if (so->n_pivot_idx > 0) { FREE(so->pivot_idx); }
      } else {
	if (so->n_pivot_idx != len) {
	  FREE(so->pivot_idx);
	  so->pivot_idx = (int*) MALLOC(len * sizeof(int));
	}
	for (int L=0; L<len; L++) so->pivot_idx[L] = Integer(el, name, L);
      }
      so->n_pivot_idx = len;     
    }
      break; 
    case 14: so->pivot_relerror = POS0NUM; break;    
    case 15: so->max_deviation = POSNUM; break;    
    case 16: so->max_reldeviation = POS0NUM; break;    
    case 17: so->det_as_log = LOGI; break;    
    case 18: so->actual_size = POS0NUM; break;    
    case 19: so->pivot_check = USRLOG; break;    
    case 20: so->pseudoinverse = LOGI; break;
    case 21: so->AtAmode = POS0INT; break;
    case 22: so->AtAnrow = POS0INT; break;
    case 23: so->AtAncol = POS0INT; break;
    case 24: so->StrassenMin = POSINT; break;
    case 25: so->StrassenFactor = POSNUM; break;
    default: BUG;
    }}
    break;
    
    default: BUG;
  }
}

void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		       bool isList, bool local) {
  if (!local && parallel()) 
    ERR1("Option '%.25s' can be set only through 'RFoptions' at global level",
	 allOptions[i][j]);
  setoptionsRFU(i, j, el, name, isList,  WhichOptionList(local));
}



void getoptionsRFU(SEXP sublist, int i, utilsoption_type *options) {  
   int  k = 0;
 //printf("OK %d\n", i);    
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
    ADD(ScalarLogical(p->helpinfo));    
    ADD(ScalarLogical(p->efficient));
    ADD(ScalarLogical(p->bigendian));
    ADD(ScalarLogical(p->warn_parallel));
    ADD(ExtendedBooleanUsr(p->exactness));    

  }
    break;
 
  case 1 : {
    installNrun_options *p = &(options->installNrun);
    ADD(ScalarLogical(p->kahanCorrection));   
    ADD(ScalarInteger(p->warn_unknown_option));    
    ADD(ScalarString(mkChar(LA_NAMES[p->la_usr])));
    ADD(ScalarString(mkChar(INSTALL_NAMES[p->install])));
    ADD(ScalarLogical(p->installPackages));
    ADD(ScalarLogical(p->determineLAmode));
    ADD(ScalarLogical(p->mem_is_aligned));
    SET_VECTOR_ELT(sublist, k++, Int(p->gpu_devices, p->Ngpu_devices)); 
    ADD(ScalarInteger(p->maxStreams));    
  }
    break;
 
  case 2 : {
    solve_options *p = &(options->solve);
    //    printf("sparse user interface %d %d %d\n", p->sparse, NA_LOGICAL, NA_INTEGER);
    ADD(ExtendedBooleanUsr(p->sparse));
    //
      ADD(ScalarReal(p->spam_tol));    
    SET_VECTOR_ELT(sublist, k++, Num(p->spam_min_p, 2));
   ADD(ScalarReal(p->svd_tol)); 
    ADD(ScalarReal(p->eigen2zero));
      
   SET_VECTOR_ELT(sublist, k++,
		   String((int*) p->Methods, InversionNames, SOLVE_METHODS,
			  (int) NoFurtherInversionMethod));	
   //   printf("A\n");
   SET_VECTOR_ELT(sublist, k++, Int(p->spam_min_n, 2));
   ADD(ScalarInteger(p->spam_sample_n));    
    ADD(ScalarInteger(p->spam_factor));    
    ADD(ScalarInteger(p->pivotsparse));    
    ADD(ScalarInteger(p->max_chol));    
    ADD(ScalarInteger(p->max_svd)); 
    ADD(ScalarString(mkChar(PIVOT_NAMES[p->pivot_mode])));
    //if (true)
    SET_VECTOR_ELT(sublist, k++, Int(p->pivot_idx, p->n_pivot_idx));
    //  else ADD(ScalarInteger(NA_INTEGER));
    //    ADD(ScalarLogical(p->tmp_delete));
     ADD(ScalarReal(p->pivot_relerror));    
    ADD(ScalarReal(p->max_deviation));    
    ADD(ScalarReal(p->max_reldeviation));    
    ADD(ScalarLogical(p->det_as_log));
    ADD(ScalarInteger(p->actual_size));
    ADD(ExtendedBooleanUsr(p->pivot_check));
    ADD(ScalarLogical(p->pseudoinverse));
    ADD(ScalarInteger(p->AtAmode));
    ADD(ScalarInteger(p->AtAnrow));
    ADD(ScalarInteger(p->AtAncol));
   ADD(ScalarInteger(p->StrassenMin));
   ADD(ScalarReal(p->StrassenFactor));
  }
    break;
  default : BUG;
  }
  //  printf("EE A\n");
}

void getoptions(SEXP sublist, int i, bool local) {  
  getoptionsRFU(sublist, i, WhichOptionList(local));
}


#endif
