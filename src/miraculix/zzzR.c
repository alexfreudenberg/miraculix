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


#define NO_SSE2 1

#ifndef __cplusplus
#include <stdbool.h>
#endif

#include "def.h"
#include "Basic_RandomFieldsUtils.h"

#if defined compatibility_to_R_h

#include "miraculix.h"
// RFU INCLUDE #include "win_linux_aux.h"
// RFU INCLUDE #include "RandomFieldsUtils.h"
// RFU INCLUDE #include "utils_miraculix.h"
// RFU INCLUDE #include "zzz_RFU.h"


#define none 0

/*
void F 77 _SUB(specific_lines)(int *ntotal, int* nlines, int *whichlines,
			     int *datacols, int *idcols, int *datamat);
				  
#define FDEF(name, n, type) {#name,  (DL_FUNC) &F 77_SUB(name), n, type}
*/

/*
static R_NativePrimitiveArgType
  real2int2[] = {REALSXP, REALSXP, INTSXP, INTSXP},
  int5real2char[] = {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
		     STRSXP},
  int6[] = {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP}
  ;
*/




 /* 
static R_FortranMethodDef fortranMethods[] = {
   FDEF(gmatrix, 4, real2int2),
  FDEF(gmatrix_data, 8, int5real2char),
  FDEF(gmatrix_data_recoded, 8, int5real2char),
  FDEF(specific_lines, 6, int6),
   {NULL, NULL, 0, none}
};
 */
			    
#define CDEF(name, n, type) {#name, (DL_FUNC) &name, n, type}
static R_NativePrimitiveArgType
// RFU INCLUDE host_arg[] = { STRSXP, INTSXP},
  int_arg[] = { INTSXP };
//char2int7[] = { CHARSXP, CHARSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
//	      INTSXP, INTSXP}; 

static const R_CMethodDef cMethods[]  = {
  CDEF(loadoptions, 1, int_arg),
  CDEF(detachoptions, 0, none),
  // RFU INCLUDE cRFUMethods
  {NULL, NULL, 0, none}
};

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss adoption.h eingebunden sein
  
  CALLDEF(scan, 10),
  CALLDEF(sumscan, 10),
  CALLDEF(collect_scan, 10),
  CALLDEF(collect_scan2, 13),
  
  CALLDEF(windower, 9),
  CALLDEF(copyoptions, 0),

  CALLDEF(vector012matrix, 2),
  CALLDEF(matrixvector012, 2),
  CALLDEF(crossprod, 1),
  
  CALLDEF(codeOrigins, 1),
  CALLDEF(decodeOrigins, 2),
  CALLDEF(Transform, 4),
  CALLDEF(zeroGeno, 4),
  CALLDEF(rhaplomatrix, 4),

  CALLDEF(createSNPmatrix, 2),
  CALLDEF(fillSNPmatrix, 3),
  CALLDEF(vectorGeno, 2),
  CALLDEF(genoVector, 2),
  CALLDEF(computeSNPS, 8),
  CALLDEF(compute, 9),
  CALLDEF(allele_freq, 1),
  CALLDEF(Debug, 0),
  CALLDEF(StopDebug, 0),
  CALLDEF(solveRelMat, 5),
  CALLDEF(substract_centered, 1),
  CALLDEF(get_centered, 0),
  CALLDEF(transpose, 1),
  CALLDEF(crossprodInt, 3),
  CALLDEF(existsVariant, 3),
  CALLDEF(existsTiling, 3),
  CALLDEF(existsCrossprod, 1),
  CALLDEF(existsAllelefreq, 1),
  CALLDEF(existsCoding, 2),

  //  CALLDEF(),
 //  CALLDEF(codeSNPs, 2),

  // RFU INCLUDE callRFUMethods  
  {NULL, NULL, 0}
};


#define EXTDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_ExternalMethodDef extMethods[] = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  EXTDEF(MiraculixOptions, -1), 
// RFU INCLUDE extRFUMethods  
  {NULL, NULL, 0} 
};


#define CALLABLE(FCTN) R_RegisterCCallable("miraculix",#FCTN,(DL_FUNC) FCTN)
void R_init_miraculix(DllInfo  *dll) {
  R_registerRoutines(dll, cMethods, // .C
		     callMethods,
		     NULL, //fortranMethods, // .Fortran
		     extMethods); // ext
  R_useDynamicSymbols(dll, FALSE); // OK
  R_forceSymbols(dll, TRUE); // OK
}


#if defined __GNUC__ && __GCC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
//#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
void R_unload_miraculix(DllInfo *info) { }
#if defined __GNUC__ && __GCC__
//#pragma GCC diagnostic warning "-Wcast-function-type"
#endif
  


#endif  // !standalone
