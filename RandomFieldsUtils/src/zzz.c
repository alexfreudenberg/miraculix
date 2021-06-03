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

//#include "Basic_utils.h" // must be before anything else

#include "RandomFieldsUtils.h"
#include "zzz_RandomFieldsUtils.h"
#include "Utils.h"

#define none 0

#if defined(__clang__)
//# pragma clang diagnostic ignored "-Wcast-function-type"
#endif

#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
// GCC diagnostic ignored "-Wcast-function-type"
#endif

static R_NativePrimitiveArgType 
    int_arg[] = { INTSXP },
    host_arg[] = { STRSXP, INTSXP};
  //  static R_NativeArgStyle argin[] = {R_ARG_IN},
  //    argout[] = {R_ARG_OUT},
  //   hostarg[] = {R_ARG_OUT, R_ARG_OUT};

#define CDEF(name, n, type) {#name, (DL_FUNC) & name, n, type}
static const R_CMethodDef cMethods[]  = {
  CDEF(sleepMilli,  1, int_arg),
  CDEF(sleepMicro, 1, int_arg),
  CDEF(pid, 1, int_arg),
  CDEF(hostname, 2, host_arg),
  CDEF(setCPUs, 1, int_arg),
  CDEF(recompilationNeeded, 1, int_arg),
  CDEF(loadoptions, 0, none),
  CDEF(detachoptions, 0, none),
  {NULL, NULL, 0, NULL}
};


#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF(AVXmessages, 1),
  CALLDEF(DebugCall, 0),
  CALLDEF(Chol, 1),
  CALLDEF(SolvePosDef, 3),
  CALLDEF(struve, 4),
  CALLDEF(besselk_simd, 2),
  CALLDEF(I0ML0, 1),
  CALLDEF(gaussr, 2),
  CALLDEF(WMr, 4),
  CALLDEF(logWMr, 4),
  CALLDEF(sortX, 4),
  CALLDEF(orderX, 4), 
  CALLDEF(getChar, 0),
  CALLDEF(DivByRow, 2),
  CALLDEF(colMaxs, 1),
  CALLDEF(quadratic, 2),
  CALLDEF(dotXV, 2),
  CALLDEF(rowMeansX, 2),
  CALLDEF(rowProd, 1),
  CALLDEF(dbinorm, 2),
  CALLDEF(chol2mv, 2),
  CALLDEF(tcholRHS, 2),
  CALLDEF(crossprodX, 3),
  CALLDEF(getPackagesToBeInstalled, 1),
  //  CALLDEF(),
  {NULL, NULL, 0}
};


 
#define EXTDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_ExternalMethodDef extMethods[] = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  EXTDEF(RFoptions, -1), 
  {NULL, NULL, 0} 
};



#define CALLABLE(FCTN)  R_RegisterCCallable("RandomFieldsUtils", #FCTN, (DL_FUNC)  FCTN)
void R_init_RandomFieldsUtils(DllInfo  *dll) {
  CALLABLE(utilsoption_DELETE);
  CALLABLE(utilsoption_NULL);

  CALLABLE(solve_DELETE);
  CALLABLE(solve_NULL);
  CALLABLE(solvePosDef);
  CALLABLE(invertMatrix);
  
  CALLABLE(solvePosDefSp);
  CALLABLE(sqrtPosDefFree); 
  CALLABLE(sqrtRHS);
  
  CALLABLE(detPosDef);
  CALLABLE(detPosDefsp);
  CALLABLE(XCinvXdet);
  CALLABLE(XCinvYdet);
  CALLABLE(is_positive_definite);
  CALLABLE(chol2inv);
  CALLABLE(chol);

  CALLABLE(StruveH);
  CALLABLE(StruveL);
  CALLABLE(I0mL0);

  CALLABLE(WM);
  CALLABLE(DWM);
  CALLABLE(DDWM);
  CALLABLE(D3WM);
  CALLABLE(D4WM);
  CALLABLE(logWM);
  
  CALLABLE(Gauss);
  CALLABLE(DGauss);
  CALLABLE(DDGauss);
  CALLABLE(D3Gauss);
  CALLABLE(D4Gauss);
  CALLABLE(logGauss);
  
  CALLABLE(getUtilsParam);
  CALLABLE(attachRFoptions);
  CALLABLE(detachRFoptions);

  CALLABLE(ordering);
  CALLABLE(orderingInt);
  CALLABLE(sorting);
  CALLABLE(sortingInt);
  CALLABLE(scalarX);
  //  CALLABLE(scalarInt);

  CALLABLE(pid);
  CALLABLE(sleepMicro); // problem?
 
  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     extMethods);
  R_useDynamicSymbols(dll, FALSE); //
}


#ifdef SCHLATHERS_MACHINE
#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
// GCC diagnostic push
// GCC diagnostic ignored "-Wcast-function-type"
#endif
#endif
void R_unload_RandomFieldsUtils(DllInfo *info) { }
#ifdef __GNUC__
// GCC diagnostic pop
#endif

