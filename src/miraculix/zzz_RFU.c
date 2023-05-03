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

#include "def_rfu.h"
#include "Basic_RandomFieldsUtils.h"

#if defined compatibility_to_R_h

#include "win_linux_aux.h"
#include "utils.h"
#include "zzz_RFU.h"
#include "RandomFieldsUtils.h"

#define none 0

#if defined(__clang__)
//# pragma clang diagnostic ignored "-Wcast-function-type"
#endif

#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
//#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
 

static R_NativePrimitiveArgType 
    int_arg[] = { INTSXP },
    host_arg[] = { STRSXP, INTSXP};
  //  static R_NativeArgStyle argin[] = {R_ARG_IN},
  //    argout[] = {R_ARG_OUT},
  //   hostarg[] = {R_ARG_OUT, R_ARG_OUT};

#define CDEF(name, n, type) {#name, (DL_FUNC) & name, n, type}
static const R_CMethodDef cMethods[]  = {
  cRFUMethods
  {NULL, NULL, 0, NULL}
};


#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  //  CALLDEF(),
  callRFUMethods
  {NULL, NULL, 0}
};


 
#define EXTDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_ExternalMethodDef extMethods[] = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  extRFUMethods
  {NULL, NULL, 0} 
};

#define CALLABLE(FCTN)  R_RegisterCCallable("RandomFieldsUtils", #FCTN, (DL_FUNC)  FCTN)

void R_init_RandomFieldsUtils(DllInfo  *dll) {
  RFU_CALLABLE;

  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     extMethods);
  R_useDynamicSymbols(dll, FALSE); //
}


#ifdef SCHLATHERS_MACHINE
#ifdef __GNUC__ 
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
#endif
void R_unload_RandomFieldsUtils(DllInfo *info) { }
#ifdef __GNUC__ 
//#pragma GCC diagnostic pop
#endif


#endif  // !standalone
