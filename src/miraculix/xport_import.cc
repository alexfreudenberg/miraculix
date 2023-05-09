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

#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "MXinfo.h"
//#include "2bit.h"
//#include "3bit.h"
//#include "5codes.h"
#include "mmagpu.h"
#include "intrinsics_specific.h"


#define importfrom "RandomFieldsUtils"

#if defined CALL
#undef CALL
#endif
#define CALL(what) what##_type Ext_##what = NULL
UTILSCALLS;


#if defined CALL
#undef CALL
#endif

#if defined compatibility_to_C_h
#define CALL(what) Ext_##what = (what##_type) what; 
#else 
#define CALL(what) Ext_##what = (what##_type) R_GetCCallable(importfrom, #what)
#endif
void includeXport() {
  UTILSCALLS;
} 
  


coding_type getAutoCodingIntern() { // used also in options.h
  return TwoBitGeno;
}


void load_utilsoptions(utilsoption_type *S, int local) {
  int params[N_UTILS_PARAM];
  Ext_params_utilsoption(local, params);
  S->solve.n_pivot_idx = params[PIVOT_IDX_N];
  S->solve.pivot_idx = S->solve.n_pivot_idx == 0 ? NULL
    : (int*) MALLOC(S->solve.n_pivot_idx * sizeof(*(S->solve.pivot_idx)));
  Ext_get_utilsoption(S, local);
}




THIS_FILE_ANYSIMD;
EXTERN_SIMD_CHECK(2bit128);
EXTERN_SIMD_CHECK(2bit256);
EXTERN_SIMD_CHECK(2bit512);
EXTERN_SIMD_CHECK(mma_61);
EXTERN_SIMD_CHECK(mma_75);
Uint check_intrinsics() { 
  //  if (BytesPerUnit != sizeof(Uint)) // OK
  //    ERR1("The programme relies on the assumption that an unsigned integer has %u Bytes.", (Uint) BytesPerUnit);// OK
  
  if (sizeof(int) != sizeof(Uint))
    ERR2("The programme relies on the assumption that a signed integer the same lenth than an unsigned integer. Found %u and %u bytes respectively.",
	 (Uint) sizeof(int),  (Uint) sizeof(Uint));
  
  if (SizeOfDouble != sizeof(double))
    ERR2("The programme relies on the assumption that size of a 'double' is %d. Found %u bytes.", SizeOfDouble, (Uint) sizeof(double));
  
   if (SizeOfInt != sizeof(Uint)) 
    ERR2("The programme relies on the assumption that size of a 'int' is %d. Found %u bytes.", SizeOfInt, (Uint) sizeof(Uint));

   if (sizeof(unit_t) != BytesPerUnit) 
     ERR2("The programme relies on the assumption that size of 'unit_t' is %d. Found %u bytes.", BytesPerUnit, (Uint) sizeof(unit_t));
 
  

  Uint maxsize = MaxUnitsPerAddress * BytesPerUnit; // OK
  if (sizeof(void*) > maxsize)
    ERR2("The programme relies on the assumption that an 'void*' has at most %u Bytes. Found %u bytes.", 
	 maxsize, (Uint) sizeof(void*));

  if (!sseAvail)
    RFERROR("programm does not run on machine that old (even not having sse)\n");

  CHECK_THIS_FILE_ANYSIMD;
  CHECK_FILE(2bit128); // sufficient to check 2bit* only
  CHECK_FILE(2bit256);
  CHECK_FILE(2bit512);
  CHECK_FILE(mma_75);
  CHECK_FILE(mma_61);
  return SIMD_INFO;
}


void InitPlink();
void Init2();
void Init3();
void Init5();


void startRFU();
void startMiraculix(int n) {
   check_intrinsics();
#if defined compatibility_to_C_h
  startRFU();
#endif  
  includeXport();
  OPTIONS_MIRACULIX.genetics.coding = getAutoCodingIntern();
  const bool local=false;
  utilsoption_type S;
  load_utilsoptions(&S, local);
  S.solve.max_chol = 50000;
  S.solve.max_svd = 6555;
  S.solve.Methods[0]=S.solve.Methods[1]=Cholesky;
  S.solve.Methods[2]=NoFurtherInversionMethod;
  S.basic.cores = GreaterZero(n);

  //   printf("start bigendian = %d\n", S.basic.bigendian);
 
  Ext_push_utilsoption(&S, local);
  Ext_del_utilsoption(&S);

  Init2();
  Init3();
  Init5();
  InitPlink();
  check_cuda(); // uses  KEY_M
}



#if defined compatibility_to_C_h  
extern utilsoption_type OPTIONS;
#endif
void WhichOptionList(bool VARIABLE_IS_NOT_USED local, option_type **global,
		     utilsoption_type **utils // readonly
		     ) {
  //WhichOptionsList considers 'local' whenever possible, otherwise alternatives
  option_type *G = &OPTIONS_MIRACULIX;
  utilsoption_type *U = NULL;
#if defined compatibility_to_R_h
  KEY_type *KT = KEYT_M();
  if (KT == NULL) BUG;
  U = &(KT->global_utils);  // readonly!
  if (local) G = &(KT->global);
#else
  U = &OPTIONS;
#endif
  if (global != NULL) *global = G;
  if (utils != NULL) *utils = U;
}


