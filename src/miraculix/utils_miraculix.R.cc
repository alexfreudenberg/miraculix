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



#define NO_OMP 1
#include "Basic_miraculix.h"
#if defined compatibility_to_R_h

#include "xport_import.h"
#include "utils_miraculix.h"
#include "MXinfo.h"
#include "kleinkram.h"
#include "extern.h"
#include "options.h"


#define UTILSINFO_M(M) if (!KEYT_M()->global_utils.basic.helpinfo) {} else PRINTF("%s\n(Note that you can unable this information by 'RFoptions(helpinfo=FALSE)'.)\n", M) // OK

										    
int *ToIntI(SEXP X, bool *create, bool round) {
  int len = LENGTH(X);
  if (len == 0) {*create=false; return NULL;}
  if (TYPEOF(X) == INTSXP) {
    *create = false;
    return INTEGER(X);
  }
  if (TYPEOF(X) == LGLSXP) {
    *create = false;
    return LOGICAL(X);
  }
  int *y;
  KEY_type *KT = KEYT_M();
  utilsoption_type *global_utils = &(KT->global_utils);
  
  if (len > 100 || global_utils->basic.Rprintlevel > 1) {
    UTILSINFO_M("Better use 'integer' as storage mode (for one of the arguments)."); }

  if (*create || KT->ToIntN < len) {
    y = (int *) MALLOC(sizeof(*y) * len);    
    if (y == NULL) ERR1("not enough memory for an %d vector of integers", len);
    if (!*create) {
      FREE(KT->ToIntDummy);
      KT->ToIntDummy = y;
      KT->ToIntN = len;
    }
  } else {
    y = KT->ToIntDummy;
  }
  double *x = (double *) REAL(X);
  if (round) for (int i=0; i<len; i++) y[i] = (int) ROUND(x[i]);
  else for (int i=0; i<len; i++) y[i] = (int) x[i];
  return y;
}



#endif
