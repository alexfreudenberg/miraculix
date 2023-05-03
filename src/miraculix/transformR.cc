
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

#if defined compatibility_to_R_h
#include "miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "utils_miraculix.h"
#include "MXinfo.h"
#include "transform.h"
#include "Files.h"


SEXP Transform(SEXP SxI, SEXP Coding, SEXP SelSNPS, SEXP SelIndiv) {
  // from user: Coding == oldcoding
  // to user: Coding == newcoding
  //
  // (i) AutoCoding or given SelSNPS/SelIndiv makes always a copy
  //     AutoCoding takes the coding from options$snpcoding
  // (ii) UnknownSNPcoding or (identical old coding and new coding)
  //            completes information without copying
  
   // if (TYPEOF(SxI) == INTSXP)  printf("hier %d %d\n", INTEGER(SxI)[0], INTEGER(SxI)[1]); else if (TYPEOF(SxI) == STRSXP) printf("file %d %s\n", LENGTH(SxI), CHAR(STRING_ELT(SxI, 0))); else printf("Aaahh %f %f\n", REAL(SxI)[0], REAL(SxI)[1]);
  Long
    lenSnps = LENGTH(SelSNPS),
    lenIndiv = LENGTH(SelIndiv);
  int
    *selSnps = lenSnps == 0 ? NULL : INTEGER(SelSNPS),
    *selIndiv = lenIndiv == 0 ? NULL : INTEGER(SelIndiv);
  const coding_type origCoding = (coding_type) INTEGER(Coding)[0];
  
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  Long *info = GetInfoUnchecked(SxI);
  unit_t *code = NULL;
  bool doFree = false;
  if (info == NULL) {
    doFree = true;
    ToIntGlobal(SxI);
    code = (unit_t*) SxIint;
  }
  SEXP Ans;
  PROTECT(Ans = Transform(SxI, code, origCoding, selSnps, lenSnps,
			  selIndiv, lenIndiv,
			  global, utils));  
  if (doFree) FREEglobal();
  UNPROTECT(1);
  return Ans;
}


SEXP transpose(SEXP SxI){
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  return transpose(SxI, global, utils);
 }

#endif
