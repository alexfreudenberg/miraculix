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



#ifndef Miraculixxport_H
#define Miraculixxport_H 1

#include "def.h"
#include "General_utils.h"
#include "zzz_RFU.h"



#define UTILS_ALWAYS \
  CALL(startRFU);				\
  CALL(solve_NULL);				\
  CALL(solve_DELETE);				\
  CALL(orderingL);				\
  CALL(orderingInt);				\
  CALL(SolvePosDefSp);				\
  CALL(sleepMicro);				\
    CALL(get_utilsoption);			\
  CALL(get_utils_basic);			\
  CALL(push_utilsoption);			\
   CALL(params_utilsoption);			\
  CALL(del_utilsoption);			\
  CALL(scalarX);				\
  CALL(scalarXint);				\
  CALL(parallel);			      	\
  CALL(pid)

#if defined compatibility_to_R_h 
#define UTILSCALLS				\
  CALL(setoptionsRFU);		      		\
  CALL(getoptionsRFU);		      		\
  CALL(attachRFUoptions);			\
  CALL(detachRFUoptions);			\
  CALL(RFUoptions);		     		\
  CALL(attachSetNGet);		      		\
  UTILS_ALWAYS
extern SEXP Information, Filecoding, Filename, Next, Missings, Precise, Doubled;
#else
#define UTILSCALLS  UTILS_ALWAYS
#endif
  
#if defined CALL
#undef CALL
#endif
#define CALL(what) extern what##_type Ext_##what
UTILSCALLS;

extern const char *R_TYPE_NAMES[LAST_R_TYPE_NAME + 1];

void includeXport();
void PIDKEY_M_DELETE();
typedef
struct option_type option_type;
void WhichOptionList(bool local, option_type **global,
			     utilsoption_type **utils);
void load_utilsoptions(utilsoption_type *S, int local);

void startMiraculix(int n);

void FREEglobal();

Uint check_intrinsics();

#define check_cuda check_7_5 // the lowest number for which code exists


#endif
