
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


#ifndef compatibility_general_h 
#define compatibility_general_h 1

#ifndef __cplusplus
#include <stdbool.h>
#endif
#include "compatibility.header.h"
#include <inttypes.h>

typedef unsigned int Uint;
typedef uint64_t Ulong;
typedef int64_t Long;
typedef unsigned char Uchar;

#define F77call F77_CALL // rename to control that USE_FC_LEN_T has been called
#ifdef __cplusplus
#define F77name extern "C" void F77_NAME // rename to control that USE_FC_LEN_T has been called
#else
//#warning not cplusplus
#define F77name void F77_NAME 
#endif
#define F77dgesdd F77call(dgesdd)
#define F77dgemv F77call(dgemv)
#define F77ddot F77call(ddot)
#define F77dsyrk F77call(dsyrk)


void stopIfNotIntI(Long i, Long line, const char *file);
#define stopIfNotInt(i) stopIfNotIntI(i, __LINE__, __FILE__);

void stopIfNotUIntI(Long i, Long line, const char *file);
#define stopIfNotUInt(i) stopIfNotUIntI(i, __LINE__, __FILE__);

void stopIfNotAnyIntI(Long i, Long line, const char *file);
#define stopIfNotAnyInt(i) stopIfNotAnyIntI(i, __LINE__, __FILE__);

#define stopIfNotSame(i,j)\
  if (sizeof(i) != sizeof(j)) { ERR6("'%s' (%ld) and '%s' (%ld) do not have the same size at line %d in '%s'\n", #i, sizeof(i), #j, sizeof(j), __LINE__, __FILE__) }


#if defined STAND_ALONE
//#warning STAND_ALONE
#include "compatibility.C.h"
#else
//#warning R_VERSION
#include "compatibility.R.h"
#endif


#if defined compatibility_to_R_h
typedef double LongDouble;
#else
typedef long double LongDouble;
#endif

#if !defined STAND_ALONE || defined DO_PARALLEL
#define ASSERT_SOLVE(sp) assert(sp != NULL);
#else
#define ASSERT_SOLVE(sp) if (sp == NULL) sp = &(OPTIONS.solve); else {}
#endif


#endif
