
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


#ifndef randomfieldsutils_extern_H
#define randomfieldsutils_extern_H 1


#include "AutoRandomFieldsUtilsLocal.h"
#include "zzz_RFU.h"

extern int PLoffset, PL;
extern utilsoption_type OPTIONS; 

#define prefixN 3
extern const char * prefixlist[prefixN], **allOptions[prefixN];
extern int allOptionsN[prefixN];


// AutoRFU
extern const char *LA_NAMES[LA_LAST + 1], *PIVOT_NAMES[PIVOT_LAST + 1],
  *INSTALL_NAMES[INSTALL_LAST + 1];

//extern const char *basic[basicN];
extern const char * InversionNames[nr_InversionMethods];
//extern const char * solve[solveN];
// extern bool ToFalse[1];
//

extern int parentpid;


#endif
