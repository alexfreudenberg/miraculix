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


#include "AutoRandomFieldsUtils.h"

const char

*R_TYPE_NAMES[LAST_R_TYPE_NAME + 1] = { // never change ! see kleinkram.cc
  "NILSXP" /* 0 */,
  "SYMSXP", "LISTSXP", "CLOSXP", "ENVSXP", "PROMSXP",
  "LANGSXP", "SPECIALSXP", "BUILTINSXP", "CHARSXP", "LGLSXP" /* 10 */,
  "??", "??", "INTSXP", "REALSXP", "CPLXSXP",
  "STRSXP", "DOTSXP", "ANYSXP", "ECSXP", "EXPRSXP" /*20 */,
  "BCODESXP", "EXTPTRSXP", "WEAKREFSXP", "RAWSXP", "S4SXP" /* 25 */,
  "", "", "", "", "NEWSXP" /* 30 */,
  "FREESXP", "??SXP"},

 *LA_NAMES[LA_LAST + 1] = { "intern", "R", "auto", "GPU", "query"},

*PIVOT_NAMES[PIVOT_LAST + 1] = {"no privoting",
  "do", "auto", "idx", "undefined"},

*INSTALL_NAMES[INSTALL_LAST + 1] = {"no installation", "install", "ask", 
				    "sse", "sse2", "sse3", "ssse3", "avx",
				    "avx2", "avx512f", "gpu"};
