/*
 Authors 
 Alexander Freudenberg, afreuden@mail.uni-mannheim.de

 Copyright (C) 2022-2023 Martin Schlather, Alexander Freudenberg

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


#include "def_rfu.h"
#include "Basic_RandomFieldsUtils.h"
#if defined compatibility_to_R_h
#include "RandomFieldsUtils.h"
#include "errors_messages.h"

#if defined USEGPU
SEXP gpu_info_61(SEXP DEVICES);
#endif


SEXP gpu_info(SEXP DEVICES){
#if defined USEGPU
  return gpu_info_61(DEVICES);
#else  
  ERR0("No CUDA devices found");
  return R_NilValue;
#endif
}

#endif
