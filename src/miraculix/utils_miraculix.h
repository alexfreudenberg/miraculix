
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


#ifndef miraculix_utils_H
#define miraculix_utils_H 1

int *ToIntI(SEXP X, bool *create, bool round);

#define ToIntDo(X, CREATE)					\
  bool create##X = CREATE;					\
  Uint *X##int = (Uint*) ToIntI(X, &create##X, false)

#define ToIntLocal(X) ToIntDo(X, true)	
#define ToIntGlobal(X) ToIntDo(X, false)	

#define FREElocal(X) if (!create##X) {} else FREE(X##int)

#endif
