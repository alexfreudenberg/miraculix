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

#ifndef miraculix_avx_miss_H
#define miraculix_avx_miss_H 1

// !defined appropriate AVX
void static AVXmissing() { ERR3("ERROR: '%sv%d' needs appropriate AVX features (%s).\n", CODING_NAMES[MY_CODING], MY_VARIANT, __FILE__ ); }
#define Sm { AVXmissing(); return R_NilValue; }
#define Su { AVXmissing(); return 0; }
#define Sv { AVXmissing(); }
#if defined V 
#undef V
#endif
#if defined VARIABLE_IS_NOT_USED
#define V VARIABLE_IS_NOT_USED
#else
#define V
#endif



#endif
