/*
 Authors 
 Alexander Freudenberg, alexander.freudenberg@stads.de

 Copyright (C) 2022-2023 Alexander Freudenberg

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



/*
    This file replaces the R-specific types and functions in mmagpu.cu to allow the utilities to be used outside of miraculix's R version.

*/
#ifndef AUX_MMAGPU 
#define AUX_MMAGPU 1

#define CodesPerByte 4L
#define BytesPerUnit 4L
#define Uint unsigned int
#define Ulong unsigned long

#define SEXP int
void PROTECT(SEXP a){}
void UNPROTECT(SEXP a){}

#define PRINTF printf
#define MEMCOPY memcpy
void ERR(const char* s){
    printf("%s",s);
    exit(1);
}
long static inline UnitsPerIndiv256(Ulong snps){
    return (1L + (snps - 1L)/((4L * 8L / 2) * (32L / 4L))) * (32L / 4L);
}  

#endif