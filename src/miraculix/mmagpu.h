/*
 Authors 
 Alexander Freudenberg, alexander.freudenberg@uni-mannheim.de
 
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



#ifndef miraculix_mmagpu_H
#define miraculix_mmagpu_H 1

void crossprod_mmagpu(Uint* code, Long snps, Long individuals,
		      Uint Lda, int cores, double *ans);
void check_7_5();

void crossprod_PlainGPU(Uint* M, Long snps, Long individuals,
			Uint VARIABLE_IS_NOT_USED lda, int cores, double* A );

void err_check(const char* string);

#endif
