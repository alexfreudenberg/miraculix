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

#ifndef miraculix_auto_H
#define miraculix_auto_H 1

Long LdaPlain(Long snps, Long ldAbitalign);
Long LdaPlainDoubled(Long snps, Long ldAbitalign);
sumGeno_header(Plain);
coding_header(_4Byte, Plain);
coding_header(_1Byte, Plain);
get_matrix_header(_4Byte,Plain);
get_matrix_header(_1Byte,Plain);

void TwoBithaplo2genoPlain(unit_t * Code, Long snps, Long individuals, Long lda,
			   int cores, unit_t *A);
void crossprod_Plain(unit_t * Code, Long snps, Long individuals, Long lda,
		     int cores, double *A);

#endif

