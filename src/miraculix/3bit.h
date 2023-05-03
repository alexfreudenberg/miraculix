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

#ifndef miraculix_3bit_H
#define miraculix_3bit_H 1


void Init3();
Long Lda3(Long snps, Long ldAbitalign);
sumGeno_header(3);

coding_header(_4Byte, 3);
coding_header(_1Byte, 3);
get_matrix_header(_4Byte,3);
get_matrix_header(_1Byte,3);

void haplo2geno3(unit_t * Code, Long snps, Long individuals, Long ldH,
		 int cores, unit_t *A);
void allele_sum3(unit_t *Code, Long snps, Long individuals, Long lda,
		 basic_options *opt, Ulong *ans);
void zeroGeno3(SEXP SxI, Uint *Snps, Uint lenSnps, Uint *Indiv, Uint lenIndiv,
	       basic_options *opt);

void TwoBithaplo2geno3(unit_t * Code, Long snps, Long individuals,
		       Long lda, int cores, unit_t *A, Long ldAns);
void crossprod3(unit_t * Code, Long snps, Long individuals, Long lda,
		int cores, double *A);



Long CodesPerBlock3();


#endif
