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




#ifndef miraculix_OneByte_H
#define miraculix_OneByte_H 1


Long LdaOneByte(Long snps, Long ldAbitalign);
sumGeno_header(OneByte);

void OneByte2T(unit_t *tmp, Long snps, Long individuals, Long lda,
	       unit_t *Ans, Long ldAns);
void OneByteHaplo2geno64(unit_t *H, Long snps, Long individuals, Long lda,
		      int cores, unit_t *Ans, Long ldAns);

coding_header(_4Byte, OneByte);
coding_header(_1Byte, OneByte);

get_matrix_header(_4Byte,OneByte);
get_matrix_header(_1Byte,OneByte);

gV_vG_header( LongDouble, OneByte);
gV_vG_header( double, OneByte);
gV_vG_header( Ulong, OneByte);



void allele_sum_OneByte(unit_t *Code, Long snps, Long individuals, 
			Long lda, Long LDAbitalign, basic_options *opt,
			Ulong *ans);


sparseTGeno_header(OneByte);


#endif
