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




#ifndef miraculix_plink_H
#define miraculix_plink_H 1


Long LdaPlink(Long snps, Long LDAbitalign);
Long LdaOrigPlink(Long snps, Long LDAbitalign);

coding_header(_4Byte, Plink);
coding_header(_1Byte, Plink);
get_matrix_header(_4Byte,Plink);
get_matrix_header(_1Byte,Plink);
vectorGeno_header(Ulong, Ulong, Ulong, Plink);
vectorGeno_header(double, double, double, Plink);
vectorGeno_header(LongDouble, LongDouble, LongDouble, Plink);
sumGeno_header(Plink);

sparseTGeno_header(Plink);
vectorGeno_header(double, double, double, PlinkMatrix256);

void getPlinkMissings(Uchar* plink, SEXP SxI, basic_options *opt);

coding_header(_4Byte, Plinktrans);
coding_header(_1Byte, Plinktrans);

#if defined CODE_ONE
#undef CODE_LASTZERO
#undef CODE_ONE
#undef CODE_TWO
#endif
#define CODE_LASTZERO 1
#define CODE_ONE 2
#define CODE_TWO 3

#endif
