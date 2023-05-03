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


#ifndef miraculix_2bit_H
#define miraculix_2bit_H 1

//#include "MX.h"



#define DeclareVersion2(NR,SU,SV, SM)					\
  Ulong scalar_2v##NR(unit_t V *x, unit_t V *y,Long V blocks,Long V DeltaPart2) SU \
    void multi_2v##NR(unit_t V *x, unit_t V  *y, Long V individuals, \
		      Long V blocks, Long V lda,		 \
		      double V  *ans) SV		\
    void multi2x2_2v##NR(unit_t V *x, unit_t V *y, Long V individuals, \
			 Long V blocks, Long V lda,			\
			 double V *ans) SV				\
    Ulong sumGeno2v##NR(unit_t V *S, Long V snps, Long V individuals,	\
			Long V lda, int V cores, Ulong V *sums) SU	\
    /* subsequently function working only under bigendian */		\
    void TwoBithaplo2geno##NR(unit_t V * Code, Long V snps, Long V individuals, \
			      Long V ldH, int V cores,			\
			      unit_t V *A, Long V ldAns) SV		\
    void allele_sum_2v##NR(unit_t V *Code, Long V snps, Long V individuals, \
			   Long V lda, Long V LDAbitalign,		\
			   basic_options V *opt,			\
			   Ulong V *ans) SV				\
  void coding_2v##NR(unit_t V *M, Long V ldM,				\
		       Long V  start_snp, Long V end_snp,		\
		       Long V start_individual, Long V end_individual,	\
		       basic_options V *opt, double V *G, SEXP  V Ans) SV \
  void coding_2v##NR##_4Byte(unit_t V *M,  Long V ldM,			\
			      Long V start_snp, Long V end_snp,		\
			      Long V start_individual, Long V end_individual, \
			      int V cores, double V *G,			\
			      unit_t V *ans, Long V individuals,Long V lda) SV \
  void trafo2Geno1Geno##NR(unit_t V  *OLD, Long V snps, Long V individuals, \
			     Long V oldLDA, unit_t V *ANS, Long V ldAns) SV \
      void TrafoFileValues_##NR(unit_t V *CodeFile, Long V blocks) SV


#define V 
#define NICHTS ;
DeclareVersion2(512,NICHTS,NICHTS,NICHTS)
DeclareVersion2(256,NICHTS,NICHTS,NICHTS)
DeclareVersion2(128,NICHTS,NICHTS,NICHTS)
DeclareVersion2(64,NICHTS,NICHTS,NICHTS)
#if defined V 
#undef V
#undef NICHTS
#endif

void Init2();

Long Lda2Bit(Long snps, Long ldAbitalign);
sumGeno_header(2);


coding_header(_4Byte, 2);
coding_header(_1Byte, 2);
coding_header(_4Byte, 2trans);
coding_header(_1Byte, 2trans);

bool coding_2v128(unit_t *pM, Long blocks, Long rest, uint32_t* ans);

get_matrix_header(_4Byte,2);
get_matrix_header(_1Byte,2);


genoVector_header(double, double, double, 2matrix256);
genoVector_header(double, double, double, 2v256);
vectorGeno_header(double, double, double, 2v256);
gV_vG_header(LongDouble, 2);
gV_vG_header(double, 2);
gV_vG_header(Ulong, 2);

void allele_sum2(unit_t *Code, Long snps, Long individuals,
		 Long lda, 
		 basic_options *opt, Ulong *ans);


void zeroGeno2(SEXP SxI, Uint *Snps, Uint lenSnps, Uint *Indiv, Uint lenIndiv,
	       basic_options *opt);



void crossprod2v64(unit_t * M, Long snps, Long individuals, Long lda,
		   int cores, double *A);
void crossprod2v32(unit_t * M, Long snps, Long individuals, Long lda,
		   int cores, double *A);


void TwoBithaplo2geno2(unit_t * Code, Long snps, Long individuals,
		       Long lda, int cores, unit_t *A, Long ldAns);

sparseTGeno_header(2Bit);
		   
#endif
