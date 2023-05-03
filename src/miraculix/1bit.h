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

#ifndef miraculix_1bit_H
#define miraculix_1bit_H 1


#define DeclareVersion1(NR,SU,SV)					\
  Ulong scalar_1v##NR(unit_t V *x, unit_t V *y,Long V blocks,Long V DeltaPart2) SU \
  void multi_1v##NR(unit_t V *x, unit_t V  *y, Long V individuals, \
		    Long V blocks, Long V lda,		       \
		    double V  *ans) SV					\
    Ulong scalar_1v##NR##A(unit_t V *x, unit_t V *y,Long V blocks,	\
			Long V DeltaPart2) SU			  \
  void multi_1v##NR##A(unit_t V *x, unit_t V  *y, Long V individuals, \
		       Long V blocks, Long V lda,		  \
		      double V  *ans) SV		\
  Ulong scalar_1v##NR##B(unit_t V *x, unit_t V *y,Long V blocks,Long V DeltaPart2) SU \
  void multi_1v##NR##B(unit_t V *x, unit_t V  *y, Long V individuals, \
		       Long V blocks, Long V lda,		  \
		       double V  *ans) SV		\
    void OneBithaplo2geno##NR(unit_t V *X, Long V snps, Long V individuals, \
			      Long V lda,				\
			      int V cores, unit_t V *Ans, Long V ldAns) SV \
  void trafo1Geno2Geno##NR(unit_t V *OLD, Long V snps, Long V indiv, \
			   Long V oldLDA,			  \
			   unit_t V *Ans, Long V ldAns) SV		      
 

#if defined V 
#undef V
#endif
#define V
#define NICHTS ;
DeclareVersion1(512,NICHTS,NICHTS)
DeclareVersion1(256,NICHTS,NICHTS)
DeclareVersion1(128,NICHTS,NICHTS)
#if defined V 
#undef V
#undef NICHTS
#endif

Long Lda1Bit(Long snps, Long ldAbitalign);
sumGeno_header(1);

coding_header(_4Byte,1);
coding_header(_1Byte,1);

get_matrix_header(_4Byte,1);
get_matrix_header(_1Byte,1);


Ulong scalar_1v64(unit_t *x, unit_t *y, Long blocks, Long Delta);
Ulong scalar_1v32(unit_t *x, unit_t *y, Long blocks, Long Delta);
void OneBithaplo2geno64(unit_t *X, Long snps, Long individuals, Long lda,
			int cores, unit_t *Ans, Long ldAns);
void OneBithaplo2geno32(unit_t *X, Long snps, Long individuals, Long lda,
			int cores, unit_t *Ans, Long ldAns);
void multi_1v64(unit_t *x, unit_t *y, Long individuals, Long blocks, Long lda,
		double *ans);
void multi_1v32(unit_t *x, unit_t *y, Long individuals, Long blocks, Long lda,
		 double *ans);

void zeroGeno1(SEXP SxI, Uint *Snps, Uint lenSnps, Uint *Indiv, Uint lenIndiv,
	       basic_options *opt);


#endif
