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


#ifndef miraculix_vector_matrix_template_H
#define miraculix_vector_matrix_template_H 1

#include "options.h"


typedef double floatD;
typedef double longD;

typedef void (*coding_t)(unit_t *, Long, Long, Long, Long, Long,
			      int, double *, unit_t*, Long, Long);

#define coding_header(UINT,NAME)				\
  void coding##NAME##UINT(unit_t *M,  Long ldM,				\
			     Long start_snp, Long end_snp,		\
			  Long start_individual, Long end_individual,	\
			     int VARIABLE_IS_NOT_USED cores,		\
			     double VARIABLE_IS_NOT_USED *G,		\
			     unit_t* Ans, \
			     Long VARIABLE_IS_NOT_USED individuals,	\
			     Long ldAns)			

#define get_matrix_header(UINT,NAME)					\
  void get_matrix##NAME##UINT(unit_t *code, Long VARIABLE_IS_NOT_USED snps, \
				Long individuals, Long lda,		\
			      int VARIABLE_IS_NOT_USED cores,		\
			      unit_t *Ans, Long ldAns) 	     



#define vector_header(TYPE, ENDTYPE)					\
  /* !!! not snps / indivduals for FiveCodesTransposed !!! */		\
  (unit_t *code, Long rows, Long cols, Long lda,			\
   coding_type VARIABLE_IS_NOT_USED coding,				\
   int VARIABLE_IS_NOT_USED variant,					\
   LongDouble VARIABLE_IS_NOT_USED * pseudoFreq,			\
   TYPE *V, Long repetV, Long ldV,					\
   basic_options VARIABLE_IS_NOT_USED *opt,				\
   tuning_options VARIABLE_IS_NOT_USED *tuning,				\
   ENDTYPE *Ans, Long ldAns)


typedef void (*vector_double_t) vector_header(double, double);
typedef void (*vector_LongDouble_t) vector_header(LongDouble, LongDouble);
typedef void (*vector_Ulong_t) vector_header(Ulong, Ulong);
typedef void (*vector_floatD_t) vector_header(double, double);

#define vectorGeno_header(TYPE, TRDTYPE, ENDTYPE, NAME) \
  void vectorGeno##NAME##_##TYPE vector_header(TYPE, ENDTYPE)



#define genoVector_header(TYPE, TRDTYPE, ENDTYPE, NAME)		\
  void genoVector##NAME##_##TYPE vector_header(TYPE, ENDTYPE)



#define sumGeno_header0					\
  (unit_t *code, Long snps, Long individuals, Long lda,	\
   int VARIABLE_IS_NOT_USED cores, Ulong *sums)
#define sumGeno_header(NAME) Ulong sumGeno##NAME sumGeno_header0	\

typedef Ulong (*sumGeno_t) sumGeno_header0;


#define sumGeno_start(NAME)			\
  sumGeno_header(NAME) {			\
  Ulong VARIABLE_IS_NOT_USED total = 0UL;


#define gV_vG_header0(TYPE)						\
  (SEXP SxI, TYPE *V, Long repetV, Long ldV,				\
   bool orig_gV, bool SubstractMeanSxI, 				\
   option_type *global, utilsoption_type *utils,			\
   TYPE *ans, Long ldAns)

#define gV_vG_header(TYPE, NAME)					\
  void gV_vG##NAME##_##TYPE gV_vG_header0(TYPE) 

typedef void (*gV_vG_double_t) gV_vG_header0(double);
typedef void (*gV_vG_LongDouble_t) gV_vG_header0(LongDouble);
typedef void (*gV_vG_Ulong_t) gV_vG_header0(Ulong);


#define sparseTGeno_Header /* (sparse) %*% t(code) */		\
  (Uchar *code,								\
   Long nrow_code, Long ncol_code, Long ldaInByte,			\
   coding_type VARIABLE_IS_NOT_USED coding,				\
   double *valueB,							\
   int nIdx,/*2nd dimension of sparse; == len(rowIdxB)*/		\
   int *rowIdxB, int *colIdxB,						\
   bool VARIABLE_IS_NOT_USED tSparse,					\
   option_type VARIABLE_IS_NOT_USED *global,				\
   utilsoption_type VARIABLE_IS_NOT_USED *utils,			\
   double *Ans, Long ldAns)

#define sparseTGeno_header(NAME) /* (sparse) %*% t(code) */		\
  void sparseTGeno##NAME sparseTGeno_Header


#define CodesPerLong ((int) (sizeof(Long) * BitsPerByte / BitsPerCode))
#define sparseTGeno_start(NAME)						\
  sparseTGeno_header(NAME){						\
  if (nrow_code > (Long)MAXINT || nrow_code < (Long)CodesPerLong ||	\
      ncol_code > (Long)MAXINT) BUG;					\
  const basic_options *opt = &(utils->basic);				\
  const int cores = GreaterZero(opt->cores);
 
  
typedef void (*sparseTGeno_t) sparseTGeno_Header;
		 

#endif
