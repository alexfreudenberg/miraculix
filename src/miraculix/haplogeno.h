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

#ifndef miraculix_haplogeno_H
#define miraculix_haplogeno_H 1


#define TotalSumIndex 0
#define SqTotalSumIndex 1
#define SigmaSqIndex 2
#define PseudoSqTotalSumIndex 3
#define PseudoSigmaSqIndex 4
#define SumFreqIndex 5
#define PseudoSumFreqIndex 6
#define FreqIndex 7
#define SumIndex (FreqIndex + (Long) info[SNPS])
#define PseudoFreqSxIIndex (FreqIndex + 2L * (Long) info[SNPS])
#define FreqSxIIndex (FreqIndex + 3L *  (Long) info[SNPS])
#define PseudoFreqIndex (FreqSxIIndex + (Long) info[INDIVIDUALS])
#define PseudoSumIndex (FreqSxIIndex + 2L * (Long) info[INDIVIDUALS])

unit_t *Align(SEXP Code, basic_options *opt);
Long GetLDA(Long snps, Long individuals, coding_type coding, Long ldbitalign);
Long GetNatLDA(Long snps, Long individuals, coding_type coding);
Long SnpLDA(Long snps, coding_type coding, Long ldAbitalign);
Long SnpNatLDA(Long snps, coding_type coding);
Long GetLDAbitalign(coding_type coding);
Long totalMem(Long snps, Long individuals, Long lda,
	      coding_type coding, bool OneHaploOnly);
Long totalMem(Long snps, Long individuals, Long lda,
	      coding_type coding,  bool OneHaploOnly, bool relaxed);
SEXP createSNPmatrix(Long snps, Long individuals, coding_type coding,
                     int variant, Long ldAbitalign, bool OneHaploOnly, SEXP V,
                     option_type *global,  utilsoption_type *utils);  
SEXP createSNPmatrix(Long snps, Long individuals, coding_type coding,
                     int variant,  Long ldAbitalign,
                     option_type *global,  utilsoption_type *utils);
SEXP CompleteCodeVector(Long snps, Long individuals, coding_type coding,
			int variant, Long ldAbitalign, bool OneHaploOnly,
			option_type*global, utilsoption_type *utils,
			void* pointer, void* nextpointer, 
			SEXP Code);

SEXP CreateEmptyCodeVector(Long snps, Long individuals, coding_type coding,
			   int variant, Long ldAbitalign, bool OneHaploOnly,
			  option_type *global, utilsoption_type *utils);

Ulong sumGeno(unit_t * Code, Long snps, Long individuals, coding_type coding);

Long RoundUpSnps(Long snps, coding_type coding, int variant);

void allInfo(Long * info);
void allInfo(SEXP M);

Long calculateAlignedMem(Long memInU, coding_type coding);

Long GetBytesPerBlock(coding_type coding, int variant);
Long GetCodesPerBlock(coding_type coding, int variant);
Long GetCodesPerUnit(coding_type coding); // OK

void matrix_get_4Byte(unit_t *code, Long snps, Long indiv, Long lda, 
		     coding_type coding, int cores, unit_t *ans, Long ldAns);

void allele_freq_double(SEXP SxI, Long *mssngs, Long missings, bool pseudo,
			option_type *global, utilsoption_type *utils,
			double *sums, double *ans);

void allele_freq_LongDouble(SEXP SxI, Long *mssngs, Long missings, bool pseudo,
			    option_type *global, utilsoption_type *utils,
			    LongDouble *sums, LongDouble *ans);



LongDouble getTotalSum(SEXP SxI);
LongDouble *getFreq(SEXP SxI);
LongDouble *getSum(SEXP SxI);
LongDouble *getFreqSxI(SEXP SxI);
LongDouble getSqTotalSum(SEXP SxI);
LongDouble getSigmaSq(SEXP SxI);
LongDouble getSumFreq(SEXP SxI);
LongDouble *getPseudoFreq(SEXP SxI);
LongDouble *getPseudoSum(SEXP SxI);
LongDouble getPseudoSumFreq(SEXP SxI);
LongDouble *getPseudoFreqSxI(SEXP SxI);
LongDouble getPseudoSqTotalSum(SEXP SxI);
LongDouble getPseudoSigmaSq(SEXP SxI);

LongDouble *getFreq(SEXP SxI, option_type *global, utilsoption_type *utils);
LongDouble *getFreq(SEXP SxI, double *externalFreq,
		    option_type *global, utilsoption_type *utils);

void crossprod(unit_t * Code, Long snpsOrig, Long individuals, Long lda,
	       coding_type coding, int variant,
	       Long logSnpGroupSize, bool smoothSnpGroupSize,
	       Long indivGroupSize, basic_options *opt,
	       double *A);

double *crossprod(unit_t * Code, Long snps, Long individuals,  Long lda,
		  coding_type coding, int variant,		  
		  centering_type centered, normalizing_type normalized,
		  bool squared,
		  LongDouble *PP,
		  basic_options *opt, tuning_options *tuning);

coding_type switch_transposed(coding_type coding);

void sparseTGeno(SEXP SxI,
		 bool tSxI,
		double *valueB,
		int nIdx, // 2nd dimension of sparse; == length(rowIdxB)
		int *rowIdxB,
		int *colIdxB,
		bool tSparse,
		 option_type *global, utilsoption_type *utils,
		double *C,
		 Long Ldc);
void sparseTGeno(Uchar *code,
		 Long nrow_code,
		 Long ncol_code,
		 Long lda,
		 coding_type coding,
		 double *valueB,
		 int nIdx, // 2nd dimension of sparse; == length(rowIdxB)
		 int *rowIdxB,
		 int *colIdxB,
		 bool tSparse,
		 option_type *global, utilsoption_type *utils,
		 double *C,
		 Long Ldc);


#endif

