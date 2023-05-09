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



#ifndef miraculix_5codes_H
#define miraculix_5codes_H 1

#define gV5_header(NR, TYPE, SNDTYPE, TRDTYPE, ENDTYPE)			\
  void genoVector5v##NR##_##TYPE vector_header(TYPE, ENDTYPE)		\


#if defined  basic_miraculix_H

#define ORIGBITSperFIVE 10 // 5 * 2bits

Long LdaFiveCodes(Long snps, Long ldAbitalign);
//Long LdaFiveTransCodes(Long snps, Long ldAbitalign);
void Init5();


coding_header(_4Byte, 5);
coding_header(_1Byte, 5);
coding_header(_4Byte, 5trans);
coding_header(_1Byte, 5trans);

void coding5(unit_t *M, Long ldM,
	     Long start_snp, Long end_snp, 
	     Long start_individual, Long end_individual,     
	     basic_options *opt,  double * G, SEXP Ans);


gV_vG_header(LongDouble, 5);
gV_vG_header(double, 5);
gV_vG_header(Ulong, 5);


void trafo2Geno5codes32(unit_t *OLD, Long snps, Long individuals, Long oldLDA,
			coding_type coding,
			int cores, 
			unit_t *ANS, Long ldAns);
void trafo2Geno5codestrans32(unit_t *OLD, Long snps, Long individuals,
			     Long oldLDA, coding_type coding,
			     int cores, 
			     unit_t *ANS, Long ldAns);
void trafo2Geno5codes256(unit_t *OLD, Long snps, Long individuals, Long oldLDA,
			 coding_type coding,
			 int cores,
			 unit_t *ANS, Long ldAns);
void trafo2Geno5codestrans256(unit_t *OLD, Long snps, Long indiv, Long oldLDA,
			      coding_type coding,
			      int cores,
			      unit_t *ANS, Long ldAns);


gV5_header(512, floatD, LongDouble, float, double);
gV5_header(512, double, LongDouble, double, double);
gV5_header(256, floatD, LongDouble, float, double);
gV5_header(256, double, LongDouble, double, double);
gV5_header(128, floatD, LongDouble, float, double);
gV5_header(128, double, LongDouble, double, double);
gV5_header(32, floatD, LongDouble, float, double);
gV5_header(32, double, LongDouble, double, double);


#endif



#ifdef __cplusplus
extern "C" {
#endif
  
  
  void setOptions5(int gpu, int cores, int floatprecision,
		   int meanV, int meanSubstract,		   
		   int missingsFully0,
		   int centered, int normalized,
		   int use_miraculix_freq, 
		   int variant, int print_details);
  void plinkToSxI(char *plink, char *plink_transposed, 
		  long snps, long indiv, // long OK
		  int coding,
		  double *f, 
		  int max_n,
		  void**compressed);
  void free5(void **compressed) ;
  void genoVector5api(void *SxI, double *V,
		      long repetV, long LdV, double *ans, long LdAns);// long OK
  void vectorGeno5api(void *SxI, double *V,
		      long repetV, long LdV, double *ans, long LdAns);// long OK
  
  void getFreq5(void *compressed,  double *f) ;

  void getStartedOptions();
  
  
  void sparseTGenoPlinkApi(char *compressed,
			   int snps,
			   int indiv,
			   double *valueB,
			   int nIdx,// 2nd dimension of sparse; == length(rowIdxB)
			   int *rowIdxB,
			   int *colIdxB,
			   int tSp,
			   double *C,
			   int Ldc);
  
  void vectorGenoPlinkApi(char *compressed,
			   int snps, 
			   int indiv,
			   double *f,
			   double *B,
			   int n,
			   int LdB,
			   double *C,
			   int Ldc);
    
  
  void check_started();
  void plink2compressed(char *plink, char *plink_transposed,
			int snps, int indiv,
			double *f, 
			int max_n,
			void**compressed);
  void dgemm_compressed(char *t, // transposed?,
			void *compressed,
			int n, // number of columns of the matrix B
			//			double *f,
			double *B,
			int Ldb, 
			double *C, int Ldc);
  void free_compressed(void **compressed) ;
  
  void free5(void **compressed);
  
  void get_compressed_freq(void *compressed, double *f);
  
#ifdef __cplusplus
}
#endif


#endif
