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


#ifndef miraculix_Haplo_H
#define miraculix_Haplo_H 1



unit_t GetTwoBitHaplo(unit_t *code, Long snp);

typedef unit_t (*GetTwoBitHaploFun)(unit_t *code, Long snp);


void codeHaplo_4Byte(unit_t *MM, Long snps, Long individuals, Long lda, 	
		    coding_type coding,		
		    int *indiv_list, int indiv_list_N, 
		    basic_options *opt,
		    unit_t *Ans, Long NewIndiv, Long ldAns,
		    coding_type newCoding );

void codeHaplo_1Byte(unit_t *MM, Long snps, Long individuals, Long lda,    
		     coding_type coding, 			
		     int *indiv_list, int indiv_list_N,
		     basic_options *opt,
		     unit_t *Ans, Long NewIndiv,
		     Long ldAns, coding_type newCoding);

void decodeHaplo_4Byte(unit_t *code, Long snps, Long individuals ,Long lda,
		      coding_type coding, 		      
		      int *indiv_list, int indiv_list_N,
		      usr_bool set1, int cores,
		      unit_t *Ans, Long NewIndiv, Long ldAns,
		      coding_type newCoding);

void decodeHaplo_1Byte(unit_t *code, Long snps, Long individuals, Long lda,
 		      coding_type coding, 		       
		      int *indiv_list, int indiv_list_N,
		       usr_bool set1, int cores,
		       unit_t *Ans, Long NewIndiv, Long ldAns,
		       coding_type newCoding
		       );

void getHaploIncr(Long snps, Long individuals, Long lda, coding_type coding, 
		  Long compressed_lda, coding_type compressed_coding, 
		  Long currentBits, // 1 or 2
		  Long *uIncr, Long *nextHaploIncr,  // all un-compressed
		  Long *delta, Long *indivIncr,     // all un-compressed
		  Long *bitsPerCodeCompr, Long *shiftIncrCompr,//compressed
		  Long *deltaCompressed, Long *CpUcompressed);

void getHaploIncr(Long individuals, Long lda, coding_type coding, 
		  Long currentBits,
		  Long *nextHaploIncr, Long *delta, Long *indivIncr);


coding_header(_4Byte,Haplo1Byte);

void InnerDetermTwoBitHaplo(double *freq1, double *freq2,
			    Long snps,Long individuals,
			    coding_type coding,
			    unit_t *ans, Long ldAns);
void InnerRandomTwoBitHaplo(double *freq1, double *freq2,
			    Long snps,Long individuals,
			    coding_type coding,
			    unit_t *ans, Long ldAns);
#endif
