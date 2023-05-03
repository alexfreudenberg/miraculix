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



#ifndef miraculix_Bit23Uint_H
#define miraculix_Bit23Uint_H 1


#define UNIT_CODING(UINT, M, ldM, start_snp, end_snp, \
		    start_individual, end_individual, Ans, ldAns, FROM)	\
  UINT *pM = (UINT*) (M + (i - start_individual) * ldM);	\
  Uint shift = 0,						\
    counter = 0;							\
  unit_t compressed = 0;						\
  unit_t *pAns =  Ans + i * ldAns;					\
  for (Long s=start_snp; s<end_snp; s++) {				\
    compressed |= (unit_t) (geno_code[FROM] << shift);			\
    shift += BitsPerCode;						\
    if (++counter >= CodesPerPartUnit) {				\
      shift += deltaBitsPartUnit; /* 3bit */				\
      counter = 0;							\
    }									\
    if (shift >= BitsPerUnit) {						\
      *pAns = compressed;						\
      pAns++;								\
      shift = counter = 0;						\
      compressed = ZERO();						\
    }									\
  }									\
  if (shift > 0) {							\
    *pAns = compressed;							\
    pAns++;								\
  }									




typedef  char table_type; //  10.585 
//typedef unsigned char table_type; //  10.612
//typedef unsigned int table_type; //   10.325 
//typedef int table_type; //  10.217


void initiate_tableI(table_type **TABLE, int TABLESIZE,
		     int codesPerPartUnit, int bits,
		     Uint *result_code,
		     Uint *result_value, int NrResults);

inline static Uint ExtractCode(unit_t *M, Long S) {
  const unit_t MC = M[S / CodesPerUnit];
  const int s = (int) (S % CodesPerUnit);
  return (MC >> (BitsPerPartUnit * (s / CodesPerPartUnit) +
		 (s % CodesPerPartUnit) * BitsPerCode)) & CodeMask;
}

inline static void PositionCode(Long S,  Long *j, unit_t *mask) {
  *j = S / CodesPerUnit;
  int idx = (int) S % CodesPerUnit;	
  *mask = (unit_t) CodeMask << (BitsPerPartUnit * (idx / CodesPerPartUnit) +
				(idx % CodesPerPartUnit) * BitsPerCode);
}

/*
inline static void PositionCode(Long S, Uint codemask,  Long *j, Uint *mask) {
  *j = S / CodesPerUnit;
  int idx = (int) (S % CodesPerUnit);	
  *mask = codemask << ( BitsPerPartUnit * (idx / CodesPerPartUnit) +
			(idx % CodesPerPartUnit) * BitsPerCode);
}
*/



/* 2 Bit / Plink Transposed */
/* EXTREM LANGSAM !! */
#define coding2trans_start(UINT,NAME)				\
  coding_header(UINT, NAME##trans) {					\
  if (start_snp != 0) BUG;						\
  if (start_individual % CodesPerBlock != 0) BUG; /* // OK */		\
  unit_t *ans =Ans + start_individual / CodesPerUnit;			\
  Long snps = end_snp - start_snp;					\
  Long deltaIndiv = end_individual - start_individual;			\
  Long end_individualM4 = (deltaIndiv) / 4;				\
  Long ldAnsInByte = ldAns * sizeof(unit_t);


#define coding2trans_do(UINT,TRAFO)						\
  for (Long i=0; i<end_individualM4; i++) {				\
    Uchar *a = ((Uchar*) ans) + i;					\
    UINT *pM0 = (UINT*) (M + (4 * i + 0) * ldM);		\
    UINT *pM1 = (UINT*) (M + (4 * i + 1) * ldM);		\
    UINT *pM2 = (UINT*) (M + (4 * i + 2) * ldM);		\
    UINT *pM3 = (UINT*) (M + (4 * i + 3) * ldM);		\
    for (Long j=0; j<snps; j++) {					\
      a[ldAnsInByte*j] = (Uchar)(TRAFO(pM0[j]) + (TRAFO(pM1[j]) << 2) + \
				 (TRAFO(pM2[j]) << 4) + (TRAFO(pM3[j]) << 6)); \
    }									\
  }									\
  Long i = end_individualM4;						\
  Uchar *a = ((Uchar*) ans) + i;					\
  UINT *pM0 = (UINT*) (M + (4 * i + 0) * ldM);		\
  UINT *pM1 = (UINT*) (M + (4 * i + 1) * ldM);		\
  UINT *pM2 = (UINT*) (M + (4 * i + 2) * ldM);		\
  if (4 * i < deltaIndiv) {						\
    for (Long j=0; j<snps; j++) {					\
      Uchar a0 = (Uchar) TRAFO(pM0[j]);					\
      if (4 * i + 1 < deltaIndiv) {					\
	a0 += (Uchar) (TRAFO(pM1[j]) << 2);				\
	if (4 * i + 2 < deltaIndiv) {					\
	  a0  += (Uchar) (TRAFO(pM2[j]) << 4);				\
	}								\
      }									\
      a[ldAnsInByte * j] = a0;						\
      if (false) PRINTF("a[%ld, %ld] = %d; %ld %ld ldAnsInByte=%ld\n", i,j, a0, i, end_individualM4, ldAnsInByte); \
      /*    if (a[ldAnsInByte * j]) {printf("X %d %d %d %d; %ld %ld %ld %ld\n", (int)a[ldAnsInByte * j], (int)pM0[j],(int)pM1[j],(int) pM2[j],j, snps, deltaIndiv, 4 * i);  BUG;} */  \
    }									\
  }									\
  }



#endif
