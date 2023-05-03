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

#define MY_VARIANT 64
#define PlainInteger64 1

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "Template.h"
#include "utils_miraculix.h"
#include "MXinfo.h"
#include "Haplo.h"
#include "transform.h"

void getHaploIncr(Long individuals, Long lda, coding_type coding, 
		  Long currentBits,
		  Long *nextHaploIncr, // in [UINT]
		  Long *delta,         // in [unit_t] 
		  Long *indivIncr      // in [unit_t]    
		  ) {
  switch(coding) {
  case OneByteHaplo:
  case FourByteHaplo:
  case OneByteGeno :
  case FourByteGeno : 
    *indivIncr = lda;    // next individual in [unit_t]    
    *nextHaploIncr = 1L; // next haplo/geno;  in [UINT](!) steps
    *delta = lda *individuals;//jump to twin bit in [unit_t]
      break;
 case FourByteSingleBit:
    *indivIncr = lda * currentBits;
    *nextHaploIncr = 1L; // == lda / sizeof(UINT); // 
    *delta = lda  ;
     break;
  case EightByteHaplo: 
    *indivIncr = lda;
    *nextHaploIncr = currentBits;
    *delta = 1L;
    break;
  case FourByteSingleBit_Tr:
  case FourByte_Tr  :
    *indivIncr = currentBits;
    *nextHaploIncr = currentBits * individuals;
    *delta = 1L;
    break;
  case EightByteHaplo_Tr: 
    *indivIncr = 1L;
    *nextHaploIncr = currentBits * individuals;
    *delta = lda;
    break;
  default: BUG;
  }
  *delta *= (currentBits > 1);
}

void getHaploIncr(Long VARIABLE_IS_NOT_USED snps, Long individuals, Long lda,
		  coding_type coding, 
		  Long compressed_lda, coding_type compressed_coding,
		  Long currentBits, // 1 or 2
		  Long *unitIncr, Long *nextHaploIncr,  // all un-compressed
		  Long *delta, Long *indivIncr,     // all un-compressed
		  Long *bitsPerCodeCompr, Long *shiftIncrCompr,// compressed
		  Long *deltaCompressed, Long *CpUcompressed) {// compressed
  //  printf("entering gethapoincr\n");

  assert(compressed_coding == OneBitHaplo ||
	 compressed_coding == TwoBitHaplo ||
	 compressed_coding == OneByteHaplo);
  bool doubled_cols = doubledCols(compressed_coding, true);
  *CpUcompressed = GetCodesPerUnit(compressed_coding);
  *bitsPerCodeCompr = BitsPerUnit /  *CpUcompressed;
  *shiftIncrCompr = (*bitsPerCodeCompr / 2) * (!doubled_cols);
  *deltaCompressed = compressed_lda * individuals * doubled_cols;

  if (coding == UnknownSNPcoding) {
    //printf("RETURNING\n");
    assert(unitIncr == NULL);
    return;
  }
  assert(unitIncr != NULL);

  getHaploIncr(individuals, lda, coding, currentBits, // in 
	       nextHaploIncr, delta, indivIncr); //out

  switch(coding) {
  case OneByteHaplo:
  case FourByteHaplo:
  case OneByteGeno :
  case FourByteGeno :
    // stepforward [unit_t] after a unit_t has 
    //      been completely filled with compressed code
    *unitIncr = *CpUcompressed / GetCodesPerUnit(coding);
    break;
  case FourByteSingleBit:
    *unitIncr = *CpUcompressed / GetCodesPerUnit(coding);
    break;
  case EightByteHaplo: 
    *unitIncr = *CpUcompressed * currentBits;
    break;
  case FourByteSingleBit_Tr:
  case FourByte_Tr  :
    *unitIncr = *CpUcompressed * *nextHaploIncr;
    break;
  case EightByteHaplo_Tr:
    *unitIncr = *CpUcompressed * *nextHaploIncr;
    break;
  default: BUG;
  }
  
}


// Achtung!! BigEndian und LitteEndian: Uchar[] und 4-Byte-int
// sind nicht dasselbe!!. Somit kann Uchar NICHT gleichzeitig als
// internes (SSE) und externes (array) format interpretiert werden!!
// FESTLEGUNG: Uchar[] ist das fuehrende Format; soit wird bei
// interpretation 
#define codeHaplo_start(UINT)						\
  void codeHaplo##UINT(unit_t *MM, Long snps, Long OldIndiv, Long lda,	\
			coding_type coding,				\
			int *indiv_list, int indiv_list_N,		\
		        basic_options *opt,				\
			unit_t *Ans, Long NewIndiv, Long ldAns,		\
			coding_type newCoding/*2BitHaplo*/) {					\
  /* // printf("coding %s -> %s\n", CODING_NAMES[coding], CODING_NAMES[newCoding]); */  \
  Long unitIncr, nextHaploIncr, delta, indivIncr, shiftIncrCompr,	\
    bitsPerCodeCompr, deltaCompressed, CpU;				\
  Long  indiv_N = indiv_list_N > 0 ? NewIndiv : OldIndiv,		\
    nIndiv = 0;								\
  Long  total = 0;							\
  getHaploIncr(snps, OldIndiv, lda, coding, 				\
	       ldAns, newCoding, 2,					\
	       &unitIncr, &nextHaploIncr, &delta, &indivIncr,		\
	       &bitsPerCodeCompr, &shiftIncrCompr, &deltaCompressed, &CpU); \
  stopIfNotInt(bitsPerCodeCompr);					\
  int BitsPCCompr = (int) bitsPerCodeCompr;				\
	/* // printf("code lda=%u indiv=%u uincr=%lu nxtHaploIncr=%lu delta=%lu nxtIndiv=%lu\n\tldAns=%u bpcc=%lu shift=%lu deltaC=%lu CpU=%lu sizeUInt=%lu\n", \
	 lda, OldIndiv,						\
	 unitIncr, nextHaploIncr, delta, indivIncr,			\
	 ldAns,							\
	 BitsPCCompr, shiftIncrCompr, deltaCompressed, CpU, sizeof(UINT)); */  \
  Long allUnitsM1 = snps2entities(snps, CpU) - 1;			\
  int bigendianCorr = BitsPerUnit - BitsPCCompr;		\
  int rest = (int) ((snps - allUnitsM1 * CpU) * BitsPCCompr);	\
  bool bigendian = opt->bigendian & isOneByte(newCoding);	\
  MEMSET(Ans, 0, totalMem(snps, NewIndiv, ldAns, newCoding, false,	\
			  true) * BytesPerUnit)

#define codeHaplo_ende(UINT)						\
  for (Long i=0; i<indiv_N; i++) { /*// printf("-- %lu %u\n",i, indiv_N);*/ \
    Long ii=i;								\
    if (indiv_list_N > 0) {						\
      ii = indiv_list[i];						\
      if (ii >= OldIndiv) continue;					\
    }									\
    /*// printf("%d %d incr=%d %lu\n", i, ii, indivIncr, MM); */	\
    unit_t *pC = Ans + (nIndiv++) * ldAns;				\
    unit_t *M = MM + ii * indivIncr;					\
    UINT check = 0;							\
    for (Long j=0; j<=allUnitsM1; j++) {				\
      unit_t dummy1 = 0;							\
      unit_t dummy2 = 0;							\
      UINT *mm1 = (UINT*) ( M + j * unitIncr);			\
      UINT *mm2 = (UINT*) (((unit_t*) mm1) + delta); /* // ((( OK */	\
      int end = j==allUnitsM1 ? rest : BitsPerUnit;		\
      for (int shift=0; shift < end;					\
	   mm1+=nextHaploIncr, mm2 += nextHaploIncr) {			\
	Long Shft = bigendian ? bigendianCorr - shift : shift;		\
	check |= *mm1 | *mm2;						\
	dummy2 |= ((unit_t) (*mm2)) << (Shft + shiftIncrCompr);		\
 	dummy1 |= ((unit_t) (*mm1)) << Shft;				\
	shift += BitsPCCompr;					\
      }									\
      if (shiftIncrCompr) {						\
	pC[j] = dummy1 | dummy2;  /* // printf("\n");BUG; */		\
      } else {								\
	pC[j] = dummy1;							\
	pC[j + deltaCompressed] = dummy2;				\
      }									\
    }									\
    total += check > 1;							\
    if (total > 0)	{ /*//printf("B\n");*/ BUG;}			\
  }									\
  if (total > 0) ERR0("values outside 0,1");				\
  }



codeHaplo_start(_4Byte);
//#ifdef DO_PARALLEL
//#pragma omp parallel for num_threads(cores) reduction(+:total) schedule(static)
//#endif
codeHaplo_ende(_4Byte)

codeHaplo_start(_1Byte);
//#ifdef DO_PARALLEL
//#pragma omp parallel for num_threads(cores) schedule(static)
//#endif
codeHaplo_ende(_1Byte)
   

void decodeHaplo_4Byte(unit_t *code, Long snps, Long OldIndiv, Long lda,	
			coding_type coding,			
		       int *indiv_list, int indiv_list_N, 
			usr_bool usr_set1,				
			int VARIABLE_IS_NOT_USED cores, 
			unit_t *Ans, Long NewIndiv, Long ldAns,
		      coding_type newCoding) {
  //  printf("hier %lu\n", (Long) code);
  const bool set1 = usr_set1 == Nan || usr_set1 == True;		
  Long unitIncr, nextHaploIncr, delta, indivIncr, shiftIncrCompr,	
    bitsPerCodeCompr, deltaCompressed, CpUlong;				
  Long nIndiv = 0,							
    currentBits =  1 + (usr_set1 == Nan),				
    indiv_N = indiv_list_N > 0 ? NewIndiv : OldIndiv;
  Long totalUnits = totalMem(snps, NewIndiv, ldAns, newCoding,
				 usr_set1 != Nan); 
  if ((currentBits == 1) xor isOneSet2G(newCoding))  BUG;		
  getHaploIncr(snps, OldIndiv, ldAns, newCoding,
	       lda, coding,			
	       currentBits,						
	       &unitIncr, &nextHaploIncr, &delta, &indivIncr,		
	       &bitsPerCodeCompr, &shiftIncrCompr, &deltaCompressed, &CpUlong);
  int CpU = (int) CpUlong;
  //printf("Ans=%lu\n", (Long) Ans);
  /*//printf("decode lda=%u indiv=%u uincr=%lu %lu delta=%lu %lu bpcc=%lu %lu deltaC=%lu %lun", 
	 lda, OldIndiv,							
	 unitIncr, nextHaploIncr, delta, indivIncr,			
	 bitsPerCodeCompr, shiftIncrCompr, deltaCompressed, CpU); */	
  MEMSET(Ans, 0, BytesPerUnit * totalUnits); /* avoids unset memory */	
  Long allUnitsM1 = snps2entities(snps, CpU) - 1;			
  int rest = (int) (snps - allUnitsM1 * CpU);					
  /* !! wg nIndiv: KEIN  # p r a g m a !!  */				
  for (Long i=0; i<indiv_N; i++) {  /*//printf("-- %lu %un",i, indiv_N);*/ 
    //    printf("i=%d\n", i);						
    Long ii = i;							
    if (indiv_list_N > 0) {						
      ii = indiv_list[i];						
      if (ii >= OldIndiv) continue;					
    }									
    unit_t *pC = code + lda * ii;						
    unit_t *M = Ans + (nIndiv++) * indivIncr;		
    // printf("M=%lu\n", (Long) M);
    for (Long j=0; j<=allUnitsM1; j++) {				
      // printf("i=%d (%d) %d %d (%d)\n", i, indiv_N, ii, j, allUnitsM1);	
      unit_t *mm1 =  M + j * unitIncr;			
      unit_t *mm2 =  mm1 + delta;			
      unit_t C1 = pC[j],							
	C2 = pC[j + deltaCompressed];
      int
	shift = 1,							
	end = j == allUnitsM1 ? rest : CpU;				
      //   printf("XXXC1=%u %u pC=%lu %lu\n", C1, C2, (Long) pC, (Long) (pC+deltaCompressed)); 
   for (int u=0; u<end; u++, mm1+=nextHaploIncr, mm2 += nextHaploIncr) { 
	/* ordering important !! as mm2 might be overwritten by mm!! */ 
     // printf("unitInce=%d\n", unitIncr);	
     //	printf("X %lu\n", (Long) mm1);	
	//    printf("ZX %lu\n", (Long) mm2);	
        *mm2 = (C2 & (shift << shiftIncrCompr)) > 0; /* never change as mm */ 
	                                        /* could be mm2!! */
       	if (set1) *mm1 = (C1 & shift) > 0;				
	/*// printf("(%u, %u) ", *mm1, *mm2);	*/			
	shift <<= bitsPerCodeCompr;					
   }									
   //  printf("half done\n");				
    }									
   } 								
  }





#define decodeHaplo(UINT)						\
  void decodeHaplo##UINT(unit_t *code, Long snps, Long OldIndiv, Long lda, \
			coding_type coding,		\
			int *indiv_list, int indiv_list_N, \
			usr_bool usr_set1,				\
			int VARIABLE_IS_NOT_USED cores, \
			unit_t *Ans, Long NewIndiv, Long ldAns,	\
			coding_type newCoding) {			\
  const bool set1 = usr_set1 == Nan || usr_set1 == True;		\
  Long unitIncr, nextHaploIncr, delta, indivIncr, shiftIncrCompr,	\
    bitsPerCodeCompr, deltaCompressed, CpUlong;				\
  Long nIndiv = 0,							\
    currentBits =  1 + (usr_set1 == Nan),				\
    indiv_N = indiv_list_N > 0 ? NewIndiv : OldIndiv;			\
  Long totalUnits = totalMem(snps, NewIndiv, ldAns, newCoding,		\
			     usr_set1 != Nan);				\
  if ((currentBits == 1) xor isOneSet2G(newCoding))  BUG;		\
  getHaploIncr(snps, OldIndiv, ldAns, newCoding,			\
	       lda, coding,						\
	       currentBits,						\
	       &unitIncr, &nextHaploIncr, &delta, &indivIncr,		\
	       &bitsPerCodeCompr, &shiftIncrCompr, &deltaCompressed, \
	       &CpUlong);						\
  MEMSET(Ans, 0, BytesPerUnit * totalUnits); /* avoids unset memory */	\
  int CpU = (int) CpUlong;						\
  Long allUnitsM1 = snps2entities(snps, CpU) - 1;			\
  int rest = (int) (snps - allUnitsM1 * CpU);				\
  /* !! wg nIndiv: KEIN  # p r a g m a !!  */				\
  for (Long i=0; i<indiv_N; i++) {  /*//printf("-- %lu %u\n",i, indiv_N);*/ \
    /*//printf("i=%lu\n", i);*/						\
    Long ii = i;							\
    if (indiv_list_N > 0) {						\
      ii = indiv_list[i];						\
      if (ii >= OldIndiv) continue;					\
    }									\
    unit_t *pC = code + lda * ii;						\
    unit_t *M =  Ans + (nIndiv++) * indivIncr;		\
    for (Long j=0; j<=allUnitsM1; j++) {				\
      /*// printf("i=%lu (%d) %lu %lu (%d)\n", i, indiv_N, ii, j, allUnitsM1); */ \
      UINT *mm1 = (UINT*) (M + j * unitIncr);			\
      UINT *mm2 =  (UINT*) (((unit_t*) mm1) + delta);  /* // ((( OK */	\
      unit_t C1 = pC[j],						\
	C2 = pC[j + deltaCompressed];					\
      int shift = 1,							\
	end = j == allUnitsM1 ? rest : CpU;				\
      /*//printf("C1=%u %u pC=%lu %lu\n", C1, C2, (Long) pC, (Long) (pC+deltaCompressed));*/ \
   for (int u=0; u<end; u++, mm1+=nextHaploIncr, mm2 += nextHaploIncr) { \
	/* ordering important !! as mm2 might be overwritten by mm!! */ \
	*mm2 = (C2 & (shift << shiftIncrCompr)) > 0; /* never change as mm */ \
	                                        /* could be mm2!! */	\
	if (set1) *mm1 = (C1 & shift) > 0;				\
	/*// printf("(%u, %u) ", *mm1, *mm2);	*/			\
	shift <<= bitsPerCodeCompr;					\
   }									\
   /*//printf("half done\n");*/						\
    }									\
   } 								\
  }



      
      //decodeHaplo(_4Byte);
decodeHaplo(_1Byte)




#define coding_Haplo1Byte_start(UINT)				\
  coding_header(UINT,Haplo1Byte) {				\
  Long delta = ldAns * individuals;					\
  Long DeltaM = ldM * (end_individual - start_individual + 1);	

#define coding_Haplo1Byte_end(UINT)				\
  /* // printf("coding_Haplo1Byte\n"); */			\
  for (Long i=start_individual; i<end_individual; i++) {		\
    UINT *pM = (UINT*) (M + (i - start_individual) * ldM);	\
    UINT *pM2 = (UINT*) (((unit_t*) pM) + DeltaM); /* // ((( OK */	\
    unit_t *A = Ans + i * ldAns;					\
    _1Byte *pAns = (_1Byte*) A;						\
    _1Byte *pAns2 = (_1Byte*) (A + delta);				\
    for (Long s=start_snp; s<end_snp; pM++, pM2++, s++) {		\
      pAns[s] = (_1Byte) *pM;						\
      pAns2[s] = (_1Byte) *pM2;						\
    }									\
  }									\
}

coding_Haplo1Byte_start(_4Byte)				     
//#ifdef DO_PARALLEL
  //#pragma omp parallel for num_threads(cores) schedule(static)
//#endif
coding_Haplo1Byte_end(_4Byte)	

coding_Haplo1Byte_start(_1Byte) // geht schneller nur mit MEMCOPY -- ToDo 
//#ifdef DO_PARALLEL
  //#pragma omp parallel for num_threads(cores) schedule(static)
//#endif
coding_Haplo1Byte_end(_1Byte)

