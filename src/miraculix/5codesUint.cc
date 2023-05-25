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
// sudo tomte reconfigure all



#include "5codesDef.h"
#define MY_VARIANT 32


#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "TemplateUint.h"
#include "Haplo.h"
#include "MX.h"
//#include "1bit.h"
#include "5codes.h"
#include "5codesIntern.h"





Uchar CODING_TABLE5[3][3][3][3][3],
  TWOBIT2FIVE[1024], PLINK2FIVE[1024];
Uint ALLELESUM[243];

Long LdaFiveCodes(Long snps, Long VARIABLE_IS_NOT_USED LDAbitalign) {
  //  printf("snps=%ld %ld %ld\n", snps, LDAalignBitSnps(snps * BitsPerByte), MY_LDABITALIGN);  if (LDAbitalign) BUG;
  return LDAalignBitSnps(snps * BitsPerByte); 
}



void initiate_table5I() { // hash table  
  for (Uint i4=0; i4<3; i4++) {
    Uint V4 = i4;
    for (Uint i3=0; i3<3; i3++) {
      Uint V3 = 3U * V4 + i3;
      for (Uint i2=0; i2<3; i2++) {
	Uint V2 = 3U * V3 + i2;
	for (Uint i1=0; i1<3; i1++) {
	  Uint V1 = 3U * V2 + i1;
	  for (Uint i0=0; i0<3; i0++) {
	    Uint V0 = 3U * V1 + i0;
	    CODING_TABLE5[i0][i1][i2][i3][i4] = (Uchar) V0;
	    ALLELESUM[V0] = i4 + i3 + i2 + i1 + i0;
	  }
	}
      }
    }
  }

 for (Uint i4=0; i4<4; i4++) {
    Uint V4 = i4;
    Uint twobit4 = TWOBIT2HUMAN(i4);
    Uint plink4 = PLINK2HUMAN(i4);
    for (Uint i3=0; i3<4; i3++) {
    Uint V3 = 4U * V4 + i3;
      Uint twobit3 = 3U * twobit4 + TWOBIT2HUMAN(i3);
      Uint plink3 = 3U * plink4 + PLINK2HUMAN(i3);
    for (Uint i2=0; i2<4; i2++) {
	Uint V2 = 4U * V3 + i2;
	Uint twobit2 = 3U * twobit3 + TWOBIT2HUMAN(i2);
	Uint plink2 = 3U * plink3 + PLINK2HUMAN(i2);
 	for (Uint i1=0; i1<4; i1++) {
	  Uint V1 = 4U * V2 + i1;
	  Uint twobit1 = 3U * twobit2 + TWOBIT2HUMAN(i1);
	  Uint plink1 = 3U * plink2 + PLINK2HUMAN(i1);
	  for (Uint i0=0; i0<4; i0++) {
	    Uint V0 = 4U * V1 + i0;
	    Uint twobit0 = 3U * twobit1 + TWOBIT2HUMAN(i0);
	    Uint plink0 = 3U * plink1 + PLINK2HUMAN(i0);
	    TWOBIT2FIVE[V0] = (Uchar) twobit0;
	    PLINK2FIVE[V0] = (Uchar) plink0;
	  }
	}
    }
    }
 }
}


void Init5() {
  //  assert(PartUnitsPerBlock == 2); // OK
  //printf("CpB %f %f\n", (double) CpB, (double) CodesPerByte); 
  assert(CpB == CodesPerByte);
  initiate_table5I();
}


#define coding5_start(UINT)						\
  coding_header(UINT, 5) {						\
  if (start_snp != 0) BUG;						\
  if (start_individual != 0) BUG;					\
  unit_t *ans = Ans + start_snp / BytesPerUnit;				\
  Long end_indiM1_CpB = (end_individual - start_individual - 1L) / CpB, \
    ldMUINT = ldM * BytesPerUnit / sizeof(UINT);			


#define coding5_do(UINT)						\
  for (Long i=0; i<=end_indiM1_CpB; i++) {				\
    Long rest = CpB;							\
    UINT *pM = (UINT*) (M + ldM * CpB * i);				\
    Uchar *a = ((Uchar*) (ans + ldAns * i));				\
    Uint x[CpB] = {0, 0, 0, 0, 0};					\
    if (i < end_indiM1_CpB) {} else {	 				\
      rest = (end_individual - start_individual) - end_indiM1_CpB * rest; \
      assert(rest > 0);							\
      for (Long k=rest; k<CpB; k++) x[k] = 0;				\
    }									\
     for(Long j=start_snp; j<end_snp; j++, pM++) {			\
      for (Long k=0; k<rest; k++) {	 				\
	x[k] = (Uint) pM[k * ldMUINT];					\
      }									\
      assert(a + j < (Uchar*) (ans + ldAns * individuals));		\
      a[j] = (Uchar) CODING_TABLE5[x[0]][x[1]][x[2]][x[3]][x[4]];	\
      /* a[i] = (Uchar) (x[0] + 3 * x[1] + 9 * x[2] + 27 * x[3] + 81 * x[4]); */ \
    }									\
  }									\
  }



coding5_start(_4Byte)
#ifdef DO_PARALLEL 
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding5_do(_4Byte)


coding5_start(_1Byte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding5_do(_1Byte)




/* Five Codes Transposed */
/* EXTREM LANGSAM !! */
#define coding5trans_start(UINT)			\
  coding_header(UINT, 5trans) {				\
  if (start_snp != 0) BUG;				\
  if (start_individual != 0) BUG;			\
  unit_t *ans =Ans + start_snp / CodesPerUnit;		\
  Long end_snpM1_CpB = (end_snp - start_snp -1L) / CpB,	\
    ldaByte = ldAns * BytesPerUnit;


#define coding5trans_do(UINT)						\
  for (Long i=start_individual; i<end_individual; i++) {		\
    Uchar *a0 = ((Uchar*) ans) + i;					\
    UINT *pM0 = (UINT*) (M + i * ldM);					\
    for(Long j=0; j<=end_snpM1_CpB; j++) {				\
      Long rest = CpB;							\
      UINT *pM = pM0 + CpB * j;						\
      Uchar *a = a0 + ldaByte * j;					\
      Uint x[CpB] = {0, 0, 0, 0, 0};					\
      if (j < end_snpM1_CpB) {} else {					\
	rest = (end_snp - start_snp) - end_snpM1_CpB * rest;		\
	assert(rest > 0);						\
	for (Long k=rest; k<CpB; k++) x[k] = 0;				\
      }									\
      for (Long k=0; k<rest; k++) {					\
	x[k] = (Uint) pM[k]; /* HIER SCHEINEN SICH DIE HINTEREINANDERLIEGENDEN BYTES IM WEG ZU STEHEN */ \
      }									\
      *a = (Uchar) CODING_TABLE5[x[0]][x[1]][x[2]][x[3]][x[4]];		\
    }									\
  }									\
  }




coding5trans_start(_4Byte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif								
coding5trans_do(_4Byte)

coding5trans_start(_1Byte)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
coding5trans_do(_1Byte)


void coding5(unit_t *M, Long ldM,
	     Long start_snp, Long end_snp, 
	     Long start_individual, Long end_individual,     
	     basic_options *opt,  double * G, SEXP Ans) {
  int cores =  GreaterZero(opt->cores);
  Long *info =  GetInfo(Ans);  
  coding5_4Byte(M, ldM,
		start_snp, end_snp, start_individual, end_individual, 
		cores, G,
		Align(Ans, opt), info[INDIVIDUALS], info[LDA]);
  SEXP next = getAttribPointer(Ans, Next);
  assert(next != R_NilValue);
  Long *infoNext =  GetInfo(next);  
  coding5trans_4Byte(M, ldM,
		     start_snp, end_snp, 
		     start_individual,  end_individual,	
		     cores, G,
		     Align(next, opt), infoNext[INDIVIDUALS], infoNext[LDA]);
}          


gV5_start(32, longD, LongDouble, LongDouble, LongDouble, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_CreateHash(longD, LongDouble, LongDouble, LongDouble, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_MainLoop(longD, LongDouble, LongDouble, LongDouble, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_SumUp(longD, LongDouble, LongDouble, LongDouble, double)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_Sort(longD, LongDouble, LongDouble, LongDouble, double)



gV5_start(32, LongDouble, LongDouble, LongDouble, LongDouble, LongDouble)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_CreateHash(LongDouble, LongDouble, LongDouble, LongDouble, LongDouble)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif 
gV5_MainLoop(LongDouble, LongDouble, LongDouble, LongDouble, LongDouble)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_SumUp(LongDouble, LongDouble, LongDouble, LongDouble, LongDouble)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_Sort(LongDouble, LongDouble, LongDouble, LongDouble, LongDouble)


   
gV5_start(32, Ulong, Ulong, Ulong, Ulong, Ulong)
#ifdef DO_PARALLEL 
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif 
gV5_CreateHash(Ulong, Ulong, Ulong, Ulong, Ulong) 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif 
gV5_MainLoop(Ulong, Ulong, Ulong, Ulong, Ulong)
#ifdef DO_PARALLEL 
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif 
gV5_SumUp(Ulong, Ulong, Ulong, Ulong, Ulong)
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif
gV5_Sort(Ulong, Ulong, Ulong, Ulong, Ulong)




#define gV_vG5_start(TYPE)	 					\
  gV_vG_start(TYPE, 5);							\
  if (!gV) {								\
    if (!global->tuning.addtransposed) {				\
      ERR0("option 'addtransposed' must be true.");			\
    } else BUG;								\
  }

#define gV5_Select(TYPE)			\
  gV_vG5_start(TYPE);							\
 /* vD = variant < 256 ? genoVector5v32_##TYPE : genoVector5v256_##TYPE;*/ \
  vD = genoVector5v32_##TYPE;						\
  gV_vG_end;    		 					\
  } 
  
gV5_Select(LongDouble)     
gV5_Select(Ulong)  
 
  
gV_vG5_start(double);  {
  // variant = 32;  printf("variant v%d of 5codes with loop=%d\n", variant, tuning->floatLoop);
  if (tuning->floatLoop < 0) vD = genoVector5v32_longD;	
  else if (tuning->floatLoop == 0) {
    vD = (variant < 128 ? genoVector5v32_double : 
	  variant < 256 ? genoVector5v128_double :
	  variant < 512 ? genoVector5v256_double :
	  genoVector5v512_double);
  }
  else vD = (variant < 128 ? genoVector5v32_floatD :
	     variant < 256 ? genoVector5v128_floatD :
	     variant < 512 ? genoVector5v256_floatD :
	     genoVector5v512_floatD);
  gV_vG_end;
  //  printf("end gV_vG5 longD=%d gv=%d\n", vD == genoVector5v32_longD, gV);
}}  

 
#define HANGING_TYPE Uint
#define TBM 0x00000003
void trafo2Geno5codes32(unit_t *M, Long snps, Long individuals, Long ldM,
			coding_type coding,
			int VARIABLE_IS_NOT_USED cores,
			unit_t *Ans, Long ldAns) {
  if (sizeof(HANGING_TYPE) * BitsPerByte != MY_VARIANT) BUG;
  Uchar *table = coding == TwoBitGeno ? TWOBIT2FIVE
      : coding == Plink ? PLINK2FIVE : NULL;
  if (table == NULL) BUG;
  Long
    end_indiM1_CpB = (individuals - 1L) / CpB,
    CpUorig = GetCodesPerUnit(coding),
    Munits = 1L + (snps - 1L) / CpUorig;

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<=end_indiM1_CpB; i++) {
    Long rest = CpB;
    unit_t *pM =  M + ldM * CpB * i;
    Uchar *a0 = ((Uchar*) (Ans + ldAns * i ));
    HANGING_TYPE h[CpB];
 
    if (i < end_indiM1_CpB) {} else {
       rest = individuals - end_indiM1_CpB * rest;		
      assert(rest > 0);
      for (Long ll=rest; ll<CpB; ll++) h[ll] = 0;
     }
    
    for (Long j=0; j<Munits; j++, pM++) {
      Uchar *a = a0 + j * CpUorig;
      for (Long z=0; z<rest; z++) h[z] = (HANGING_TYPE) pM[z * ldM];
      for (Long k=0; k<CpUorig; k++) {
	const HANGING_TYPE m = (h[0] & TBM) | ((h[1] & TBM) << 2) |
	  ((h[2] & TBM) << 4) | ((h[3] & TBM) << 6) | ((h[4] & TBM) << 8);
	a[k] = table[m];
	for (Long ll=0; ll<CpB; ll++) h[ll] >>= 2;
      }
    }
  } 
}



#define ORIGBITS_MASK 0x000003FF
#define LONG_HANGING_TYPE Uint
void trafo2Geno5codestrans32X(unit_t *M, Long snps, Long individuals, Long ldM,
			      coding_type coding,
			      int VARIABLE_IS_NOT_USED cores,
			      unit_t *Ans, Long ldAns){
  Uchar *table = coding == TwoBitGeno ? TWOBIT2FIVE
    : coding == Plink ? PLINK2FIVE : NULL;
  if (table == NULL) BUG;
  Long
    end_snpM1_CpB = (snps -1L) / CpB,			
    ldaByte = ldAns * BytesPerUnit;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) {
    int availableBits = 0;
    LONG_HANGING_TYPE hanging = 0,
      m = 0;
    int hangingBits = 0;
    Uchar *a0 = ((Uchar*) Ans) + i;					
    unit_t *pM = M + i * ldM;
    for(Long j=0; j<=end_snpM1_CpB; j++) {      
      Uchar *a = a0 + ldaByte * j;					
      while (availableBits < ORIGBITSperFIVE) { 
	if (hangingBits == 0) {
	  hanging =  (HANGING_TYPE) *(pM++);
	  hangingBits = MY_VARIANT;
	}
	m |= hanging << availableBits;
	int  movedBits =  hangingBits <= MY_VARIANT - availableBits
	  ? hangingBits : MY_VARIANT - availableBits;
	availableBits += movedBits;
	hangingBits -= movedBits;
	hanging >>= movedBits; 
      }
      *a = table[m & ORIGBITS_MASK];
      m >>= ORIGBITSperFIVE;
      availableBits -= ORIGBITSperFIVE;
    }	
  }	
}


#if defined  LONG_HANGING_TYPE
#undef LONG_HANGING_TYPE
#endif
//#define LONG_HANGING_TYPE Long
//#define HTT_BITS 64
#define LONG_HANGING_TYPE Long
#define HTT_BITS 64
#define ORIGBITS_TRANS_MASK 0x00000000000003FF
// irgendwo fehler drin -- sollte eigentlich bisschen schneller gehen
void trafo2Geno5codestrans32(unit_t *M, Long snps, Long individuals, Long ldM,
			     coding_type coding,
			     int VARIABLE_IS_NOT_USED cores,
			     unit_t *Ans,  Long ldAns) {
  Uchar *table = coding == TwoBitGeno ? TWOBIT2FIVE
    : coding == Plink ? PLINK2FIVE : NULL;
  if (table == NULL) BUG;
  unit_t *ans =Ans; /* */
  Long
    end_snpM1_CpB = (snps -1L) / CpB,			
    ldaByte = ldAns * BytesPerUnit;

  if (MY_LDABITALIGN_2BIT < sizeof(Long)) BUG;
  if(HTT_BITS != sizeof(LONG_HANGING_TYPE) * BitsPerByte) BUG;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (Long i=0; i<individuals; i++) {
    Long availableBits = 0;
    LONG_HANGING_TYPE hanging = 0,
      m = 0;
    Long hangingBits = 0;
    Uchar *a0 = ((Uchar*) ans) + i;					
    LONG_HANGING_TYPE *pM = (LONG_HANGING_TYPE*) (M + i * ldM);
    for(Long j=0; j<=end_snpM1_CpB; j++) {      
      Uchar *a = a0 + ldaByte * j;					
      while (availableBits < ORIGBITSperFIVE) {
	if (hangingBits == 0) {
	  hanging = *(pM++);
	  hangingBits = HTT_BITS;
	}
	m |= hanging << availableBits;
	Long movedBits = hangingBits <= HTT_BITS-availableBits
	  ? hangingBits : HTT_BITS - availableBits;
	//	printf("i=%ld j=%ld %d %ld %d  %d<%d-%d\n", i, j, hangingBits, (Long) M - (Long) pM, movedBits, hangingBits, sizeof(LONG_HANGING_TYPE), availableBits);
	availableBits += movedBits;
	hangingBits -= movedBits;
	hanging >>= movedBits;
      }
      *a = table[m & ORIGBITS_TRANS_MASK];
      m >>= ORIGBITSperFIVE;
      availableBits -= ORIGBITSperFIVE;
    }
  }	
}
 
 
