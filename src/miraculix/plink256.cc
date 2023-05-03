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

#define BitsPerCode 2
#define MY_CODING Plink
#define MY_LDABITALIGN MY_LDABITALIGN_PLINK

#define MY_VARIANT 256

#define CORE_FACTOR 2

#include "Basic_miraculix.h"
#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "Files.h"
#include "Template.h"
#include "plink.h"
#include "MX.h"

#if defined AVX2

ASSERT_SIMD(plink256, avx2);

#include "haplogeno.h"
#include "intrinsicsIntern.h"

/*

IMPORTANT FUNCTIONS IN 2bit256.cc

 */




#define maxVPatonce 8
#define BaseM 4
// Vergleichszahlen BaseM
// 1: 8.0
// 2: 3.9
// 3: 3.3
// 4: 2.6
// 5: 2.5
// 8: 3.8 // somit code nicht auf richtigkeit ueberprueft!


vectorGeno_header(double, double, double, PlinkMatrix256) {
  // calculates  v * t(geno)  for some geno matrix given by "code"
  // In :    V : snps x repetV
  //      code : snps x indiv
  // Out:  Ans : indiv x repetV
  // most inner matrix: (m "doubles of V") x (indivAtOnce1 indiv)
  const int cores = GreaterZero(opt->cores);
  const Long individuals = cols;
  const Long snps = rows;

  
 const Long VDoublesBlocks = DIV_GEQ(repetV, 4);  
  Long Vblocks[maxVPatonce + 1] = {0};
  switch(VDoublesBlocks) {
  case 0 : BUG;
  case 1 : case 2: case 3: case 4: Vblocks[VDoublesBlocks] = 1; break;
  case 5 : Vblocks[2] = Vblocks[3] = 1; break;
  default :
#if BaseM == 3 || BaseM == 4
#define notBaseM (3 + 4 - BaseM)
    Vblocks[notBaseM] = VDoublesBlocks % BaseM;
#if BaseM == 4
    if (Vblocks[notBaseM] > 0) Vblocks[notBaseM] = BaseM - Vblocks[notBaseM];
#endif
    Vblocks[BaseM] = (VDoublesBlocks - notBaseM * Vblocks[notBaseM]) / BaseM;
#else // not  BaseM == 3 || BaseM == 4
    Vblocks[BaseM] = VDoublesBlocks  / BaseM;
#endif
  }

  //printf("VDoublesBlocks = %ld\n", VDoublesBlocks);
  if (false) for (int m=0; m<=maxVPatonce; m++) PRINTF("%d: %ld\n", m, Vblocks[m]);

  //  if (Vblocks[0] > 0 || Vblocks[1] > 0 || Vblocks[2] > 0 || Vblocks[3] > 0) BUG;
  
  if (Vblocks[0] < 0 || Vblocks[1] < 0 || Vblocks[2] < 0 ||
      Vblocks[3] < 0 || Vblocks[4] < 0 ||
      1 * Vblocks[1] + 2 * Vblocks[2] + 3 * Vblocks[3] + 4 * Vblocks[4] +
      5 * Vblocks[5] + 6 * Vblocks[6] + 7 * Vblocks[7] + 8 * Vblocks[8]
      !=
      VDoublesBlocks) BUG;

  //make standalone &&  Wageningen/run_gcc snps=1280 ind=1369 cores=10 coding=4  mode=3 SNPmode=3 repetV=16 trials=100 centered=0 variant=257 cmp=7

 
#define lenSNPtype1 8
#define lenSNPtype2 8
#define lenSNPtype3 8
#define lenSNPtype4 8
#define lenSNPtype5 8
#define lenSNPtype8 8
#define snpType2 Ulong
#define snpType1 Ulong
#define snpType3 Ulong
#define snpType4 Ulong  // int schlechter als Long
#define snpType5 Ulong
#define snpType8 Ulong

#define indivAtOnce1 8 // <= 15
#define indivAtOnce2 7 // <= 7
#define indivAtOnce3 4 // <= 4
#define indivAtOnce4 3 // <= 3 
#define indivAtOnce5 2 // <= 3
#define indivAtOnce8 1 // <= 3
  const int indivAtO[maxVPatonce + 1] =
    {0, indivAtOnce1, indivAtOnce2, indivAtOnce3, indivAtOnce4, indivAtOnce5,
     0, 0, indivAtOnce8};
  Long indivAtO_blocks[maxVPatonce + 1] = { 0 };
  for (int m = 1; m <= maxVPatonce; m++) {
    if ( (indivAtO[m] + 1) * m > 16) BUG;
    if (indivAtO[m] == 0) continue;
    indivAtO_blocks[m] = DIV_GEQ(individuals, indivAtO[m]);
    //   printf("m=%d, indiv %ld / %d = %ld\n", m, individuals, indivAtO[m],indivAtO_blocks[m]);
  }


#define miniSNP1 (lenSNPtype1 * CodesPerByte)
#define miniSNP2 (lenSNPtype2 * CodesPerByte)
#define miniSNP3 (lenSNPtype3 * CodesPerByte)
#define miniSNP4 (lenSNPtype4 * CodesPerByte)
#define miniSNP5 (lenSNPtype5 * CodesPerByte)
#define miniSNP8 (lenSNPtype8 * CodesPerByte)
  const int miniSNP[maxVPatonce + 1] = // guaranteed alignment of TV
    {0, miniSNP1, miniSNP2, miniSNP3, miniSNP4, miniSNP5, 0, 0, miniSNP8};
  
  Long miniSNP_blocks[maxVPatonce + 1] = { 0 };
  Long DoublesTV[maxVPatonce + 1] = { 0 };
  Long DoublesAns[maxVPatonce + 1] = { 0 };
  Long bytesTV = 0,
    bytesTans =0;
  for (int m = 1; m <= maxVPatonce; m++) {
    if (indivAtO[m] == 0) continue;
    miniSNP_blocks[m] = DIV_GEQ(snps, miniSNP[m]);
    DoublesTV[m] = m * miniSNP[m] * miniSNP_blocks[m];
    DoublesAns[m] = m * indivAtO[m] * indivAtO_blocks[m];
    //   printf("m=%d %ld %ld\n %d %ld %d %ld\n\n", m, DoublesTV[m], DoublesAns[m], miniSNP[m],miniSNP_blocks[m],indivAtO[m],indivAtO_blocks[m]);
    
    bytesTV = MAX(DoublesTV[m], bytesTV);
    bytesTans = MAX(DoublesAns[m], bytesTans);
  }
  bytesTV *= sizeof(Doubles);
  bytesTans *= sizeof(Doubles);
  //  printf("bytesTV = %ld %ld\n", bytesTV, bytesTans);
  Doubles *TV = (Doubles*) MALLOC(bytesTV);
  Doubles *Tans = (Doubles*) MALLOC(bytesTans);
  assert(TV != NULL);
  assert(Tans != NULL);
 

#define miniIndiv1 8
#define miniIndiv2 8
#define miniIndiv3 8
#define miniIndiv4 8 // besser als 16 & 2
#define miniIndiv5 8
#define miniIndiv8 8


#define littleSNP1 4
#define littleSNP2 4
#define littleSNP3 4
#define littleSNP4 4 // 4 besser als 16 -- 1 moeglicerweise optimum
#define littleSNP5 4
#define littleSNP8 4


#define littleIndiv1 (4 * miniIndiv1)  
#define littleIndiv2 (4 * miniIndiv2)  
#define littleIndiv3 (4 * miniIndiv3)    
#define littleIndiv4 (4 * miniIndiv4) // 1 & 64 schlechter
#define littleIndiv5 (4 * miniIndiv5)  
#define littleIndiv8 (4 * miniIndiv8)  


#define BpC BitsPerCode
	
#define V_LOAD(VM,DUMMY) const Doubles V##VM = LOADuDOUBLE((double*) (pTV + VM)); if (false ) PRINTF("get %ld %f %f %f %f\n", (double*) (pTV + VM) - (double*) TV, ((double*) &V##VM)[0],((double*) &V##VM)[1],((double*)& V##VM)[2], ((double*)& V##VM)[3]);
#define LOAD_V MULTI_V(V_LOAD,0,;)

#define JJ_CONDLOAD(IN) const SNP_TYPE c##IN = (SNP_TYPE*) (C + (IN) * lda) >= (SNP_TYPE*) codeEnd ? 0 : *((SNP_TYPE*) (C + (IN) * lda))
#define JJ_LOAD(IN) const SNP_TYPE c##IN = *((SNP_TYPE*) (C + (IN) * lda)); if (false) PRINTF("&c%s=%ld; %lu\n", #IN,(Uchar*) (C+(IN) * lda) - (Uchar*) code, (Ulong) c##IN); 
  
#define JJ_SHR(IN) SNP_TYPE d##IN = c##IN >> miniJ
#define SHR_JJ MULTI_INDIV(JJ_SHR,;)

#define JJ_ANDMASK(IN) d##IN &= CodeMask
#define ANDMASK_JJ MULTI_INDIV(JJ_ANDMASK,;)

#define ANS_LOAD(VM,IN) if (false) PRINTF("&Ans[%s,%s] = %lu <= %ld; %ld\n", #VM, #IN, (double*) (tAns + (VM) + (IN) * M) - (double*) Tans, bytesTans / 8, bytesTans);if (false && (double*) (tAns + (VM) + (IN) * M) - (double*) Tans > bytesTans / 8) BUG; Doubles A##VM##_##IN = LOADuDOUBLE((double*) (tAns + (VM) + (IN) * M))
#define LOAD_ANS MULTI_ANS(ANS_LOAD,;)

#define ANS_STORE(VM,IN) STOREuDOUBLE((double*) (tAns + (VM) + (IN) * M), A##VM##_##IN)
#define STORE_ANS MULTI_ANS(ANS_STORE,;)


#define SCALAR_PROD_SINGLE(VM,IN) A##VM##_##IN += V##VM
#define SCALAR_PROD_MINI(IN, M, Z)					\
  if (false)  PRINTF("d%s=%d \n", #IN, (int) d##IN);			\
  if (d##IN <= CODE_LASTZERO) goto zero##Z##M##IN;			\
  MULTI_V(SCALAR_PROD_SINGLE,IN,;);					\
  if (d##IN == CODE_ONE) goto zero##Z##M##IN; 				\
  MULTI_V(SCALAR_PROD_SINGLE,IN,;);					\
 zero##Z##M##IN:							       


  #define MINI_SCALAR_PROD(X) MULTI_M(SCALAR_PROD_MINI,M,X,;)

      

#define LOOP(NR,MM)							\
  for (Long iB=0; iB<IndivAtO_blocks; iB += LITTLE_INDIV) {   /* I indep*/\
    Long iLittleEnd = MIN(iB + LITTLE_INDIV, IndivAtO_blocks);		\
    Doubles *tAns1 = Tans + ANS_MINI * iB;				\
    Doubles *pTV1 = TV;							\
									\
    for (Long jjB=0; jjB<MiniSNP_blocks; jjB += LITTLE_SNP,  /* S finish*/ \
	   pTV1 += TV_MINI * LITTLE_SNP) {				\
      unit_t* C2 = (unit_t*) (((SNP_TYPE*) code) + jjB);		\
      Doubles *tAns0 = tAns1;						\
      Doubles *pTV0end = pTV1 +						\
	TV_MINIINDIV * MIN(LITTLE_SNP, MiniSNP_blocks - jjB);		\
									\
      for (Long iLittle=iB; iLittle < iLittleEnd; iLittle+=MINI_INDIV,	\
	     tAns0 += ANS_MINIINDIV) {                       /* I interm*/ \
	SNP_TYPE *C1 = (SNP_TYPE*) (C2 + LdaAtonce * iLittle);		\
	Doubles *tAnsEnd = tAns0 + MIN(ANS_MINIINDIV, ANS_MINI *	\
				       (int) (IndivAtO_blocks - iLittle));\
									\
	if (iLittle+MINI_INDIV <= iLittleEnd){				\
	  if (false) PRINTF("iB=%ld<%ld jjB=%ld<%ld litt=%ld,%ld<%ld; ", \
		   iB,IndivAtO_blocks, jjB, MiniSNP_blocks,iLittle,	\
		   iLittle+MINI_INDIV, iLittleEnd);			\
	  if (false) PRINTF("tAns: %ld %ld end= %ld / %d = %ld; min=%d<=%d\n", \
		   (double*) tAns1 - (double*) Tans,			\
		   (double*) tAns0 - (double*) Tans,			\
		   (double*) tAnsEnd - (double*) Tans,			\
		   ANS_MINI,						\
		   ((double*) tAnsEnd - (double*) Tans) / ANS_MINI,	\
		   MIN(ANS_MINIINDIV,ANS_MINI*(int)(IndivAtO_blocks-iLittle)),\
		   ANS_MINIINDIV);					\
	  INNERLOOP##NR(,MM);						\
	} else {							\
	  INNERLOOP##NR(COND,MM);						\
	}								\
      }									\
    }									\
  }
  
  

#define INNERLOOP(COND,MM)						\
  for (Doubles *pTV0 = pTV1; pTV0 < pTV0end;        /*S: interm*/	\
       pTV0 += TV_MINIINDIV,						\
	 C1 ++) {							\
    unit_t* C = (unit_t*) C1;						\
    for (Doubles *tAns = tAns0;  tAns < tAnsEnd;   /* I: V recyl*/	\
	 tAns += ANS_MINI, C += LdaAtonce	) {			\
      /* Umsortieren von code so dass C += MINI_SNP o.ae. */		\
      /* bringt ca 5%; lohnt Aufwand nicht */				\
      MULTI_INDIV(JJ_##COND##LOAD,;); /* C[e + (IN) * lda] */		\
      Doubles *pTV = pTV0;						\
      LOAD_ANS; /* tAns + (VM) + (IS) * M) */				\
      /* memory addressed  below is guaranteed to exist */		\
      /* no matter what the user input is */				\
      for (int miniJ=0; miniJ < (int)(BpC * MINI_SNP);  /*S no Ans*/	\
	   miniJ += BpC) { /*unroll bringt nix*/			\
	LOAD_V;/* (pTV + VM) */						\
	pTV += M;							\
	SHR_JJ;/* >> miniJ */						\
	ANDMASK_JJ;							\
	MINI_SCALAR_PROD(COND##MM);/*(repetV x SNPS) +=LOAD_V %*% t(LOAD_JJ)*/ \
      }									\
      STORE_ANS;							\
    }									\
  }




#define INNERLOOP1(COND,MM)						\
  for (Doubles *pTV0 = pTV1; pTV0 < pTV0end;        /*S: interm*/	\
       pTV0 += TV_MINIINDIV,						\
	 C1 ++) {							\
    unit_t* C = (unit_t*) C1;						\
    for (Doubles *tAns = tAns0;  tAns < tAnsEnd;   /* I: V recyl*/	\
	 tAns += ANS_MINI, C += LdaAtonce	) {			\
       MULTI_INDIV(JJ_##COND##LOAD,;); /* C[e + (IN) * lda] */		\
       Doubles *pTV = pTV0;						\
       LOAD_ANS; /* tAns + (VM) + (IS) * M) */				\
       for (int miniJ=0; miniJ < (int)(BpC * MINI_SNP);  /*S no Ans*/	\
	   miniJ += BpC, pTV += M) { /*unroll bringt nix*/		\
	SHR_JJ;/* >> miniJ */						\
	ANDMASK_JJ;							\
	LOAD_V;/* (pTV + VM) */	/*repetV  +=LOAD_V */			\
	if (d0 <= CODE_LASTZERO) continue;				\
	MULTI_V(SCALAR_PROD_SINGLE,0,;);				\
	if (d0 == CODE_ONE) continue;					\
	MULTI_V(SCALAR_PROD_SINGLE,0,;);				\
       }								\
      STORE_ANS;							\
    }									\
  }

  
#define TV_MINI (M * INDIV_ATONCE)
#define TV_MINIINDIV  (TV_MINI * MINI_INDIV)

#define ANS_MINI (M * INDIV_ATONCE)
#define ANS_MINIINDIV (ANS_MINI * MINI_INDIV)
#define ANS_LITTLEINDIV (ANS_MINI * LITTLE_INDIV)

  BUG; //nachfolgend fueyhrt zu grossen abweichungen:
  //  make standalone &&  Wageningen/run_gcc snps=1280 ind=1369 cores=10 coding=4  mode=3 SNPmode=3 repetV=16 trials=100 centered=0 variant=257 cmp=7
  
  MEMSET(Ans, 0, sizeof(double) * ldAns * repetV);
  Long cum_m = 0;
  double *pV = V;
  double *pA = Ans;

  unit_t *codeEnd = code + lda * individuals;
  
  assert(ldAns == individuals);
	
  for (int m = 1; m <= maxVPatonce; m++) {
    if (indivAtO[m] == 0) continue;

    const int M4      =  m * doubles;
    //    const int MiniSNP = miniSNP[m];
    const Long MiniSNP_blocks  = miniSNP_blocks[m];
    const Long IndivAtO_blocks = indivAtO_blocks[m];
    const int  Indiv_atonce    = indivAtO[m];
    const Long LdaAtonce       = lda * Indiv_atonce;
   
    for (int iVB=0; iVB < Vblocks[m]; iVB++,
	   pV += M4 * ldV, pA += M4 * ldAns, cum_m += M4) {
      
       /// set TV
      MEMSET(TV, 0, bytesTV);
      MEMSET(Tans, 0, bytesTans);     
      const int m4 = MIN(M4, (int) (repetV - cum_m));
      //      printf("B m=%d iVB=%d M4=%d m4=%d\n", m, iVB, M4, m4);
      
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) 
#endif 
     for (int r = 0; r < m4; r++) {
	double *rrV = pV + r * ldV;
	double *pTV = ((double*) TV) + r;
	for (int j=0; j<snps; j++) {
	  //	  if (false)  printf("rrV[%d + %ld = %ld]=%f -> %d max=%ld\n", j, r * ldV,j + r * ldV, rrV[j], m4 *j+r, bytesTV/8);
	  pTV[m4* j] = rrV[j];
	}
      }
     
      
      switch(m) {
      case 0: BUG; break;
      case 1: {
#if defined M
#undef M
#undef LEN_SNP_TYPE 
#undef SNP_TYPE     
#undef INDIV_ATONCE 
#undef MINI_SNP
#undef MINI_INDIV
#undef LITTLE_SNP   
#undef LITTLE_INDIV
#undef MULTI_V
#undef MULTI_INDIV
#undef MULTI_M
#undef MULTI_ANS
#endif

#define M 1
#define LEN_SNP_TYPE  lenSNPtype1
#define SNP_TYPE         snpType1
#define INDIV_ATONCE indivAtOnce1
#define MINI_SNP         miniSNP1
#define MINI_INDIV     miniIndiv1
#define LITTLE_SNP     littleSNP1
#define LITTLE_INDIV littleIndiv1
#define MULTI_V(X,N,SZ) X(0,N) 
#define MULTI_INDIV(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) SZ X(4) SZ X(5) SZ X(6) SZ X(7) SZ X(8) 
#define MULTI_M(X,M,Z,SZ) X(0,M,Z) SZ X(1,M,Z) SZ X(2,M,Z) SZ X(3,M,Z) SZ X(4,M,Z) SZ X(5,M,Z) SZ X(6,M,Z) SZ X(7,M,Z) SZ X(8,M,Z)
#define MULTI_ANS(X,SZ) MULTI_V(X,0,SZ) SZ  MULTI_V(X,1,SZ) SZ MULTI_V(X,2,SZ) SZ MULTI_V(X,3,SZ) SZ MULTI_V(X,4,SZ) SZ MULTI_V(X,5,SZ) SZ MULTI_V(X,6,SZ) SZ  MULTI_V(X,7,SZ) SZ MULTI_V(X,8,SZ) 
    
#ifdef DO_PARALLEL
	#pragma omp parallel for num_threads(cores) schedule(static)
#endif
	LOOP(,M);
       
	
      }	break;
      case 2: {
#if defined M
#undef M
#undef LEN_SNP_TYPE 
#undef SNP_TYPE     
#undef INDIV_ATONCE 
#undef MINI_SNP
#undef MINI_INDIV
#undef LITTLE_SNP   
#undef LITTLE_INDIV
#undef MULTI_V
#undef MULTI_INDIV
#undef MULTI_M
#undef MULTI_ANS
#endif
	
#define M 2
#define LEN_SNP_TYPE  lenSNPtype2
#define SNP_TYPE         snpType2
#define INDIV_ATONCE indivAtOnce2
#define MINI_SNP         miniSNP2
#define MINI_INDIV     miniIndiv2
#define LITTLE_SNP     littleSNP2
#define LITTLE_INDIV littleIndiv2
#define MULTI_V(X,N,SZ) X(0,N) SZ X(1,N) 
#define MULTI_INDIV(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3)  SZ X(4) SZ X(5) SZ X(6) 
#define MULTI_M(X,M,Z,SZ) X(0,M,Z) SZ X(1,M,Z) SZ X(2,M,Z) SZ X(3,M,Z) SZ X(4,M,Z) SZ X(5,M,Z) SZ X(6,M,Z)
#define MULTI_ANS(X,SZ) MULTI_V(X,0,SZ) SZ  MULTI_V(X,1,SZ) SZ MULTI_V(X,2,SZ) SZ MULTI_V(X,3,SZ) SZ MULTI_V(X,4,SZ) SZ MULTI_V(X,5,SZ) SZ MULTI_V(X,6,SZ)

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
	LOOP(,M);
 	

      }	break;
      case 3: {
	
#if defined M
#undef M
#undef LEN_SNP_TYPE 
#undef SNP_TYPE     
#undef INDIV_ATONCE 
#undef MINI_SNP     
#undef MINI_INDIV
#undef LITTLE_SNP   
#undef LITTLE_INDIV
#undef MULTI_V
#undef MULTI_INDIV
#undef MULTI_M
#undef MULTI_ANS
#endif

#define M 3
#define LEN_SNP_TYPE  lenSNPtype3
#define SNP_TYPE         snpType3
#define INDIV_ATONCE indivAtOnce3
#define MINI_SNP         miniSNP3
#define MINI_INDIV     miniIndiv3
#define LITTLE_SNP     littleSNP3
#define LITTLE_INDIV littleIndiv3
#define MULTI_V(X,N,SZ) X(0,N) SZ X(1,N) SZ X(2,N) 
#define MULTI_INDIV(X,SZ) X(0) SZ X(1) SZ X(2) SZ X(3) 
#define MULTI_M(X,M,Z,SZ) X(0,M,Z) SZ X(1,M,Z) SZ X(2,M,Z) SZ X(3,M,Z) 
#define MULTI_ANS(X,SZ) MULTI_V(X,0,SZ) SZ  MULTI_V(X,1,SZ) SZ MULTI_V(X,2,SZ) SZ MULTI_V(X,3,SZ)

    
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
	LOOP(,M)
       
      }	break;
      case 4: {

#if defined M
#undef M
#undef LEN_SNP_TYPE
#undef SNP_TYPE     
#undef INDIV_ATONCE 
#undef MINI_SNP     
#undef MINI_INDIV
#undef LITTLE_SNP   
#undef LITTLE_INDIV
#undef MULTI_V
#undef MULTI_INDIV
#undef MULTI_M
#undef MULTI_ANS
#endif
      
#define M 4
#define LEN_SNP_TYPE  lenSNPtype4
#define SNP_TYPE         snpType4
#define INDIV_ATONCE indivAtOnce4
#define MINI_SNP         miniSNP4
#define MINI_INDIV     miniIndiv4
#define LITTLE_SNP     littleSNP4
#define LITTLE_INDIV littleIndiv4
#define MULTI_V(X,N,SZ) X(0,N) SZ X(1,N) SZ X(2,N) SZ X(3,N)
#define MULTI_INDIV(X,SZ) X(0) SZ X(1) SZ X(2) 
#define MULTI_M(X,M,Z,SZ) X(0,M,Z) SZ X(1,M,Z) SZ X(2,M,Z) 
#define MULTI_ANS(X,SZ) MULTI_V(X,0,SZ) SZ  MULTI_V(X,1,SZ) SZ MULTI_V(X,2,SZ)
	

	assert(sizeof(SNP_TYPE) == LEN_SNP_TYPE);
	
	// printf("cores=%d indivAtO_blocks=%ld\n", cores, IndivAtO_blocks);

	// NOTE:
	// * indiv is treated in indiv-miniblocks
	// * snps  is treated in SNP-miniblocks (except most inner loop)
	// * repetV is fixes as a m-block of Doubles

	// loop character:
	// I block (independent Ans blocks)
	// S block (finishes Ans block)
	//   I block (intermediate)
	//   S block (intermediate)	
	//       I block avoids too many V loads; need load of Ans
	//        S block without load of Code nor Ans; needs load of V only
	//                load Ans / load V efficient
	
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
	// schedule(static, n), n < 25 ist wurscht
	//ordered // ist wurscht
#endif
	LOOP(,M)
 	
      } break;
      case 5: {

#if defined M
#undef M
#undef LEN_SNP_TYPE
#undef SNP_TYPE     
#undef INDIV_ATONCE 
#undef MINI_SNP     
#undef MINI_INDIV
#undef LITTLE_SNP   
#undef LITTLE_INDIV
#undef MULTI_V
#undef MULTI_INDIV
#undef MULTI_M
#undef MULTI_ANS
#endif

#define M 5
#define LEN_SNP_TYPE  lenSNPtype5
#define SNP_TYPE         snpType5
#define INDIV_ATONCE indivAtOnce5
#define MINI_SNP         miniSNP5
#define MINI_INDIV     miniIndiv5
#define LITTLE_SNP     littleSNP5
#define LITTLE_INDIV littleIndiv5
#define MULTI_V(X,N,SZ) X(0,N) SZ X(1,N) SZ X(2,N) SZ X(3,N) SZ X(4,N)
#define MULTI_INDIV(X,SZ) X(0) SZ X(1) 
#define MULTI_M(X,M,Z,SZ) X(0,M,Z) SZ X(1,M,Z) 
#define MULTI_ANS(X,SZ) MULTI_V(X,0,SZ) SZ  MULTI_V(X,1,SZ) 
    
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
	LOOP(,M);
 	
      } break;
  

      case 8: {
  
#if defined M
#undef M
#undef LEN_SNP_TYPE
#undef SNP_TYPE     
#undef INDIV_ATONCE 
#undef MINI_SNP     
#undef MINI_INDIV
#undef LITTLE_SNP   
#undef LITTLE_INDIV
#undef MULTI_V
#undef MULTI_INDIV
#undef MULTI_M
#undef MULTI_ANS
#endif
     
#define M 8
#define LEN_SNP_TYPE  lenSNPtype8
#define SNP_TYPE         snpType8
#define INDIV_ATONCE indivAtOnce8
#define MINI_SNP         miniSNP8
#define MINI_INDIV     miniIndiv8
#define LITTLE_SNP     littleSNP8
#define LITTLE_INDIV littleIndiv8
#define MULTI_V(X,N,SZ) X(0,N) SZ X(1,N) SZ X(2,N) SZ X(3,N) SZ X(4,N) SZ X(5,N) SZ X(6,N) SZ X(7,N) 
#define MULTI_INDIV(X,SZ) X(0) 
#define MULTI_M(X,M,Z,SZ) X(0,M,Z) 
#define MULTI_ANS(X,SZ) MULTI_V(X,0,SZ) 

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif 
	LOOP(1,M);

      } break;
  
	
      default: BUG;
      }

      //      printf("Near end\n");

      // transpose back to Ans
      const int M4Indiv = M4 * Indiv_atonce;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif 
      for (int i=0; i<individuals; i++) {
	int iAtonceBlock = i / Indiv_atonce;
	int iAtO = i % Indiv_atonce;
	double *rrAns = pA + i;
	double *pTans = ((double*) Tans) + iAtonceBlock * M4Indiv + iAtO * M4;
	for (int r = 0; r < m4; r++) {
	  if (ldAns * r + i >= repetV * individuals) 
	    //	    printf("%ld * %d + %d >=  %ld * %ld\n", ldAns, r , i , repetV , individuals);
	  //
	  assert(ldAns * r + i < repetV * individuals);
	  rrAns[ldAns * r] = pTans[r];
	}
      } // i

      //printf("end iVB\n");
      
    } // iVB
  } // m
   
  
  FREE(TV);
  FREE(Tans);
 

  // printf("Plink matrix done\n");
}



#else // !defined AVX
#include "avx_miss.h"

vectorGeno_header(double, double, double, 2matrix256) Sv


// DeclareVersion(256,Su,Sv,Sm)
SIMD_MISS(plink256, avx2);
#endif  


// 4 & 5 : 5codes und plink ungefaehr gleichgut. Darunter 5codes besser.

// make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=10 coding=5   mode=3 SNPmode=3 repetV=32 trials=1 centered=0 variant=257
// 5codes : 7.136 ( 54.504)
// plink  : 4.153 ( 37.615)
// 4, SNPmode=0,2 : gleich. Schraeg.

//make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=20 coding=5   mode=3 SNPmode=3 repetV=32 trials=1 centered=0 variant=257
// 5 : 7.709 (109.217)
// 4 : 3.192 ( 57.394) 

//  make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=20 coding=4   mode=3 SNPmode=3 repetV=8 trials=1 centered=0 variant=257
// 5 : 1.892 ( 25.512)
// 4 : 0.739 ( 12.754) sec for gV_vG (plink; G * v; no means)

//  make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=20 coding=4   mode=3 SNPmode=3 repetV=6 trials=1 centered=0 variant=257
// 5 : 1.458 ( 19.341)
// 4 : 0.990 ( 16.696)


//  make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=20 coding=4   mode=3 SNPmode=3 repetV=5 trials=1 centered=0 variant=257
// 5 : 1.177 ( 15.924)
// 4 : 0.996 ( 16.743)


//  make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=20 coding=4   mode=3 SNPmode=3 repetV=4 trials=1 centered=0 variant=257
// 5 : 0.955 ( 12.762) 
// 4 : 0.762 ( 13.330)


//  make standalone && Wageningen/run_gcc snps=128000 ind=136960 cores=20 coding=4   mode=3 SNPmode=3 repetV=2 trials=1 centered=0 variant=257
// 5 : 0.511 (  6.419)
// 4 : 0.942 ( 16.091) 
