/*
 Authors 
 Guido Moerkotto
 maintained and modified by Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Guido Moerkotte, Martin Schlather

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


#ifndef miraculix_1bitIntern_H
#define miraculix_1bitIntern_H 1

#include "bitBase.h"

static Ulong scalar(BlockType0 *x, BlockType0 *y, Long blocks,
		      Long DeltaPart2) {
  // Moerkotte BarXXXa, BarXXXb
 // const bool      lTrace = false;
  // lookup table for popcnt of half-byte
  
  // mask for lower 4 bits
  const BlockType lZero = ZERO();

  const BlockType* p1i = (BlockType0*) x;
  const BlockType* p1j = (BlockType0*) y;
  const BlockType* p2i = (BlockType0*) x + DeltaPart2;
  const BlockType* p2j = (BlockType0*) y + DeltaPart2;
  
  BlockType lSumOuter1 = lZero;
  BlockType lSumOuter2 = lZero;
  BlockType lSumOuter4 = lZero;
#if ! defined PlainInteger32 && ! defined PlainInteger64
  BlockType lSumInner1 = lZero;
  BlockType lSumInner2 = lZero;
  BlockType lSumInner4 = lZero;
  const BlockType lLookupTable =
    SETREV8(
	    //  0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F
	    /**/0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
#endif
  
  
 for (char k=0; k<=1; k++) {
   Long bundles;
   int minibundle;
    if (k) {
      minibundle = (int) (blocks % loop_1bit);
      bundles = minibundle > 0;
    } else {
      bundles = blocks / loop_1bit;
      minibundle = (int) loop_1bit;
    }

    for(Long lCurr = 0; lCurr < bundles; lCurr ++) { // loop_1bit1: 16 * 4 uint64_t

 #if ! defined PlainInteger32 && ! defined PlainInteger64
      lSumInner1 = lZero;
      lSumInner2 = lZero;
      lSumInner4 = lZero;
 #endif
     
      // Units: minibundle * UnitsPerBlock // OK
      for(int i = 0; i < minibundle; ++i) { // loop_1bit 2 16 uint64_t
	const BlockType a1i   = LOAD(p1i++);
	//printf("lcurr=%d  %d %d\n", lCurr, i, k);
 	const BlockType a1j   = LOAD(p1j++);
	const BlockType a2i   = LOAD(p2i++);
	const BlockType a2j   = LOAD(p2j++);

	//	printf("XXXlcurr=%d  %d  %d\n", lCurr, i, DeltaPart2);
	// return a1i + a1j + a2i + a2j;
	
	const BlockType b1i1j = AND(a1i, a1j); //
	const BlockType b1i2j = AND(a1i, a2j);
	const BlockType b1j2i = AND(a1j, a2i);
	const BlockType b2i2j = AND(a2i, a2j); //
	
	const BlockType abmx  = OR(b1i2j, b1j2i); //
	
#if ! defined PlainInteger32 && ! defined PlainInteger64
	const BlockType l1i1j = AND(b1i1j, O1F1);
	const BlockType l12mx = AND(abmx , O1F1);
	const BlockType l2i2j = AND(b2i2j, O1F1);
	const BlockType h1i1j = AND(SHR16(b1i1j, 4), O1F1);
	const BlockType h12mx = AND(SHR16(abmx , 4), O1F1);
	const BlockType h2i2j = AND(SHR16(b2i2j, 4), O1F1);
	
	const BlockType pl1   = SHUFFLE8(lLookupTable, l1i1j);
	const BlockType ph1   = SHUFFLE8(lLookupTable, h1i1j);
	const BlockType plmx  = SHUFFLE8(lLookupTable, l12mx);
	const BlockType phmx  = SHUFFLE8(lLookupTable, h12mx);
	const BlockType pl4   = SHUFFLE8(lLookupTable, l2i2j);
	const BlockType ph4   = SHUFFLE8(lLookupTable, h2i2j);
	lSumInner1 = ADD8(lSumInner1, ADD8(pl1, ph1));
	lSumInner2 = ADD8(lSumInner2, ADD8(plmx, phmx));
	lSumInner4 = ADD8(lSumInner4, ADD8(pl4, ph4));
      }
      lSumOuter1 = ADD64(lSumOuter1, SAD8(lSumInner1, lZero));
      lSumOuter2 = ADD64(lSumOuter2, SAD8(lSumInner2, lZero));
      lSumOuter4 = ADD64(lSumOuter4, SAD8(lSumInner4, lZero));
#else       
      lSumOuter1 += POPCNT64(b1i1j);
      lSumOuter2 += POPCNT64(abmx);
      lSumOuter4 += POPCNT64(b2i2j);
    }
 #endif
    }
  }

  lSumOuter2   = SHL64(lSumOuter2, 1);
  lSumOuter4   = SHL64(lSumOuter4, 2);


  const BlockType0 Sum = ADD64(ADD64(lSumOuter1, lSumOuter2), lSumOuter4);
  return HORIZ_ADD64(Sum);
}



#if defined AVX512VPOPCNTDQ || defined AVX512VL
#define CSA(H, L, A, B, C) \
  H = _mm512_ternarylogic_epi32(A,B,C, 0xE8); \
  L = _mm512_ternarylogic_epi32(A,B,C, 0x96);
#else
#define CSA(OVL, RES, IN1, IN2, IN3, TMP) \
   const BlockType TMP = XOR(IN1,IN2); \
   OVL = OR( AND(IN1 , IN2), AND(TMP , IN3) );	\
   RES = XOR(TMP , IN3) ;
#endif

// load two sets A and B aligned
#define LOAD8 \
      a1i_a   = LOAD(p1i++); \
      a1j_a   = LOAD(p1j++); \
      a2i_a   = LOAD(p2i++); \
      a2j_a   = LOAD(p2j++); \
      a1i_b   = LOAD(p1i++); \
      a1j_b   = LOAD(p1j++); \
      a2i_b   = LOAD(p2i++); \
      a2j_b   = LOAD(p2j++); \

// prepare two sets A and B
#define PREP \
      in0_a   = AND(a1i_a, a1j_a); \
      b1i2j_a = AND(a1i_a, a2j_a); \
      b1j2i_a = AND(a1j_a, a2i_a); \
      in2_a   = AND(a2i_a, a2j_a); \
      in1_a   = OR(b1i2j_a, b1j2i_a); \
      in0_b   = AND(a1i_b, a1j_b); \
      b1i2j_b = AND(a1i_b, a2j_b); \
      b1j2i_b = AND(a1j_b, a2i_b); \
      in2_b   = AND(a2i_b, a2j_b); \
      in1_b   = OR(b1i2j_b, b1j2i_b); \


#define POPCNT4(Y,I,X) \
  const BlockType x##X##_lo = AND(                  I##X    , O1F1); \
  const BlockType x##X##_hi = AND(SHR16(I##X, 4), O1F1); \
  const BlockType p##X##_lo = SHUFFLE8(PopCountTable, x##X##_lo); \
  const BlockType p##X##_hi = SHUFFLE8(PopCountTable, x##X##_hi); \
  Y##X##_lo = ADD8(Y##X##_lo, p##X##_lo); \
  Y##X##_hi = ADD8(Y##X##_hi, p##X##_hi); 



#define ADD2BASE(X) \
  /**/ X##_0  = XOR(in0_a , in0_b); \
       o1     = AND(in0_a , in0_b); \
       s1     = XOR(in1_a , in1_b); \
       o2     = AND(in1_a , in1_b); \
       s2     = XOR(in2_a , in2_b); \
  /**/ X##_3  = AND(in2_a , in2_b); \
  /**/ X##_1  = OR(o1 , s1);   \
  /**/ X##_2  = OR(o2 , s2);   

#define CSA4X0(A,B,C,T) \
  CSA(A##_0, lSum0, lSum0, B##_0, C##_0, T##_0) \
  CSA(A##_1, lSum1, lSum1, B##_1, C##_1, T##_1) \
  CSA(A##_2, lSum2, lSum2, B##_2, C##_2, T##_2) \
  CSA(A##_3, lSum3, lSum3, B##_3, C##_3, T##_3) 

#define CSA4X1(A,B,C,T) \
  CSA(A##_0, lSum1, lSum1, B##_0, C##_0, T##_0) \
  CSA(A##_1, lSum2, lSum2, B##_1, C##_1, T##_1) \
  CSA(A##_2, lSum3, lSum3, B##_2, C##_2, T##_2) \
  CSA(A##_3, lSum4, lSum4, B##_3, C##_3, T##_3)

#define CSA4X2(A,B,C,T) \
  CSA(A##_0, lSum2, lSum2, B##_0, C##_0, T##_0 ) \
  CSA(A##_1, lSum3, lSum3, B##_1, C##_1, T##_1 ) \
  CSA(A##_2, lSum4, lSum4, B##_2, C##_2, T##_2 ) \
  CSA(A##_3, lSum5, lSum5, B##_3, C##_3, T##_3 )


#define GET_INPUT LOAD8 PREP


#if ! defined PlainInteger32 && ! defined PlainInteger64

static Ulong scalarB(BlockType0 *x, BlockType0 *y, Long blocks,
		      Long DeltaPart2) {
 // const bool      lTrace = false;
  // lookup table for popcnt of half-byte

 
  // mask for lower 4 bits
 
  const BlockType* p1i = (BlockType0*) x;
  const BlockType* p1j = (BlockType0*) y;
  const BlockType* p2i = (BlockType0*) x + DeltaPart2;
  const BlockType* p2j = (BlockType0*) y + DeltaPart2;

  BlockType a1i_a;   // read from A_1[i,*] representing 2^0
  BlockType a1j_a;   // read from A_1[j,*] representing 2^0
  BlockType a2i_a;   // read from A_2[i,*] representing 2^1
  BlockType a2j_a;   // read from A_2[j,*] representing 2^1
  BlockType b1i2j_a; // intermediate in prep
  BlockType b1j2i_a; // intermediate in prep
  BlockType a1i_b;   // read from A_1[i,*] representing 2^0
  BlockType a1j_b;   // read from A_1[j,*] representing 2^0
  BlockType a2i_b;   // read from A_2[i,*] representing 2^1
  BlockType a2j_b;   // read from A_2[j,*] representing 2^1
  BlockType b1i2j_b; // intermediate in prep
  BlockType b1j2i_b; // intermediate in prep


  BlockType in0_a; // 0. input 0: each bit represents a value of 2^0
  BlockType in0_b; // 1. input 0: each bit represents a value of 2^0
  BlockType in1_a; // 0. input 1: each bit represents a value of 2^1
  BlockType in1_b; // 1. input 1: each bit represents a value of 2^1
  BlockType in2_a; // 0. input 2: each bit represents a value of 2^2
  BlockType in2_b; // 1. input 2: each bit represents a value of 2^2

  // s_i_j: for input 2^i s_i_j contains the j-th bit of the sum, i.e., representing a value of 2^i * 2^j
  BlockType s0_0 = ZERO();  //  1
  BlockType s0_1 = ZERO();  //  2
  BlockType s0_2 = ZERO();  //  3
  BlockType s0_3 = ZERO();  //  4
  BlockType s0_4 = ZERO();  //  5
  BlockType s1_0 = ZERO();  //  6
  BlockType s1_1 = ZERO();  //  7
  BlockType s1_2 = ZERO();  //  8
  BlockType s1_3 = ZERO();  //  9
  BlockType s1_4 = ZERO();  // 10
  BlockType s2_0 = ZERO();  // 11
  BlockType s2_1 = ZERO();  // 12
  BlockType s2_2 = ZERO();  // 13
  BlockType s2_3 = ZERO();  // 14
  BlockType s2_4 = ZERO();  // 15

  // registers for intermediate results
  BlockType aa0_1 = ZERO();  // 
  BlockType aa0_2 = ZERO();  // 
  BlockType aa0_3 = ZERO();  // 
  BlockType aa1_1 = ZERO();  // 
  BlockType aa1_2 = ZERO();  // 
  BlockType aa1_3 = ZERO();  // 
  BlockType aa2_1 = ZERO();  // 
  BlockType aa2_2 = ZERO();  // 
  BlockType aa2_3 = ZERO();  // 
  BlockType bb0_1 = ZERO();  // 
  BlockType bb0_2 = ZERO();  // 
  BlockType bb0_3 = ZERO();  // 
  BlockType bb1_1 = ZERO();  // 
  BlockType bb1_2 = ZERO();  // 
  BlockType bb1_3 = ZERO();  // 
  BlockType bb2_1 = ZERO();  // 
  BlockType bb2_2 = ZERO();  // 
  BlockType bb2_3 = ZERO();  // 

  BlockType lSum0_4 = ZERO();
  BlockType lSum1_4 = ZERO();
  BlockType lSum2_4 = ZERO();

  Long bundles = blocks / loop_1bit;
  #define GET_INPUT LOAD8 PREP
  for(Long lCurr = 0; lCurr < bundles; lCurr ++) { // 
       // 0-3
       GET_INPUT
       CSA(aa0_1, s0_0, s0_0, in0_a, in0_b, tmp0);  // Wertigkeit 2^0
       CSA(aa1_1, s1_0, s1_0, in1_a, in1_b, tmp1);  // Wertigkeit 2^1
       CSA(aa2_1, s2_0, s2_0, in2_a, in2_b, tmp2);  // Wertigkeit 2^2

       GET_INPUT
       CSA(bb0_1, s0_0, s0_0, in0_a, in0_b, tmp3);  // Wertigkeit 2^0
       CSA(bb1_1, s1_0, s1_0, in1_a, in1_b, tmp4);  // Wertigkeit 2^2
       CSA(bb2_1, s2_0, s2_0, in2_a, in2_b, tmp5);  // Wertigkeit 2^2

       CSA(aa0_2, s0_1, s0_1, aa0_1, bb0_1, tmp6);  // Wertigkeit 2^0
       CSA(aa1_2, s1_1, s1_1, aa1_1, bb1_1, tmp7);  // Wertigkeit 2^2
       CSA(aa2_2, s2_1, s2_1, aa2_1, bb2_1, tmp8);  // Wertigkeit 2^2

       // 4-7;
       GET_INPUT
       CSA(aa0_1, s0_0, s0_0, in0_a, in0_b, tmp9);  // Wertigkeit 2^0
       CSA(aa1_1, s1_0, s1_0, in1_a, in1_b, tmp10);  // Wertigkeit 2^1
       CSA(aa2_1, s2_0, s2_0, in2_a, in2_b, tmp11);  // Wertigkeit 2^2

       GET_INPUT
       CSA(bb0_1, s0_0, s0_0, in0_a, in0_b, tmp12);  // Wertigkeit 2^0
       CSA(bb1_1, s1_0, s1_0, in1_a, in1_b, tmp13);  // Wertigkeit 2^1
       CSA(bb2_1, s2_0, s2_0, in2_a, in2_b, tmp14);  // Wertigkeit 2^2

       CSA(bb0_2, s0_1, s0_1, aa0_1, bb0_1, tmp15);  // Wertigkeit 2^0
       CSA(bb1_2, s1_1, s1_1, aa1_1, bb1_1, tmp16);  // Wertigkeit 2^1
       CSA(bb2_2, s2_1, s2_1, aa2_1, bb2_1, tmp17);  // Wertigkeit 2^2

       CSA(aa0_3, s0_2, s0_2, aa0_2, bb0_2, tmp18);  // finally    2^0
       CSA(aa1_3, s1_2, s1_2, aa1_2, bb1_2, tmp19);  // finally    2^1
       CSA(aa2_3, s2_2, s2_2, aa2_2, bb2_2, tmp20);  // finally    2^2

       // 8-11
       GET_INPUT
       CSA(aa0_1, s0_0, s0_0, in0_a, in0_b, tmp21);  // Wertigkeit 2^0
       CSA(aa1_1, s1_0, s1_0, in1_a, in1_b, tmp22);  // Wertigkeit 2^1
       CSA(aa2_1, s2_0, s2_0, in2_a, in2_b, tmp23);  // Wertigkeit 2^2

       GET_INPUT
       CSA(bb0_1, s0_0, s0_0, in0_a, in0_b, tmp24);  // Wertigkeit 2^0
       CSA(bb1_1, s1_0, s1_0, in1_a, in1_b, tmp25);  // Wertigkeit 2^1
       CSA(bb2_1, s2_0, s2_0, in2_a, in2_b, tmp26);  // Wertigkeit 2^2

       CSA(aa0_2, s0_1, s0_1, aa0_1, bb0_1, tmp27);  // Wertigkeit 2^0
       CSA(aa1_2, s1_1, s1_1, aa1_1, bb1_1, tmp28);  // Wertigkeit 2^1
       CSA(aa2_2, s2_1, s2_1, aa2_1, bb2_1, tmp29);  // Wertigkeit 2^2

       // 12-15
       GET_INPUT
       CSA(aa0_1, s0_0, s0_0, in0_a, in0_b, tmp30);  // Wertigkeit 2^0
       CSA(aa1_1, s1_0, s1_0, in1_a, in1_b, tmp31);  // Wertigkeit 2^1
       CSA(aa2_1, s2_0, s2_0, in2_a, in2_b, tmp32);  // Wertigkeit 2^2

       GET_INPUT
       CSA(bb0_1, s0_0, s0_0, in0_a, in0_b, tmp33);  // Wertigkeit 2^0
       CSA(bb1_1, s1_0, s1_0, in1_a, in1_b, tmp34);  // Wertigkeit 2^1
       CSA(bb2_1, s2_0, s2_0, in2_a, in2_b, tmp35);  // Wertigkeit 2^2

       CSA(bb0_2, s0_1, s0_1, aa0_1, bb0_1, tmp36);  // Wertigkeit 2^0
       CSA(bb1_2, s1_1, s1_1, aa1_1, bb1_1, tmp37);  // Wertigkeit 2^1
       CSA(bb2_2, s2_1, s2_1, aa2_1, bb2_1, tmp38);  // Wertigkeit 2^2

       CSA(bb0_3, s0_2, s0_2, aa0_2, bb0_2, tmp39);  // finally    2^0
       CSA(bb1_3, s1_2, s1_2, aa1_2, bb1_2, tmp40);  // finally    2^1
       CSA(bb2_3, s2_2, s2_2, aa2_2, bb2_2, tmp41);  // finally    2^2

       CSA( s0_4, s0_3, s0_3, aa0_3, bb0_3, tmp42);  // 16ner fuer 2^0
       CSA( s1_4, s1_3, s1_3, aa1_3, bb1_3, tmp43);  // 16ner fuer 2^1
       CSA( s2_4, s2_3, s2_3, aa2_3, bb2_3, tmp44);  // 16ner fuer 2^2

       // lSum0_4 += number_of_bits_set<uint>(s0_4); // popcount for 16 (=2^4) * 2^0 // 8-bit cnt
       // lSum1_4 += number_of_bits_set<uint>(s1_4); // popcount for 16 (=2^4) * 2^1 // 8-bit cnt
       // lSum2_4 += number_of_bits_set<uint>(s2_4); // popcount for 16 (=2^4) * 2^2 // 8-bit cnt
       lSum0_4 =ADD64(lSum0_4, POPCNT64(s0_4));
       lSum1_4 =ADD64(lSum1_4, POPCNT64(s1_4));
       lSum2_4 =ADD64(lSum2_4, POPCNT64(s2_4));
  }

  // add up final stuff
  int minibundle = (int) (blocks % loop_1bit);
  // delta = bundles * loop_1bit * UnitsPerBlock; // OK
  BlockType0 lSum;
  lSum = ZERO();
  //  printf("mini %d %d %d %d\n",  minibundle , blocks, loop_1bit, BitsPerBlock);
  for(Long lCurr = 0; lCurr < minibundle; lCurr ++) {
    a1i_a   = LOAD(p1i++);  
    a1j_a   = LOAD(p1j++); 
    a2i_a   = LOAD(p2i++); 
    a2j_a   = LOAD(p2j++);           
    b1i2j_a = AND(a1i_a, a2j_a);
    b1j2i_a = AND(a1j_a, a2i_a);
    in0_a   = AND(a1i_a, a1j_a);
    in2_a   = AND(a2i_a, a2j_a);
    in1_a   = OR(b1i2j_a, b1j2i_a);
    lSum   = ADD64(lSum, POPCNT64(in0_a));
    lSum   = ADD64(lSum, SHL64(POPCNT64(in1_a), 1));
    lSum   = ADD64(lSum, SHL64(POPCNT64(in2_a), 2));
  }

  
  lSum =ADD64(lSum, POPCNT64(s0_0));
  lSum =ADD64(lSum, SHL64(POPCNT64(s0_1), 1));
  lSum =ADD64(lSum, SHL64(POPCNT64(s0_2), 2));
  lSum =ADD64(lSum, SHL64(POPCNT64(s0_3), 3));
  lSum =ADD64(lSum, SHL64(POPCNT64(s1_0), 1));
  lSum =ADD64(lSum, SHL64(POPCNT64(s1_1), 2));
  lSum =ADD64(lSum, SHL64(POPCNT64(s1_2), 3));
  lSum =ADD64(lSum, SHL64(POPCNT64(s1_3), 4));
  lSum =ADD64(lSum, SHL64(POPCNT64(s2_0), 2));
  lSum =ADD64(lSum, SHL64(POPCNT64(s2_1), 3));
  lSum =ADD64(lSum, SHL64(POPCNT64(s2_2), 4));
  lSum =ADD64(lSum, SHL64(POPCNT64(s2_3), 5));
  lSum =ADD64(lSum, SHL64(lSum0_4, 4));
  lSum =ADD64(lSum, SHL64(lSum1_4, 5));
  lSum =ADD64(lSum, SHL64(lSum2_4, 6));

  return HORIZ_ADD64(lSum);
}






static Ulong scalarA(BlockType0 *x, BlockType0 *y, Long blocks,
		      Long DeltaPart2) {
  // Moerkotte Bar256f, Bar256g
  
 // const bool      lTrace = false;
  // lookup table for popcnt of half-byte

  // mask for lower 4 bits
  //  const BlockType lZero = ZERO();

  const BlockType* p1i = (BlockType0*) x;
  const BlockType* p1j = (BlockType0*) y;
  const BlockType* p2i = (BlockType0*) x + DeltaPart2;
  const BlockType* p2j = (BlockType0*) y + DeltaPart2;

  // two sets of variables holding the input and preparing it
  BlockType a1i_a;   // read from A_1[i,*] representing 2^0
  BlockType a1j_a;   // read from A_1[j,*] representing 2^0
  BlockType a2i_a;   // read from A_2[i,*] representing 2^1
  BlockType a2j_a;   // read from A_2[j,*] representing 2^1
  BlockType b1i2j_a; // intermediate in prep
  BlockType b1j2i_a; // intermediate in prep
  BlockType a1i_b;   // read from A_1[i,*] representing 2^0
  BlockType a1j_b;   // read from A_1[j,*] representing 2^0
  BlockType a2i_b;   // read from A_2[i,*] representing 2^1
  BlockType a2j_b;   // read from A_2[j,*] representing 2^1
  BlockType b1i2j_b; // intermediate in prep
  BlockType b1j2i_b; // intermediate in prep

  // add2base temp vars
  BlockType o1, s1, o2, s2;


  BlockType in0_a; // 0. input 0: each bit represents a value of 2^0
  BlockType in0_b; // 1. input 0: each bit represents a value of 2^0
  BlockType in1_a; // 0. input 1: each bit represents a value of 2^1
  BlockType in1_b; // 1. input 1: each bit represents a value of 2^1
  BlockType in2_a; // 0. input 2: each bit represents a value of 2^2
  BlockType in2_b; // 1. input 2: each bit represents a value of 2^2

  // intermediate results to sum up
  BlockType c_0, c_1, c_2, c_3;
  BlockType d_0, d_1, d_2, d_3;
  BlockType e_0, e_1, e_2, e_3;
  BlockType f_0, f_1, f_2, f_3;
  BlockType g_0, g_1, g_2, g_3;
  BlockType h_0, h_1, h_2, h_3;
  BlockType i_0, i_1, i_2, i_3;

  // variables for sums
  BlockType lSum0 = ZERO();
  BlockType lSum1 = ZERO();
  BlockType lSum2 = ZERO();
  BlockType lSum3 = ZERO();
  BlockType lSum4 = ZERO();
  BlockType lSum5 = ZERO();

  // popcounts for bitvectors with a value of 2^4 per bit
  BlockType lSumInner_0_lo = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_0_hi = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_1_lo = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_1_hi = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_2_lo = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_2_hi = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_3_lo = ZERO(); // 32 x 8 bit counters
  BlockType lSumInner_3_hi = ZERO(); // 32 x 8 bit counters

  BlockType lSumOuter_0    = ZERO(); // 64-bit counter summing up lSumInner_0_*
  BlockType lSumOuter_1    = ZERO(); // 64-bit counter summing up lSumInner_1_*
  BlockType lSumOuter_2    = ZERO(); // 64-bit counter summing up lSumInner_2_*
  BlockType lSumOuter_3    = ZERO(); // 64-bit counter summing up lSumInner_3_*

  Long lCurr  = 0; // index of current pointers
  int lCount = 0; // number of executions of inner loop
  Long bundles = blocks / loop_1bit;
  for(; lCurr < bundles; lCurr+=loop_1bit) { // loop_1bit1: 16 * 4 Uintd64_t

       LOAD8 PREP ADD2BASE(c)
       LOAD8 PREP ADD2BASE(d)
       CSA4X0(e,c,d,t0)

       LOAD8 PREP ADD2BASE(c)
       LOAD8 PREP ADD2BASE(d)
       CSA4X0(f,c,d,t1)
       CSA4X1(g,e,f,t2)

       LOAD8 PREP ADD2BASE(c)
       LOAD8 PREP ADD2BASE(d)
       CSA4X0(e,c,d,t3)

       LOAD8 PREP ADD2BASE(c)
       LOAD8 PREP ADD2BASE(d)
       CSA4X0(f,c,d,t4)
       CSA4X1(h,e,f,t5)
       CSA4X2(i,g,h,t6)

       POPCNT4(lSumInner, i, _0);////////
       POPCNT4(lSumInner, i, _1);
       POPCNT4(lSumInner, i, _2);
       POPCNT4(lSumInner, i, _3);

     ++lCount;
     if(60 < lCount) {
       //printf("lcount =%d \n", lCount);
        lSumOuter_0 = ADD64(lSumOuter_0,
			    ADD64(SAD8(lSumInner_0_lo, ZERO()),
				  SAD8(lSumInner_0_hi, ZERO())));
        lSumOuter_1 = ADD64(lSumOuter_1,
			    ADD64(SAD8(lSumInner_1_lo, ZERO()),
				  SAD8(lSumInner_1_hi, ZERO())));
        lSumOuter_2 = ADD64(lSumOuter_2,
			    ADD64(SAD8(lSumInner_2_lo, ZERO()),
				  SAD8(lSumInner_2_hi, ZERO())));
        lSumOuter_3 = ADD64(lSumOuter_3,
			    ADD64(SAD8(lSumInner_3_lo, ZERO()),
				  SAD8(lSumInner_3_hi, ZERO())));
        lSumInner_0_lo = ZERO();
        lSumInner_0_hi = ZERO();
        lSumInner_1_lo = ZERO();
        lSumInner_1_hi = ZERO();
        lSumInner_2_lo = ZERO();
        lSumInner_2_hi = ZERO();
        lSumInner_3_lo = ZERO();
        lSumInner_3_hi = ZERO();
        lCount = 0;
     }     
  }

  // printf("end lcount =%d \n", lCount);
    if(0 != lCount) {
    lSumOuter_0 = ADD64(lSumOuter_0,
                        ADD64(SAD8(lSumInner_0_lo, ZERO()),
			      SAD8(lSumInner_0_hi, ZERO())));
    lSumOuter_1 = ADD64(lSumOuter_1,
                        ADD64(SAD8(lSumInner_1_lo, ZERO()),
			      SAD8(lSumInner_1_hi, ZERO())));
    lSumOuter_2 = ADD64(lSumOuter_2,
                        ADD64(SAD8(lSumInner_2_lo, ZERO()),
			      SAD8(lSumInner_2_hi, ZERO())));
    lSumOuter_3 = ADD64(lSumOuter_3,
                        ADD64(SAD8(lSumInner_3_lo, ZERO()),
			      SAD8(lSumInner_3_hi, ZERO())));
    }

    BlockType lSum256 = ZERO();

  // remaining _mm256i not enough for 16 x 256 bits
  for(; lCurr < blocks; lCurr ++) {
    a1i_a   = LOAD(p1i++);  
    a1j_a   = LOAD(p1j++);  
    a2i_a   = LOAD(p2i++);
    a2j_a   = LOAD(p2j++);
    b1i2j_a = AND(a1i_a, a2j_a);
    b1j2i_a = AND(a1j_a, a2i_a);
    in0_a   = AND(a1i_a, a1j_a);
    in2_a   = AND(a2i_a, a2j_a);
    in1_a   = OR(b1i2j_a, b1j2i_a);
    lSum256 = ADD64(lSum256,      POPCNT64(in0_a));
    lSum256 = ADD64(lSum256,SHL64(POPCNT64(in1_a), 1));
    lSum256 = ADD64(lSum256,SHL64(POPCNT64(in2_a), 2));
  }

  BlockType0 lSum;
  lSum = ADD64(lSum256,    POPCNT64(lSum0));
  lSum = ADD64(lSum, SHL64(POPCNT64(lSum1), 1));
  lSum = ADD64(lSum, SHL64(POPCNT64(lSum2), 2));
  lSum = ADD64(lSum, SHL64(POPCNT64(lSum3), 3));
  lSum = ADD64(lSum, SHL64(POPCNT64(lSum4), 4));
  lSum = ADD64(lSum, SHL64(POPCNT64(lSum5), 5));
  lSum = ADD64(lSum, SHL64(lSumOuter_0, 3));
  lSum = ADD64(lSum, SHL64(lSumOuter_1, 4));
  lSum = ADD64(lSum, SHL64(lSumOuter_2, 5));
  lSum = ADD64(lSum, SHL64(lSumOuter_3, 6));

  return HORIZ_ADD64(lSum);
  
}

#else // 32 or 64
#endif




static void OneBithaplo2genoIntern(unit_t *X, Long snps, Long individuals,
				   Long lda,
				   int VARIABLE_IS_NOT_USED cores,
				   unit_t *Ans, Long ldAns) {
  // genotypes 0 coded as 0,0
  //           1          1,0
  //           2          0,1
  // Achtung!! Nach Aufruf muss sumGeno() berechnet werden!!
  Long blocks = Blocks(snps);
  unit_t *Y = X + lda * individuals,
    *B = Ans + ldAns * individuals;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static)
#endif	  
  for (Long j=0; j<individuals; j++) {
    BlockType0 *xx = (BlockType0*) (X  + j * lda),
      *yy = (BlockType0*) (Y  + j * lda),
      *aa = (BlockType0 *) (Ans  + j * ldAns),
      *bb = (BlockType0 *) (B  + j * ldAns);
    for (Long b=0; b<blocks; b++) {
      const BlockType
	x = LOAD(xx + b),
	y = LOAD(yy + b);
      const BlockType
	first = XOR(x, y),
	second = AND(x, y);
      STORE(aa + b, first);
      STORE(bb + b, second);
    }
  }
}




#define SCALAR012(NR)							\
  Long scalar_1v##NR(unit_t *x, unit_t *y, Long blocks, Long Delta) {	\
    return scalar((BlockType0*) x, (BlockType0*) y, blocks, Delta);	\
  }									\
  void multi_1v##NR(unit_t *x, unit_t *y, Long individuals, Long blocks, Long ldc, \
		   double *ans) {			\
    *ans += (double) scalar((BlockType0*) x, (BlockType0*) y, blocks,	\
			    ldc * individuals / UnitsPerBlock);	/*// OK*/ \
 }									\
  Long scalar_1v##NR##A(unit_t *x, unit_t *y, Long blocks, Long Delta) {	\
    return scalarA((BlockType0*) x, (BlockType0*) y, blocks, Delta);	\
  }									\
  void multi_1v##NR##A(unit_t *x, unit_t *y, Long individuals, Long blocks,\
		       Long ldc, double  *ans) {			\
    *ans += (double) scalarA((BlockType0*) x, (BlockType0*) y, blocks,	\
			     ldc * individuals / UnitsPerBlock); /*// OK*/ \
  }									\
  Long scalar_1v##NR##B(unit_t *x, unit_t *y, Long blocks, Long Delta) {	\
    return scalarB((BlockType0*) x, (BlockType0*) y, blocks, Delta);	\
  }									\
  void multi_1v##NR##B(unit_t *x, unit_t *y, Long individuals,		\
		       Long blocks, Long ldc, double *ans) {		\
    *ans += (double) scalarB((BlockType0*) x, (BlockType0*) y, blocks,	\
		   ldc * individuals / UnitsPerBlock);	/*// OK*/	\
  }									\
  void OneBithaplo2geno##NR(unit_t *X, Long snps, Long individuals,	\
			    Long ldc, int cores, unit_t *ans, Long ldAns) { \
    /* Achtung!! Nach Aufruf muss sumGeno() berechnet werden!! */	\
    OneBithaplo2genoIntern(X, snps, individuals, ldc, cores,\
			   ans, ldAns);					\
  }									\


#endif // 2bitIntern.h
