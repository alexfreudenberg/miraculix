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



#ifndef miraculix_MXINFO_H
#define miraculix_MXINFO_H 1



typedef int Rint;
typedef  char table_type;
				       

Uint Inti(SEXP X, Uint i);
#define Int0(X) Inti(X, 0U)

	
Long snps2entities(Long snps, Long CodesPerEntity);


Long *GetInfoUnchecked(SEXP Code);
Long *GetInfo(SEXP M);
Long *GetInfo(SEXP Code, int nullOK);
void printINFO(SEXP S);


#define DIV_GEQ(X,R) (1 + ((X) - 1) / (R))
#define DIV_LESS(X,R) (((X) - 1) / (R))
#define DIV_LEQ(X,R) ((X) / (R))

#define ROUND_GEQ(X,R) (DIV_GEQ(X,R) * (R))
#define ROUND_LESS(X,R) (DIV_LESS(X,R) * (R))
#define ROUND_LEQ(X,R) (DIV_LEQ(X,R) * (R))





// ACHTUNG ! uintptr_t ist ein Zeiger auf Byte nicht auf uint!
#if defined BitsPerByte
static inline  uintptr_t algn_generalL(int *X, Long LDAbitalign) {
  Long f = LDAbitalign / BitsPerByte;
  assert(ROUND_GEQ((uintptr_t) X, f) >= (uintptr_t) X);
  return ROUND_GEQ((uintptr_t) X, f);
}
static inline void *algn_generalS(SEXP Code, Long LDAbitalign) {
  assert(getAttribPointer(Code, Information) == R_NilValue);
  uintptr_t X = (uintptr_t) INTEGER(Code);
  Long f = LDAbitalign / BitsPerByte;
  assert(ROUND_GEQ(X, f) >= X);
  return (void*) ROUND_GEQ(X, f);
}
#endif

#define FROMINPUT  ((int) (*(pM++)))
#define FROMHAPLO GetTwoBitHaplo(pM, s)


#define ASSERT_ENDIAN0(Bigendian)		\
  if ((Bigendian) == info[BIGENDIAN]) {} else				\
    ERR2("endians do not match in line %d if %s", __LINE__, __FILE__)

#define ASSERT_STORED_ENDIAN ASSERT_ENDIAN0(opt->bigendian)
 
#define ASSERT_LITTLE_ENDIAN ASSERT_STORED_ENDIAN; if (info[BIGENDIAN]) { BUG; } else { }

#define isRoughlyGeno(X) ((X) <= LastGenoCoding || ((X) >= CorrespondingGeno && (X) <= LastTechnicalGenoCoding))
#define isNarrowGeno(X) (isRoughlyGeno(X) && (X) != AutoCoding)
#define anyByteGeno(X) ((X) == OneByteGeno || (X) == FourByteGeno)
#define isHaplo(X) ((UnknownSNPcoding == 31) /*no contents --just an assert*/ \
		    && (X) >= FirstHaploCoding && (X) <= LastHaploCoding)
#define userHaplo(X) ((X) >= FirstUserHaploCoding &&	\
		      (X) <= LastUserHaploCoding)
#define isCompressedHaplo(X) (isHaplo(X) && !userHaplo(X))
#define transposedGeno(X) ((X) >= FiveCodesTransposed && (X) <= LastTransposedGenoCoding)
#define hasTransposed(X) ((X) == TwoBitGeno || (X) == FiveCodes || (X) == Plink|| (X) == OrigPlink)
#define transposedHaplo(X) ((X) ==  EightByteHaplo_Tr ||	\
			    (X) == FourByteSingleBit_Tr ||	\
			    (X) == FourByte_Tr)
#define transposed(X) (transposedGeno(X) || transposedHaplo(X))

#define is2Bit(X) ((X)==OneBitGeno  || (X)==TwoBitGeno ||	\
		   (X)==TwoBitHaplo || (X)==OneBitHaplo)
#define isOneByte(X) ((X) == OneByteGeno || (X) == OneByteHaplo)
#define isOneByteHaplo(X, OneHaploOnly) ((X) == OneByteHaplo || \
					 ((OneHaploOnly) && (X) == OneByteGeno))
#define isFourByteHaplo(X, OneHaploOnly) \
  ((X) == FourByteHaplo || (X) == FourByteSingleBit ||	\
   ((OneHaploOnly) && (X) == FourByteGeno))
#define isSomeByteHaplo(X, OneHaploOnly) (isOneByteHaplo(X, OneHaploOnly) || \
					  isFourByteHaplo(X, OneHaploOnly))

#define isCompressed(X) ( ((X) > AutoCoding && (X) <= OneByteGeno) ||	\
			  ((X) >= FiveCodesTransposed &&		\
			   (X) <= LastTransposedGenoCoding)        ||	\
			  ((X) >= OneBitHaplo && (X) <= OneByteHaplo) )
#define uprightlyPacked(X,VARIANT)					\
  ((X) == FiveCodes || (X) == FiveCodesTransposed)

#define is32BitAligned(X) ((X) == FourByteGeno || (X) == FourByte_Tr ||	 \
			   ((X)>=FirstUserHaploCoding &&(X)<=EightByteHaplo_Tr))

#define UnCompressedDoubleCols(X, RELAXED)				\
  ((X) == FourByteHaplo || (X) == OneByteHaplo || (X) == FourByteSingleBit || \
   ((RELAXED) && ((X) == FourByteGeno || (X) == OneByteGeno)))

#define UnCompressedDoubleRows(X) ((X) == EightByteHaplo)

#define doubledCols(X, RELAXED) ((X) == OneBitHaplo || (X) == OneBitGeno || \
				 (X) == EightByteHaplo_Tr ||		\
				 UnCompressedDoubleCols(X, RELAXED) )

#define doubledSnps(X) ((X) == EightByteHaplo_Tr || (X) == EightByteHaplo)

#define isOneSet2G(X) ((X) == OneByteGeno || (X) == FourByte_Tr || \
		       (X) == FourByteGeno)


#define has_sexp_coding(X) (((X) >= TwoBitGeno && (X) <= ThreeBit) ||	\
  ((X) >= DotFile && (X) <= FileDot) ||	(X) == OneByteHaplo ||	    \
			    (X) == OneByteGeno || (X) == FourByteGeno)

#define has_lda_in_Byte(X) ((X) == OrigPlink || (X) == OrigPlinkTransposed)

#define PLINK2HUMAN(X) ((X) <= 1 ? 0 : ((X) - 1))
#define HUMAN2PLINK(X) ((X) == 0 ? 0 : (X) == 1 ? CODE_ONE : (X) == 2 ? 3 : 1)
#define TWOBIT2HUMAN(X) ((X) <= 2 ? (X) : 0)
#define HUMAN2TWOBIT(X) ((X) <= 2 ? (X) : 0)


#endif
 
