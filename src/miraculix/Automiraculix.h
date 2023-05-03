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


#ifndef AutoMiraculix_H
#define AutoMiraculix_H 1


#define CURRENT_IMPLEMENTATION 7

// all are _4Byte ordered except OneByteGeno/OneByteHaplo which are on byte level
// So within _4Byte the big/little endian store differently, but not among _4Byte
// Hence, if little endian, the ordering in genoVector of geno and Vector do
// not match in the RAM, but only after 32-bit-loading.


// Ordering:  snps, indiv, lda, coding, variant, ldabitalign

typedef enum coding_type {AutoCoding,  // 0K
   OneBitGeno, // (256 or) 512 bit alignement
   TwoBitGeno, // (256 or) 512 bit alignement
   ThreeBit,
   Plink, // only uised in VectorGeno
   //  00 -> 00
   //  01 : missing
   //  10 -> 01
   //  11 -> 10
   FiveCodes,  // 5
   OrigPlink, // here, lda in Bytes!!
   
   unused7,
   unused8,
   unused9,

   FourBit, // format implicitely used on GPU
   OneByteGeno, // 8; 0,1,2 in 1 Bytes. This is a mixture between the
   // user defined unaligned 32-bit representation & an sse representation
   // The leading orderering of the data is through vector/matrix _1Byte[].
   // This implies that on big endian machines, the data must be reordered
   // for the leading 32bit storing _4Byte for SSE.
   // it has a standard 256-bit alignment
   TwoByte, // unused
   FourByteGeno, // 13--- 4-Byte-(R i.e. non)-alignment
   
   FourByte_Tr, // to indicate that matrix contains half haplo, transposed;
   //             this type is not processed! Non-transposed half-haplo
   //             is treated as geno-type
   unused15, unused16, 
   OneBitHaplo,  // 
   TwoBitHaplo,
   
   OneByteHaplo,
   
   FourByteHaplo,
   FourByteSingleBit, // 21, IndividualsPerColumn=TRUE, DoubledIndividuals=TRUE
   EightByteHaplo, // IndividualsPerColumn=TRUE, DoubledIndividuals=FALSE
   FourByteSingleBit_Tr, // IndividualsPerColumn=FALSE, DoubledIndividuals=TRUE
   EightByteHaplo_Tr, // IndividualsPerColumn=FALSE, DoubledIndividuals=FALSE
   DotFile,// 25 no coding, but read from the file
   FileDot,

   // technical Codings start:
   CorrespondingGeno, // not a coding on its own
   FiveCodesTransposed, // must be first transposed coding
   TwoBitGenoTransposed,
   OrigPlinkTransposed,// here, lda in Bytes!!
   PlinkTransposed,
  
   UnknownSNPcoding // 32
} coding_type;

#define nr_coding_type (UnknownSNPcoding + 1)

#define FirstUserCoding AutoCoding
#define LastUserCoding EightByteHaplo_Tr // with many gaps inbetween, see options.cc

#define FirstMoBPScoding TwoBitGeno
#define LastMoBPScoding TwoBitGeno

#define FirstUserGenoCoding OneByteGeno // used by user
#define LastUserGenoCoding FourByteGeno
#define FirstGenoCoding OneBitGeno
#define LastGenoCoding FourByteGeno
#define LastTransposedGenoCoding PlinkTransposed
#define LastTechnicalGenoCoding LastTransposedGenoCoding

#define FirstUserHaploCoding FourByteHaplo // used by user
#define LastUserHaploCoding EightByteHaplo_Tr
#define FirstHaploCoding OneBitHaplo
#define LastHaploCoding EightByteHaplo_Tr

 

// Haplo nach geno!!


// CHECK haplogeno.cc VARIANTS when changing!
#define VARIANT_DEFAULT 0
#define VARIANT_GPU 768 // never change ! see options.cc
#define VARIANT_R 896 // dito
#define VARIANT_A 32 // see options.cc, variantsSHR5
#define VARIANT_B 64 
#define VARIANT_C 96 // maximum!! Otherwise change main_variant

typedef enum centering_type {NoCentering = 0,
  RowMeans=1, 
  ColMeans, User} centering_type;
#define nr_centering_type (User + 1)

typedef enum normalizing_type {NoNormalizing=0, 
  AccordingCentering = 1, 
  GlobalNormalizing}  normalizing_type;
#define nr_normalizing_type (User + 1)



// WHAT just doubles the information that is also available through
// the class information, but easier to access on the C level for
// historical reasons. To be deleted maybe somewhen.


// Coding of the attribute "information"
// !!!! ACHTUNG ZAHLEN MUESSEN DIE GLEICHEN BLEIBEN !!!!


#define IMPLEMENTATION 0
#define SNPS 1 // Wert darf auf keinen Fall geaendert werden
#define INDIVIDUALS 2 // Wert darf auf keinen Fall geaendert werden
#define CODING 3
#define VARIANT 4
#define LDABITALIGN 5
#define LDA 6
#define ORIGINALLY_HAPLO 7
#define MISSINGS 8
#define EXPECTED_REPET_V 9

#define ADDR 10
#define ALIGNADDR 11
#define MEMinUNITS 12 // total genuine memory used to store all matrices --
// neither memory for allignment is not included, nor memory for info
#define ALIGNEDUNITS 13  // memory needed including alignment
#define BIGENDIAN 14

#define CURRENT_SNPS 20 // used in MOBPS

#define FILE_SNPxIND_READ_BYROW 30 // INFO_INDIV_PER_COL, i.e. 'percolumn'
#define FILE_HEADER 31
#define FILE_DOUBLEINDIV 32
#define FILE_LEADINGCOL 33

#define INFO_GENUINELY_LAST FILE_LEADINGCOL

#define RECENTALIGNADDR 45
#define BLOCKEDINFO 46

#define ZAEHLER 62
#define INFO_LAST 63

#define FIRST_FILE_INFO FILE_SNPxIND_READ_BYROW
#define LAST_FILE_INFO FILE_LEADINGCOL 

#define GENOMICMATRIX "genomicmatrix"
#define HAPLOMATRIX "haplomatrix"
#define ORIGINVECTOR "origindata"

#define MAX_LDA_BITALIGN 1024


#define N_MIRACULIX_NAMES_OF_NAMES 4

extern const char *CODING_NAMES[nr_coding_type],
  *NORMALIZING_NAMES[nr_normalizing_type],
  *CENTERING_NAMES[nr_centering_type],
  *INFO_NAMES[INFO_LAST + 1],
  *MIRACULIX_NAMES_OF_NAMES[N_MIRACULIX_NAMES_OF_NAMES];

#endif

