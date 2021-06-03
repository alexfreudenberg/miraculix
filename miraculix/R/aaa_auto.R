# This file has been created automatically by 'rfGenerateConstants'


 ## from  ../../RandomFieldsUtils/RandomFieldsUtils/src/AutoRandomFieldsUtils.h

 MAXUNITS 	<- as.integer(4)
 MAXCHAR 	<- as.integer(18)
 RFOPTIONS 	<- "RFoptions"
 CLASS_TRYERROR 	<- "try-error"

 WARN_UNKNOWN_OPTION_ALL 	<- as.integer(4)
 WARN_UNKNOWN_OPTION_SINGLE 	<- as.integer(3)
 WARN_UNKNOWN_OPTION_CAPITAL 	<- as.integer(2)
 WARN_UNKNOWN_OPTION_NONE1 	<- as.integer(1)
 WARN_UNKNOWN_OPTION_NONE 	<- as.integer(0)



 ## from  src/AutoMiraculix.h

 AutoMiraculix_H 	<- as.integer(1)

 AutoCoding 	<- as.integer(0)
 NoSNPcodingR 	<- as.integer(1)
 ThreeBit 	<- as.integer(2)
 Hamming2 	<- as.integer(3)
 Hamming3 	<- as.integer(4)
 NoSNPcoding 	<- as.integer(5)
 TwoBit 	<- as.integer(6)
 Multiply 	<- as.integer(7)
 Packed 	<- as.integer(8)
 Shuffle 	<- as.integer(9)
 Multiply256 	<- as.integer(10)
 Packed256 	<- as.integer(11)
 Shuffle256 	<- as.integer(12)
 MMAGPU 	<- as.integer(13)
 CaseCount 	<- as.integer(14)
 unused15 	<- as.integer(15)
 unused16 	<- as.integer(16)
 unused17 	<- as.integer(17)
 unused18 	<- as.integer(18)
 unused19 	<- as.integer(19)
 unused20 	<- as.integer(20)
 unused21 	<- as.integer(21)
 unused22 	<- as.integer(22)
 unused23 	<- as.integer(23)
 unused24 	<- as.integer(24)
 unused25 	<- as.integer(25)
 unused26 	<- as.integer(26)
 unused27 	<- as.integer(27)
 unused28 	<- as.integer(28)
 unused29 	<- as.integer(29)
 Haplo 	<- as.integer(30)
 UnknownSNPcoding 	<- as.integer(31)


 nr_snpcoding 	<- as.integer((UnknownSNPcoding+1))
 FirstMoBPSmethod 	<- as.integer(TwoBit)
 LastMoBPSmethod 	<- as.integer(MMAGPU)
 FirstGenuineMethod 	<- as.integer(ThreeBit)
 LastGenuineMethod 	<- as.integer(LastMoBPSmethod)

 CURRENT_VERSION 	<- as.integer(3)

 VERSION 	<- as.integer(0)
 SNPS 	<- as.integer(1)
 INDIVIDUALS 	<- as.integer(2)
 ADDR0 	<- as.integer(3)
 ADDR1 	<- as.integer(4)
 ALIGNADDR0 	<- as.integer(5)
 ALIGNADDR1 	<- as.integer(6)
 SUMGENO 	<- as.integer(7)
 SUMGENO_E9 	<- as.integer(8)
 METHOD 	<- as.integer(9)
 ALIGNMENT 	<- as.integer(10)
 SNPxIND 	<- as.integer(11)
 BITSPERCODE 	<- as.integer(12)
 BYTESPERBLOCK 	<- as.integer(13)
 CODESPERBLOCK 	<- as.integer(14)
 HEADER 	<- as.integer(15)
 DOUBLEINDIV 	<- as.integer(16)
 LEADINGCOL 	<- as.integer(17)
 MEMinUNITS0 	<- as.integer(18)

 MEMinUNITS1 	<- as.integer(19)
 ALIGNEDUNITS0 	<- as.integer(20)
 ALIGNEDUNITS1 	<- as.integer(21)

 UNITSPERINDIV 	<- as.integer(22)
 INFO_GENUINELY_LAST 	<- as.integer(UNITSPERINDIV)
 RECENTALIGNADDR0 	<- as.integer(23)
 RECENTALIGNADDR1 	<- as.integer(24)

 BLOCKEDINFO 	<- as.integer(61)
 ZAEHLER 	<- as.integer(62)
 INFO_LAST 	<- as.integer(63)

 CURRENT_SNPS 	<- as.integer(HEADER)

 GENOMICMATRIX 	<- "genomicmatrix"
 HAPLOMATRIX 	<- "haplomatrix"
 ORIGINVECTOR 	<- "origindata"



 ## from  src/AutoMiraculix.cc

 SNPCODING_NAMES <-
c( "AutoCoding","NoSNPcodingR","ThreeBit","Hamming2","Hamming3","NoSNPcoding","TwoBit","Multiply","Packed","Shuffle","Multiply256","Packed256","Shuffle256","MMAGPU","CaseCount","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","Haplo","unknown" )


 INFO_NAMES <-
c( "version","snps","individuals","addr0","addr1","align0","align1","sumgeno","sumgenoE9","method","alignment","isSNPxInd","bitspercode","bytesperblock","codesperblock","header","DoubledIndividuals","leadingcolumns","memInUnits0","meminUnits1","AlignedUnits0","AlignedUnits1","unitsperindiv","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused" )

