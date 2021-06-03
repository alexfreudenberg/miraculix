
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2019  Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/



#define BitsPerCode 2L


#include "IntrinsicsBase.h"
#include "intrinsics.h"
#include <General_utils.h>
#include "xport_import.h"
#include "haplogeno.h"
#include "utils.h"

#define MY_METHOD Haplo
#include "align.h"

#define BitsPerCode 2L
static Uint rev_haplo_code[4] = {0L, 1L, 1L, 2L};


Uint BytesPerBlockHaplo() { return BytesPerBlock; }
Uint CodesPerBlockHaplo() { return CodesPerBlock; }
#define UPI UnitsPerIndiv256
Uint BitsPerCodeHaplo() { return BitsPerCode; }

void assert_haplo() {
}

void assert_haplo(SEXP CM) {
  Uint *info = GetInfo(CM);
  if (info[METHOD] != Haplo) ERR("not a haplotype matrix");
}


SEXP create_codevectorHaplo(Uint snps, Uint individuals) {
  assert_haplo();
  return CreateEmptyCodeVector(snps, individuals, MY_METHOD);
}

SEXP create_codevectorHaplo(Uint snps, Uint individuals,
			     SEXP VARIABLE_IS_NOT_USED X) {
  return create_codevectorHaplo(snps, individuals);
}

#define Inner(mm12)							\
  (double *freq1, double *freq2, Uint *info, Uint *code) {		\
    Ulong snps = info[SNPS];						\
    Uint individuals = info[INDIVIDUALS],				\
      unitsPerIndiv = UPI(snps),				\
      allUnits = Units(snps);						\
    for (Uint i=0; i<individuals; i++, code += unitsPerIndiv) {	/* ok */ \
      Uint mm = 0;							\
      for (Uint j=0; j<allUnits; j++) {					\
	Uint dummy = 0;							\
	for (Uint shft=0; shft<BitsPerUnit && mm < snps; mm++) {	\
	  double f2 = freq2[mm];					\
	  mm12;								\
	  dummy |= mm1 << shft++;					\
	  dummy |= mm2 << shft++;					\
	}								\
	code[j] = dummy;						\
      }									\
    }									\
  }

void InnerDetermHaplo Inner(Uint mm1 = freq1[mm] == 1.0;	 
			    Uint mm2 = ISNA(f2) ? mm1 : (f2  == 1.0))

void InnerRandomHaplo Inner(Uint mm1 = UNIFORM_RANDOM <= freq1[mm];
			    Uint mm2 = ISNA(f2) ? mm1 : UNIFORM_RANDOM <= f2)

SEXP rhaplomatrix(SEXP Freq1, SEXP Freq2, SEXP Individuals) {
  Uint individuals = INTEGER(Individuals)[0],
    snps = length(Freq1);
  if (snps != (Uint) length(Freq2))
    ERR("'freq' and 'freq2' do not have the same length");
  return create_codevectorHaplo(snps, individuals);
}

  
SEXP rhaplomatrixPart2(SEXP Freq1, SEXP Freq2, SEXP Code) {
  double *freq1 = REAL(Freq1),
    *freq2 = REAL(Freq2);
  Uint
    *info = GetInfo(Code),
    snps = length(Freq1),
    *C =  Align256(Code, ALIGN_HAPLO);
  assert(snps == info[SNPS]);
 
  bool random = false;  
  for (Uint i=0; i < snps; i++) {
    double f = freq1[i];
    if (f < 0 || f > 1) ERR("frequencies not in [0,1]");
    random |= f != 0 && f != 1;
    f = freq2[i];
    if (!ISNA(f)) {
      if (f < 0 || f > 1) ERR("frequencies not in [0,1]");
      random |= f != 0 &&  f != 1;
    }
  }

  if (random) {
    GetRNGstate();
    InnerRandomHaplo(freq1, freq2, info, C);
    PutRNGstate();
  } else {
    InnerDetermHaplo(freq1, freq2, info, C);
  }

  return Code;
}

void codeInnerHaplo2(Uint *M1, Uint *M2, SEXP Code){ //Uint snps, Uint *code) {
  Uint j,
    *info = GetInfo(Code),
    snps = info[SNPS],
    *code = Align256(Code, ALIGN_HAPLO),
    *endM1 = M1 + snps,
    allUnits = Units(snps);					
  for (j=0; j<allUnits; j++) {
    Uint  dummy = 0;
    for (Uint shft=0; shft<BitsPerUnit && M1 < endM1; M1++, M2++) {
      assert(*M1 == 0 || *M1 == 1);
      assert(*M2 == 0 || *M2 == 1);
      dummy |= *M1 << shft++;			       
      dummy |= *M2 << shft++;
    }
    code[j] = dummy;
  }
}


SEXP codeHaplo2(SEXP M1, SEXP M2) {// 2 x Haplo - Matrix, as vector gespeichert
  if (length(M1) == 0 || length(M1) != length(M2)) {
    if (length(M1) != length(M2))
      ERR2("the two haplotype informations have unequal length (%d and %d)",
	  length(M1), length(M2))
    else ERR("both haplotype informations have length 0.");
  }
  if (!isVector(M1) || !isVector(M2)) ERR("'M' or 'N' is not a vector.");
  Uint snps = (Uint) length(M1);  
  SEXP Code;
  PROTECT(Code = create_codevectorHaplo(snps, 1));
  ToInt(M1);  
  ToInt(M2);  
  codeInnerHaplo2(M1int, M2int, Code); // snps,  Align256(Code, ALIGN_HAPLO)); 
  FREEint(M1);
  FREEint(M2);
  UNPROTECT(1);
  return Code;
}




SEXP GetHaploSequence(Uint snps, Uint individuals,
		      bool indivpercol, bool doubledIndiv,
		      Uint currentBits, // BitsPerCode
		      Uint *unitIncr, Uint *bitIncr,
		      Uint *delta, Uint *indivIncr, bool create){
  SEXP Ans = R_NilValue;
  
  if (indivpercol) {
    *indivIncr = snps * currentBits;
    if (doubledIndiv) {
      *unitIncr = CodesPerUnit; // X
      *delta = currentBits == 1 ? 0L : snps ; // X; pointing to next bit
      *bitIncr = 1L; // X
      if (create) Ans = allocMatrix(INTSXP, snps, individuals * currentBits);
   } else {
      *unitIncr = BitsPerUnit * currentBits / BitsPerCode;
      *delta = currentBits == 1 ? 0L : 1L;
      *bitIncr = currentBits;
      if (create) Ans = allocMatrix(INTSXP, (Ulong) currentBits * snps, individuals);
    }
  } else {
    *unitIncr = BitsPerUnit * individuals * currentBits / BitsPerCode;
    *bitIncr = currentBits * individuals;
    if (doubledIndiv) {
      *delta = currentBits == 1 ? 0L : 1L;
      *indivIncr = currentBits;
      if (create) Ans = allocMatrix(INTSXP, currentBits * individuals, snps);
    } else {
      *delta = currentBits == 1 ? 0L : individuals;
      *indivIncr = 1L;
      if (create) Ans=allocMatrix(INTSXP,  individuals, (Ulong) currentBits * snps);
    }
  }
  return Ans;
}




void codeInnerHaplo(Uint *MM, bool indivpercol, bool doubledIndiv, SEXP Code) {
  Uint unitIncr, bitIncr, delta, indivIncr,
    *code = Align256(Code, ALIGN_HAPLO),
    *info = GetInfo(Code),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],   
    //    bits = snps * BitsPerCode,
    unitsPerIndiv = UPI(snps),
    allUnits = Units(snps),
    allUnitsM1 = allUnits -1,
    rest = (snps - allUnitsM1 * CodesPerUnit) * BitsPerCode;

  GetHaploSequence(snps, individuals, indivpercol, doubledIndiv, BitsPerCode,
		   &unitIncr, &bitIncr, &delta, &indivIncr, false);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif
  for (Ulong i=0; i<individuals; i++) {
    Uint *cm = code + i * unitsPerIndiv,
      *M = MM + i * indivIncr;
    for (Ulong j=0; j<allUnits; j++) {
      Uint
	dummy = 0,
	*mm = M + j * unitIncr,
	*mm2 = mm + delta,
	end = j==allUnitsM1 ? rest : BitsPerUnit;
      for (Uint shft=0; shft < end; mm+=bitIncr, mm2 += bitIncr) {
	assert(*mm == 0 || *mm == 1);
	assert(*mm2 == 0 || *mm2 == 1);
	dummy |= *mm << shft++;
	dummy |= *mm2 << shft++;
      }
      cm[j] = dummy;
    }
  }
}

SEXP codeHaplo(SEXP M, SEXP IndivPerCol, SEXP DoubledIndiv) {// 2 x Haplo - Matrix, as vector gespeichert
  if (length(M) == 0) ERR("'M' has length 0.");
  Uint snps, individuals,
    nrow = (Uint) nrows(M),
    ncol = (Uint) ncols(M);
  bool doubledIndiv = LOGICAL(DoubledIndiv)[0],
    indivpercol = LOGICAL(IndivPerCol)[0];

  if (indivpercol) {
    snps = nrow;
    individuals = ncol;
  } else {
    snps = ncol;
    individuals = nrow;
  }
  if (doubledIndiv) {
    if (individuals % BitsPerCode != 0)
      ERR("information on individuals not doubled");
    individuals /= BitsPerCode;
  } else {
    if (snps % BitsPerCode != 0) ERR("information on haplotype odd");
    snps /= BitsPerCode;
  }

  SEXP Code;
  PROTECT(Code = create_codevectorHaplo(snps, individuals));
  ToInt(M);
  codeInnerHaplo(Mint, indivpercol, doubledIndiv, Code);

  // Uint *info = GetInfoUnchecked(Code);
  //  Uint *address = (Uint *) INTEGER(Code),
    //    *algnaddress = (Uint *) algn_general(INTEGER(Code), info[ALIGNMENT]),
    //    *infoaddress = (Uint*) GetAddress(info, ADDR0);
  //  int units = info[ UNITSPERINDIV ] * info[INDIVIDUALS];
  //  printf("**** %d: ", units);  for (int i=0; i<units; i++) printf("%d ", algnaddress[i]); printf("\n");

  FREEint(M);
  UNPROTECT(1);
  return Code;
}

void InitGetHaplo(SEXP CM, Uint **code, Uint *unitsPerIndiv) {
  assert_haplo(CM);
  Uint 
    *info = GetInfo(CM),
    snps = info[SNPS];
  *unitsPerIndiv = UPI(snps);
  *code = Align256(CM, ALIGN_HAPLO);
}

#define ByteCodeMask 3

Uint GetHaplo(Uint *cm, Uint snp) {
  return rev_haplo_code[(((char *) cm)[snp / CodesPerByte] >> ((snp % CodesPerByte) * BitsPerCode)) & ByteCodeMask];
}

Uint GetHaploX(Uint *cm, Uint snp) {
  return (((char *) cm)[snp / CodesPerByte] >> ((snp % CodesPerByte) * BitsPerCode)) & ByteCodeMask;
} 

SEXP decodeHaplo(SEXP CM, SEXP Individuals,
		 SEXP Sets, SEXP IndivPerCol, SEXP DoubledIndiv) {
  assert_haplo(CM);
  bool IndivGiven = length(Individuals) > 0;
  Uint unitIncr, bitIncr, delta, indivIncr,
    currentBits = length(Sets),
    *info = GetInfo(CM),
    snps = info[SNPS], 
    individuals = info[INDIVIDUALS],
    currentIndiv = IndivGiven ? length(Individuals) : individuals,
    *I = IndivGiven ? (Uint*) INTEGER(Individuals) : NULL,
    //   bits = snps * BitsPerCode,
    allUnits = Units(snps),
    allUnitsM1 = allUnits - 1,
    rest = snps - allUnitsM1 * CodesPerUnit,
    *code = Align256(CM, ALIGN_HAPLO);
  Ulong unitsPerIndiv =UPI(snps);
  bool doubledIndiv = LOGICAL(DoubledIndiv)[0],
    indivpercol = LOGICAL(IndivPerCol)[0];
  if (info[METHOD] != Haplo) ERR("not a haplotype matrix");
  if (currentBits != 1 && currentBits != BitsPerCode)
    ERR("The number of chromosome sets can only be 1 or 2.");
  bool set1 = INTEGER(Sets)[0] == 1;
  if (currentBits == 2 && (!set1 || INTEGER(Sets)[1] != 2))
    ERR("The sets must be given in the order 1:2")
  else if (!set1 && INTEGER(Sets)[0] != 2)
    ERR("The value for a chromosome set can be 1 and 2 only.");
 
  SEXP Ans;
  PROTECT(Ans = GetHaploSequence(snps, currentIndiv, indivpercol, doubledIndiv,
				 currentBits, &unitIncr, &bitIncr, &delta,
				 &indivIncr, true));
  
  Uint *MM = (Uint*) INTEGER(Ans);
  Ulong total = (Ulong) currentBits * snps * currentIndiv;
  for (Ulong i=0; i<total; i++) MM[i] = 0; // noetig?
  assert(BitsPerCode == 2);

  assert(snps > 0);
  assert(individuals > 0);
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Ulong i=0; i<currentIndiv; i++) {
    Uint *cm = code + unitsPerIndiv * (IndivGiven ? I[i] - 1 : i),
      *M = MM + i * indivIncr;
#if defined SCHLATHERS_MACHINE
    Uint snp = 0;
#endif    
    for (Ulong j=0; j<allUnits; j++) {
      Uint
	C = cm[j],
	*mm = M + j * unitIncr,
	*mm2 = mm + delta,
 	shift = 1L,
	end = j == allUnitsM1 ? rest : CodesPerUnit;
      for (Uint u=0; u<end; u++, mm+=bitIncr, mm2 += bitIncr) {
	// ordering important !!! as mm2 might be overwritten by mm!!
	*mm2 = (C & (shift << 1)) > 0; // never change as mm could be mm2!!
	if (set1) *mm = (C & shift) > 0;
#if defined SCHLATHERS_MACHINE
	assert(delta == 0 || GetHaploX(cm, snp++) == *mm + 2 * *mm2);
#endif
 	shift <<= 2;
    }
    }
  }

  UNPROTECT(1);
  return Ans;
}



void zeroNthHaplo(SEXP CM, SEXP NN) { // 0.45 bei 25000 x 25000
  assert_haplo(CM);
  Uint
    first2bits32 = 0x00000003,
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    unitsPerIndiv =  UPI(snps),
    *C0 = Align256(CM, ALIGN_HAPLO);
  if (info[METHOD] != Haplo) ERR("not a haplotype matrix");
  
  for (Ulong j=0; j<individuals; j++) {    
    Uint *C1 = C0 + j * unitsPerIndiv;
    for (Uint i=0; i<lenN; i++) {
      Uint k=N[i] / CodesPerUnit,
	n = N[i] % CodesPerUnit;
      C1[k] &= ~(first2bits32 << (BitsPerCode * n));
    }
  }
}

