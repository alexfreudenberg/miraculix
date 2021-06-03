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



// gcc -mavx -o hello_avx hello_avx.c

#define MY_METHOD TwoBit

#include "2bit.h"
#include "Bit23.intern.h"
#include "haplogeno.h"
#include "xport_import.h"
#include "align.h"


#define TWO_GENUINEBITSPERMINIBLOCK 65536L // testing only
#define TABLE_SIZE  two_genuineBitsPerMiniblock
#define nr_resultsAND 3
#define nr_resultsOR 4

static BlockType geno_code[nr_genotypes] = {0, 1, 2}, // 0, 1, 2
  result_codeAND[nr_resultsAND] = {0, 1, 2},
  result_codeOR[nr_resultsOR] = {0, 1, 2, 3};
static Uint result_valueAND[nr_resultsAND] = {0, 1, 4},
  result_valueOR[nr_resultsOR] = {0, 0, 0, 2};

static table_type *TABLE2AND = NULL,
  *TABLE2OR = NULL;

Uint BytesPerBlock2() { return BytesPerBlock; }
Uint CodesPerBlock2() { return CodesPerBlock; }
//Uint UnitsPerIndiv2(Uint snps) { return UnitsPerIndiv256(snps); }
// #define UnitsPerIndiv2 UnitsPerIndiv256
#define UPI UnitsPerIndiv256
Uint BitsPerCode2() { return BitsPerCode; }

/*
void static printbits(BlockType x) { //
  printbits(x, sizeof(BlockType) * BitsPerByte, BitsPerCode);//
}
*/

void Init2() {
  if (MiniblocksPerBlock != 4) BUG;
  assert(BytesPerBlock == sizeof(BlockType0));
  initiate_tableI(&TABLE2AND, TABLE_SIZE, CodesPerMiniblock, BitsPerCode,
		  result_codeAND, result_valueAND, nr_resultsAND);
  initiate_tableI(&TABLE2OR, TABLE_SIZE, CodesPerMiniblock, BitsPerCode,
		  result_codeOR, result_valueOR, nr_resultsOR);

}


// these defintions must match extactly as.genomicmatrix
// BODY-CODE IDENTICAL TO CODING3
SEXP matrix_start2 START23(Init2)

  
Ulong sumGeno2(Uint *S, Uint snps, Uint individuals) {	
  Ulong	unitsPerIndiv = UPI(snps),					
    sum = 0L;							
  for (Uint i=0; i<individuals; i++, S += unitsPerIndiv) { /* ok */	
    for (Ulong j=0; j<unitsPerIndiv; j++) {				
      Uint s = S[j];						
      for (Uint u=0; u<CodesPerUnit; u++) {			
	sum += s & CodeMask;					
	s >>= BitsPerCode;					
      }						
    }						  
  }						  
  return sum;						  
}


void coding2(Uint *M, Uint start_individual, Uint end_individual, 
	       Uint start_snp, Uint end_snp, Uint Mnrow, SEXP Ans,
	       double VARIABLE_IS_NOT_USED *G) {
  if (start_snp % CodesPerUnit != 0) BUG; 
  Uint
    *info = GetInfo(Ans),
    cur_snps = end_snp - start_snp,
    *pM = M,
    *endM = M + Mnrow * (end_individual - start_individual),
    allUnits = Units(cur_snps),
    allUnitsM1 = allUnits - 1,
    rest = (cur_snps - allUnitsM1 * CodesPerUnit) * BitsPerCode,
    unitsPerIndiv = UPI(info[SNPS]),
    *code = Align256(Ans, ALIGN_HAPLOGENO) + start_snp / CodesPerUnit +
               start_individual * unitsPerIndiv;

  for (; pM<endM; pM += Mnrow, code += unitsPerIndiv) {
    Uint *mm = pM;
    for (Uint j=0; j<allUnits; j++) {
      Uint C = 0;
      Uint end = j == allUnitsM1 ? rest : BitsPerUnit;
      for (Uint shft = 0; shft < end; mm++, shft+=BitsPerCode)
	C |= *mm << shft;
      code[j] = C;
    }
  }
}


void haplo2geno2(Uint *X, Uint snps, Uint individuals, Uint unitsPerIndiv,
		 Uint *ans) {
  // note that X == ans is allowed !!
  Coding23(0, individuals, 0, snps, unitsPerIndiv, FROMHAPLO);
}


 
void crossprod2(Uint * M, Uint snps, Uint individuals, double *A){
  if (TABLE2AND == NULL) Init2();
  assert(TABLE_SIZE == 65536);
  assert(TABLE2AND != NULL);
  Uint blocks = Blocks(snps),
    unitsPerIndiv = UPI(snps);
  table_type tableAnd[TABLE_SIZE], tableOr[TABLE_SIZE];

  MEMCOPY(tableAnd, TABLE2AND, TABLE_SIZE * sizeof(table_type));
  MEMCOPY(tableOr, TABLE2OR, TABLE_SIZE * sizeof(table_type));
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif	  
  for (Uint i=0; i<individuals; i++) {   
    Uint * Mi = M + i * unitsPerIndiv,
      *Mj = Mi;
    double *ptr = A + i,
      *res_ptr_r = ptr + i * individuals;
    for ( Uint j  =  i; j<individuals; j++, res_ptr_r++, Mj += unitsPerIndiv) {
      Uint sum = 0;
      BlockType *MJ = (BlockType0*) Mj,
	*MI = (BlockType0*) Mi;
      for (Uint s=0; s<blocks; s++, MJ++, MI++) {
	block_compressed x;
	x.x = *MJ & *MI;
	sum+= tableAnd[x.b[0]] + tableAnd[x.b[1]] +
	  tableAnd[x.b[2]] + tableAnd[x.b[3]];

 	x.x = *MJ | *MI;
	sum += tableOr[x.b[0]] + tableOr[x.b[1]]
	  + tableOr[x.b[2]] + tableOr[x.b[3]];
      }
      ptr[j * individuals] = *res_ptr_r = (double) sum;
    }
  }
}


double getValue2(BlockType0 *code, Uint s) {
  return (((Uint *)code)[s / genuineCodesPerUnit] >>
	  ((s % genuineCodesPerUnit) * BitsPerCode)) & CodeMask; 
}

SEXP get_matrix2(SEXP CM) {
  
 SEXP Ans = R_NilValue;
  Uint 
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    allUnits = Units(snps),
    unitsPerIndiv = UPI(snps),    
    *C0 = Align256(CM, ALIGN_HAPLO);
  if (!isGeno(info[METHOD])) ERR("not a geno coding");
   
  PROTECT(Ans = allocMatrix(INTSXP, snps, individuals));
  Uint *ans = (Uint*) INTEGER(Ans);
  for (Uint i=0, endfor= snps; i<endfor; i++) ans[i] = 0; // noetig?

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)   
#endif
  for (Ulong i=0; i<individuals; i++) {
    Uint *a = ans + i * snps,
      *endans = a + snps,
      *C = C0 + i * unitsPerIndiv;
    for (Uint j=1; j<=allUnits; j++, C++) {
      Uint c = *C,
	*end = j==allUnits ? endans : a + CodesPerUnit;
      for (; a < end;  a++) {
	*a = c & CodeMask;
	c >>= BitsPerCode; 
      }
    }
  }
 
  UNPROTECT(1);
  return Ans;
}


SEXP get_matrixN_2(SEXP CM, SEXP NN) {
    SEXP Ans = R_NilValue;
  //  BlockType first2bits32 = SET32(0x00000003);
  Uint
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    unitsPerIndiv = UPI(snps),
    *C0 = Align256(CM, ALIGN_HAPLO);
  if (!isGeno(info[METHOD]))  ERR("not a genoimatrix coding");
  
  if (lenN ==1) PROTECT(Ans = allocVector(INTSXP, individuals));
  else PROTECT(Ans = allocMatrix(INTSXP, lenN, individuals));
  Uint *M = (Uint*) INTEGER(Ans);

  for (Ulong i=0; i<individuals; i++) {    
    Uint *mm = M + i * lenN,
      *C1 = C0 + i * unitsPerIndiv;
    for (Uint ni=0; ni<lenN; ni++, mm++) {
      Uint nni = N[ni];
      Uint j=nni / CodesPerUnit,
	jj = nni % CodesPerUnit,
	*C2 = C1 + j;
      *mm = (*C2 >> (BitsPerCode * jj)) & CodeMask;
    }
  }
 
  UNPROTECT(1);  
  return Ans;

  
}



SEXP allele_freq2(SEXP SNPxIndiv) {
  if (TABLE2AND == NULL) Init2();
  Uint *M = Align256(SNPxIndiv, ALIGN_23);
  ALLELE_FREQ23(getValue2, UPI);  
}




void zeroNthGeno2(SEXP CM, SEXP NN) {					
  Uint
    first2bits32 = 0x00000003,
    *N = (Uint*) INTEGER(NN),
    lenN = length(NN),
    *info = GetInfo(CM),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    unitsPerIndiv = UPI(snps),
    *C0 = Align256(CM, ALIGN_HAPLO);
  if (!isGeno(info[METHOD]))  ERR("not a genomicmatrix coding");
  
  for (Ulong i=0; i<individuals; i++) {    
    Uint *C1 = C0 + i * unitsPerIndiv;
    for (Uint ni=0; ni<lenN; ni++) {
      Uint j=N[ni] / CodesPerUnit,
	jj = N[ni] % CodesPerUnit,
	*C = C1 + j;
      *C &= ~(first2bits32 << (BitsPerCode * jj));
    }
  }
}

  
static void genoVectorIntern(Uint *CM, double *V, Uint snps,
		      Uint individuals, double *ans) {
  if (OPTIONS_UTILS->basic.kahanCorrection) {
    ERR("'Kahan' currently not programmed.");
     //    genoVectorKahan(V, CM, snps, individuals, ans);
    // return;
  }
  Uint units = Units(snps),
    unitsM1 = units - 1,
    rest = snps - unitsM1 * CodesPerUnit,
    unitsPerIndiv = UPI(snps);
  for (Uint i=0; i<snps; ans[i++] = 0.0);

  // achtung! so nicht mit pragma parallelisierbar !
  for (Ulong i=0; i<individuals; i++) {
    Uint *cm = CM + i * unitsPerIndiv;
    double
      *a = ans,
      f0 = V[i],
      factor[4] = {0.0, f0, 2.0 * f0, RF_NA}; // nur fuer geno mit 00, 01, 10 Code

    for (Uint b=0; b<=unitsM1; b++) {
      Uint c = cm[b],
	endfor = b < unitsM1 ? CodesPerUnit : rest;      
      for (Uint j=0; j<endfor; j++) {
	*(a++) += factor[c & CodeMask];
	c >>= BitsPerCode;
      }      
    }
  }
}



void genoVector2(SEXP Z, SEXP V, double *ans) {
  Uint
    *info = GetInfo(Z),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    len = length(V);
  if (len != individuals) ERR("vector 'V' not of correct length");
  genoVectorIntern(Align256(Z, ALIGN_SHUFFLE), REAL(V), 
		    snps, individuals, ans);
}


void vectorGeno2(SEXP V, SEXP Z, double *ans) {
  Uint
    type = TYPEOF(V),
    *info = GetInfo(Z),
    snps = info[SNPS],
    individuals = info[INDIVIDUALS],
    blocks = Blocks(snps),
    rest = snps - (blocks - 1) * CodesPerBlock,
    len = length(V),
    unitsPerIndiv = UPI(snps),
    *CM = Align256(Z, ALIGN_VECTOR);
  if (len != snps) ERR("vector 'V' not of correct length");
  double *Vd = NULL;
  int *Vi = NULL;
  if (type == INTSXP) Vi = INTEGER(V);
  else if (type == LGLSXP) Vi = LOGICAL(V);
  else if (type == REALSXP) Vd = REAL(V);
  else ERR("type of V not numeric.");


#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES)
#endif	  
  for (Ulong i=0; i<individuals; i++) {
    BlockType *cm = (BlockType0*) (CM + i * unitsPerIndiv);
    double sum0 = 0.0,
      *vi = Vd;
    int *vd = Vi;
    for (Uint j=1; j<=blocks; j++, cm++) {
      BlockType cm0 = *cm;
      Uint end = j == blocks ? rest : CodesPerBlock;
      if (type == REALSXP)
	for (Uint k=0; k<end; k++, cm0 >>= BitsPerCode, vd++)
	  sum0 += *vd * (double) (cm0 & CodeMask);
      else 
	for (Uint k=0; k<end; k++, cm0 >>= BitsPerCode, vi++)
	  sum0 += (double) (*vi * (cm0 & CodeMask));
    }
    ans[i] = sum0;
  }
}
