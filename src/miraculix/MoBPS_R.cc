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
#define MY_VARIANT 32 // 


#include "Basic_miraculix.h"

#if defined compatibility_to_R_h


#include "intrinsics_specific.h"
#include "xport_import.h"
#include "options.h"
#include "haplogeno.h"
#include "miraculix.h"

#include "MX.h"
#include "Template.h"
#include "Vector.matrix.h"
#include "utils_miraculix.h"
#include "2bit.h"
#include "transform.h"
#include "intrinsicsCheck.h"

//#define DO_FLOAT 1  // must be completely re-checked if defined


//
#define show(X)
// void show(const char *s) {  if (debugging) PRINTF("%.50s", s); }

 

Uint BitMaskStart[CodesPerUnit], BitMaskEnd[CodesPerUnit]; // of PlainInteger64!

bool MoBPSNotInit = true;
void InitMoBPS() {
  if (!MoBPSNotInit) BUG;
  MoBPSNotInit = false;
  Uint Einsen = 0xFFFFFFFF;
  BitMaskEnd[0] = 0x03;

  for (int k=1; k<CodesPerUnit; k++) {
    //    XOR(BitMaskStart[k], BitMaskEnd[k-1], Einsen); // not
    BitMaskStart[k] = BitMaskEnd[k-1] ^ Einsen; // not
    BitMaskEnd[k] = BitMaskEnd[k-1];
    BitMaskEnd[k] |= 0x03 << (BitsPerCode * k);
  }
  BitMaskStart[0] = Einsen;
}


coding_type assert_MoBPS() {
  if (MoBPSNotInit) InitMoBPS();
  option_type *global = &(KEYT_M()->global);
  coding_type coding = global->genetics.coding;
  if (coding != TwoBitGeno) 
    ERR1("MoBPS works with 'TwoBitGeno', not with '%.20s'\n",
	 CODING_NAMES[coding]);
  return coding;
}


#define ORIGIN_GENERATION 0
#define ORIGIN_SEX 1
#define ORIGIN_NR 2
#define ORIGIN_HAPLO 3

#define BITS_GENE_INPUT 6
#define BITS_SEX 1
#define BITS_INDIVIDUALS 22
#define BITS_HAPLO 3

#define MAX_GENE_INPUT (1 << BITS_GENE_INPUT)
#define MAX_SEX (1 << BITS_SEX)
#define MAX_INDIVIDUALS (1 << BITS_INDIVIDUALS)
#define MAX_HAPLO (1 << BITS_HAPLO)

#define PATTERN_GENE_INPUT  (MAX_GENE_INPUT - 1)
#define PATTERN_SEX (MAX_SEX - 1)
#define PATTERN_INDIVIDUALS (MAX_INDIVIDUALS - 1)
#define PATTERN_HAPLO (MAX_HAPLO - 1)


void static inline decodeOrigin(Uint x, Uint *ans) {
  ans[ORIGIN_HAPLO] = (x & PATTERN_HAPLO);
  x >>= BITS_HAPLO;
  ans[ORIGIN_NR] = (x & PATTERN_INDIVIDUALS);
  x >>= BITS_INDIVIDUALS;  
  ans[ORIGIN_SEX] = (x & PATTERN_SEX);
  x >>= BITS_SEX;
  ans[ORIGIN_GENERATION] = x;
}


SEXP decodeOrigins(SEXP M, SEXP Line) {  
  assert_MoBPS();
  Uint line = Int0(Line);
  assert(line > 0); // R index
  SEXP Ans;
  PROTECT(Ans = allocVector(INTSXP, 4));
  Uint *ans = (Uint*) INTEGER(Ans);
  decodeOrigin(Inti(M, line - 1), ans);
  for (int j=0; j<4; ans[j++]++);
  UNPROTECT(1);
  return Ans;
}


SEXP codeOrigins(SEXP M) {
  assert_MoBPS();
  int nrow = nrows(M);
  assert(ncols(M) == 4);
  
  SEXP Ans;
  PROTECT(Ans = allocVector(INTSXP, nrow));
 
  Uint *ans = (Uint*) INTEGER(Ans);
  for (int j=0; j<nrow; j++) ans[j] = 0;  // unnoetig?!

#define DO1(TYPE, RTYPE)			\
  TYPE *gen = RTYPE(M),				\
    *sex = gen + nrow,				\
    *nr = sex + nrow,				\
    *haplo = nr + nrow					       

#ifdef SCHLATHERS_MACHINE  
#define CONTROL \
    Uint control[4];							\
    decodeOrigin(ans[j], control);					\
    if (control[0]+1 != g || control[1]+1 != s || control[2]+1 != n ||	\
	control[3]+1 != H) {						\
      /*	p rintf("i=%d: generation:%d->%d sex:%d->%d nr:%d->%d haplo:%d->%d\n", i, g, control[0], s, control[1], n, control[2], H, control[3]); */ \
      BUG;								\
    }
#else
  #define CONTROL
#endif  
  
#define DO2								\
  for (int j=0; j<nrow; j++) {						\
    Uint g = (Uint) gen[j],						\
      s = (Uint) sex[j],						\
      n = (Uint) nr[j],							\
      H = (Uint) haplo[j];						\
    if (g < 1 || g > MAX_GENE_INPUT|| s < 1 || s > MAX_SEX ||		\
	n < 1 || n > MAX_INDIVIDUALS ||	H < 1 || H > MAX_HAPLO)		\
      ERR1("some value in row %d out of bound", j + 1);			\
    ans[j] = ((((((g - 1) << BITS_SEX) +				\
		 s - 1) << BITS_INDIVIDUALS) +				\
	       n - 1) << BITS_HAPLO) +					\
      H - 1;								\
    CONTROL								\
  }

  
  if (TYPEOF(M) == REALSXP) {
    DO1(double, REAL);
    DO2   
      } else if (TYPEOF(M) == INTSXP) {
    DO1(int, INTEGER);
    DO2   
      } else BUG;
  
  UNPROTECT(1);
  return Ans;
}


// population
#define INFO 0
#define BREEDING 1

// info
#define INFO_INFOHAPLO 0
#define NSNPS 2
#define POSITION 3
#define SNP_POSITION 5
#define LENGTH_CHROM 6
#define LENGTH_TOTAL 7
#define EQUIDIST_SNPS 17
#define GENERATION_REFTABLE 18

// individual
#define RECOMBI1 0
#define RECOMBI2 1
#define MUTAT1 2
#define MUTAT2 3
#define ORIGIN1 4
#define ORIGIN2 5
//#define FATHER 6
//#define MOTHER 7
#define HAPLOTYP1 8 // in case of founder, otherwise empty
#define HAPLOTYP2 9
#define CODEDHAPLO HAPLOTYP1

#define MAX_CHROM 100


#define HaploCodeMask 0x00000001
static inline Uint GetTwoBitHaploCode(Uint *M, Uint S, Uint haplo) {
  assert(haplo <= 1);
  Uint idx = S % CodesPerUnit; 
  return (M[S / CodesPerUnit] >> (idx * BitsPerCode + haplo)) & HaploCodeMask;
}


void static inline PutTwoBitHaploOne(Uint pos, Uint haplo, Uint *code) {
  Uint idx = pos % CodesPerUnit;
  Uint *c = code + pos / CodesPerUnit; 
  *c |= HaploCodeMask << (idx * BitsPerCode + haplo);
}


  
Uint getNrSNPposition(double cur_recomb, double *position, Uint Min, Uint Max) {
  Uint min = Min,
    max = Max;
  //  print("%d %d\n", min, max);
  while (min < max) {
    Uint inbetween = (min + max + 1) / 2;
    if (position[inbetween] <= cur_recomb) min = inbetween;
    if (position[inbetween] > cur_recomb) max = inbetween - 1;
  }
  if (min > max) BUG;
  if (position[max] > cur_recomb) {
    if (max == Min) return max - 1; // links des ersten Markers auf dem
    //                                 Chromosom; somit wird letzter marker 
    //                                 vom Vor-Chromoosom zugeordnet
    BUG;
  }
  if (max < Max && position[max+1] <= cur_recomb) BUG;
  return max;
}


#define ASSERT(WHERE,WHAT,NAME)						\
  assert(STRCMP(CHAR(STRING_ELT(getAttribPointer(WHERE, R_NamesSymbol), WHAT)),\
		NAME) == 0);
SEXP IcomputeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		  int From_SNP, int To_SNP, SEXP Select, KEY_type *KT,
		  unit_t *ans){
  Uint
    cur_select = NA_INTEGER,
    individuals = LENGTH(Generation),
    len_sel = LENGTH(Select);
  bool do_select = len_sel > 0;
  if (individuals != (Uint) LENGTH(Sex) || individuals != (Uint) LENGTH(Nr))
    ERR0("'gen', 'gender' and 'nr' do not have same length");
  if (From_SNP <= 0) ERR0("'from_p' not valid.");

  ToIntLocal(Select);
   
  ASSERT(population, BREEDING, "breeding");
  SEXP breeding = VECTOR_ELT(population, BREEDING);

  ASSERT(population, INFO, "info")
    SEXP info = VECTOR_ELT(population, INFO);
  Uint len_info = LENGTH(info);

  static Uint gen_ref[MAX_GENE_INPUT] = { 0 };
  if (GENERATION_REFTABLE >= len_info ||
      LENGTH(VECTOR_ELT(info, GENERATION_REFTABLE)) == 0) {
    if (gen_ref[MAX_GENE_INPUT - 1] == 0) { // should never happen except uninit
      for (Uint j=0; j<MAX_GENE_INPUT; j++) gen_ref[j] = j;
      PRINTF("Generation reference table now initialized.\n"); 
    }
  } else {
    ASSERT(info, GENERATION_REFTABLE, "origin.gen");
    Uint *p = (Uint*) INTEGER(VECTOR_ELT(info, GENERATION_REFTABLE)); 
    for (Uint j=0; j<MAX_GENE_INPUT; j++) gen_ref[j] = (Uint) (p[j] - 1);
  }

  ASSERT(info, LENGTH_TOTAL, "length.total");
  Long n_chrom = LENGTH(VECTOR_ELT(info, LENGTH_TOTAL)) - 1;
  if (n_chrom >= MAX_CHROM) ERR1("maximum chromosoms is %d", MAX_CHROM);
  double *length_total = REAL(VECTOR_ELT(info, LENGTH_TOTAL));
  
  ASSERT(info, LENGTH_CHROM, "length");
  double *LengthChrom = REAL(VECTOR_ELT(info, LENGTH_CHROM));
  
  ASSERT(info, NSNPS, "snp");
  Uint n_snps,
    n_protect = 0;
  SEXP InfoHaplo = VECTOR_ELT(info, INFO_INFOHAPLO);
  coding_type coding = assert_MoBPS();
  

  if (TYPEOF(InfoHaplo) != INTSXP || LENGTH(InfoHaplo) != INFO_LAST + 1) {
    PROTECT(InfoHaplo = allocVector(INTSXP, INFO_LAST + 1));
    n_protect++;
    int *infoHaplo = INTEGER(InfoHaplo);
    for (int i=0; i<=INFO_LAST; infoHaplo[i++] = NA_INTEGER);
    infoHaplo[CODING] = TwoBitHaplo;
    infoHaplo[VARIANT] =
      main_variant(check_variant((coding_type) infoHaplo[CODING],
				 VARIANT_DEFAULT,
				 KT->global_utils.basic.efficient));
    infoHaplo[LDABITALIGN] = MY_LDABITALIGN;
    SET_VECTOR_ELT(info, INFO_INFOHAPLO, InfoHaplo);
  }
   
  double length_prec_chroms[MAX_CHROM + 1];
  SEXP snps_perchrom = VECTOR_ELT(info, NSNPS);
  ToIntLocal(snps_perchrom);
  Uint
    totalsnps_prec_chroms[MAX_CHROM + 1];
  length_prec_chroms[0] = 0;
  totalsnps_prec_chroms[0] = 0;
  for (Uint c=1; c<=n_chrom; c++) {
    length_prec_chroms[c] =  length_prec_chroms[c-1] + LengthChrom[c-1];
    totalsnps_prec_chroms[c] = totalsnps_prec_chroms[c-1] +
      snps_perchromint[c-1];
  }

  n_snps = totalsnps_prec_chroms[n_chrom];
  assert(n_snps > 0);
  int fromSNP = From_SNP - 1,
    toSNP = To_SNP==NA_INTEGER || To_SNP >= (int) n_snps ? MAXINT : To_SNP -1;
  if (toSNP < fromSNP) ERR0("'to_p' smaller than 'from_p'");
  if (fromSNP % CodesPerUnit != 0)
    ERR1("(from_p-1) must be devisible by %u", CodesPerUnit);
  if (toSNP < MAXINT && (toSNP + 1)  % CodesPerUnit != 0)
    ERR1("'to_p' must be devisible by %u or 'NA'", CodesPerUnit);


  Long cur_nsnps = 0;
  if (do_select) {
    cur_select = Selectint[0] - 1;
    int toSNP_P1 = toSNP + 1;
    for (Uint k=0; k<len_sel; k++)
      cur_nsnps += Selectint[k] > (Uint) fromSNP && Selectint[k] <= (Uint) toSNP_P1;
  } else cur_nsnps = 1L + MIN(n_snps - 1L, (Uint) toSNP) - fromSNP;
  int *infoHaplo = INTEGER(InfoHaplo);
  infoHaplo[CURRENT_SNPS] = cur_nsnps;
  if (infoHaplo[LDABITALIGN] == 0) infoHaplo[LDABITALIGN] = MY_LDABITALIGN;
  Long ldAns =  SnpLDA(cur_nsnps, TwoBitHaplo, infoHaplo[LDABITALIGN]);
  infoHaplo[LDA] = ldAns;
  if (ans == NULL) {
    if (n_protect) UNPROTECT(n_protect);
    return InfoHaplo;
  }
 
  ASSERT(info, EQUIDIST_SNPS, "snps.equidistant");
  bool equidist = (bool) LOGICAL(VECTOR_ELT(info, EQUIDIST_SNPS))[0];
  
  ASSERT(info, POSITION, "position");
  SEXP position = VECTOR_ELT(info, POSITION);

  ASSERT(info, SNP_POSITION, "snp.position");
  double *snp_positions = REAL(VECTOR_ELT(info, SNP_POSITION)),
    last_snp_position = snp_positions[n_snps - 1],
    first_snp_position = snp_positions[0];

  //  printf("%d %d %d\n", LENGTH(Generation), LENGTH(Sex), LENGTH(Nr));
  ToIntLocal(Generation); 
  ToIntLocal(Sex);
  ToIntLocal(Nr);
  //  for (int i=0; i<LENGTH(Generation); i++)
  //	 printf("%d %d %d\n", Generationint[i], Sexint[i], Nrint[i]);

  Uint NullEins[2] = { 0x55555555, 0xAAAAAAAA};
  bool savehaplo = KT->global.tuning.savegeno;
  //
  savehaplo = true; 
  
  // int variant =
  check_variant(coding, KT->global.tuning.variant,
		KT->global_utils.basic.efficient);

  for (Long i=0; i<individuals; i++) {
    //    printf("AC i=%d %d\n", i, individuals);
    Uint *haplo = (Uint *) (ans + i * ldAns);    
   //    print("i=%d %d\n", i, individuals);
    if (LENGTH(breeding) < (int) Generationint[i])
      ERR0("value of 'generation' too large.");
    SEXP
      g = VECTOR_ELT(breeding, Generationint[i]-1);
    if (LENGTH(g) < (int) Sexint[i]) ERR0("'sex' out of range.");
    SEXP 
      s = VECTOR_ELT(g, (int) Sexint[i]-1);
    if (LENGTH(s) < (int) Nrint[i]) ERR0("value of 'nr' too large.");  
    SEXP Indiv = VECTOR_ELT(s, Nrint[i]-1);

    //  savehaplo 0 i=0 1; G=1 S=1 N=1

    //  printf("savehaplo %d i=%d %d; G=%d S=%d N=%d CH=%d\n", savehaplo, i, individuals,  Generationint[i], Sexint[i], Nrint[i], CODEDHAPLO);



    bool check = true; // debugging only    
    SEXP origin = VECTOR_ELT(Indiv, ORIGIN1);
    if (savehaplo || origin == Indiv) {
      SEXP coded_haplo = VECTOR_ELT(Indiv, CODEDHAPLO);
      Long bytes = ldAns * BytesPerUnit;
      if (coded_haplo == R_NilValue) {
	Long mem = calculateAlignedMem(ldAns, coding);
	coded_haplo = PROTECT(allocVector(INTSXP, mem));
	MEMCOPY(algn_generalS(coded_haplo, MY_LDABITALIGN), haplo, bytes);
	SET_VECTOR_ELT(Indiv, CODEDHAPLO, coded_haplo);
	UNPROTECT(1);
      } else {
	if (check) {
	  PROTECT(coded_haplo);
	  if (MEMCMP(algn_generalS(coded_haplo, MY_LDABITALIGN), haplo, bytes))
	    BUG;
	  UNPROTECT(1);
   	}
      }
    } else if (VECTOR_ELT(Indiv, CODEDHAPLO) != R_NilValue) BUG;

    for (Uint hapnr=0; hapnr<=1; hapnr++) {
      if (TYPEOF(VECTOR_ELT(Indiv, RECOMBI1 + hapnr)) != REALSXP)
	ERR0("real valued vector for recombination information expected");
      double *recombi = REAL(VECTOR_ELT(Indiv, RECOMBI1 + hapnr));
      Uint
	n_recomb = LENGTH(VECTOR_ELT(Indiv, RECOMBI1 + hapnr)) - 1;
      origin = VECTOR_ELT(Indiv, ORIGIN1 + hapnr);
      ToIntLocal(origin);
      Uint
	n_mutat = LENGTH(VECTOR_ELT(Indiv, MUTAT1 + hapnr)),
	// separation of the following 3 allows for selection but also
	// for mutation error such as insertion and deletion
	p = 0, // position what is written if there is no selection
	q = 0, // position what is written actually
	psel = 0, pmutat = 0,
	till = 0,
	from = 0; // current recombination starting point
      bool mutations = n_mutat > 0;
      
      SEXP mutat = mutations ? VECTOR_ELT(Indiv, MUTAT1 + hapnr) : NULL;
      ToIntLocal( mutat);
      
      for (Uint index=0; index <= n_recomb; index++) {
	
	double cur_recomb = recombi[index];
	Uint cur_chrom; // starting with number 0

	if (cur_recomb <= first_snp_position) {
	  from = 0;
	  continue;
	}
	
	for (cur_chrom=0; cur_chrom<n_chrom; cur_chrom++)
	  if (length_total[cur_chrom + 1] >= cur_recomb) break;
	if (equidist) {
	  if (cur_recomb > last_snp_position) till = n_snps - 1;
	  else {
	    double halfbetween = REAL(VECTOR_ELT(position, cur_chrom))[0];
	    till = totalsnps_prec_chroms[cur_chrom] +
	      CEIL((cur_recomb - length_prec_chroms[cur_chrom] - halfbetween) /
		   (2.0 * halfbetween)) - 1;
	  }
	} else {
	  Uint min = cur_chrom == 0 ? 0 : totalsnps_prec_chroms[cur_chrom],
	    max = cur_chrom == n_chrom ? n_chrom - 1
	    : totalsnps_prec_chroms[cur_chrom + 1] - 1;
	  till = getNrSNPposition(cur_recomb, snp_positions, min, max);
	}

	if (index == 0) {
	  from = till + 1;
	  assert(from == 0);
	  continue;
	}
	
	assert(toSNP != NA_INTEGER && toSNP <= MAXINT);
	if (till >= from && from <= (Uint) toSNP && till >= (Uint) fromSNP) {
	  Uint von = MAX((Uint) fromSNP, from),
	    bis = MIN((Uint) toSNP, till),
	    o_info[4];
 	  decodeOrigin(originint[index-1], o_info);
	  //
	  assert(bis < n_snps);
	  assert(von <= bis);

	  //  SEXP Orig_gen = VECTOR_ELT(breeding, o_info[ORIGIN_GENERATION]), 
	  SEXP Orig_gen = VECTOR_ELT(breeding,
				     gen_ref[o_info[ORIGIN_GENERATION]]), 
	    O2 = VECTOR_ELT(Orig_gen, o_info[ORIGIN_SEX]),
	    Orig = VECTOR_ELT(O2, o_info[ORIGIN_NR]);
	  Uint o_haplotype = o_info[ORIGIN_HAPLO];
	  SEXP coded_haplo = VECTOR_ELT(Orig, CODEDHAPLO);
	  Uint *o_haplo = (Uint*) algn_generalS(coded_haplo, MY_LDABITALIGN);
	 
	  for (Uint op=von; op<=bis; ) {
	    if ( (op % CodesPerUnit) == (q % CodesPerUnit) && !do_select) {
	      // unit-wise copying	      
	      Uint 
		bnr = q / CodesPerUnit, 
		start = q % CodesPerUnit,
		o_rest = 1 + bis - op,
		rest = CodesPerUnit - start;
	      if (rest > o_rest) rest = o_rest;
	      Uint o_code = o_haplo[op / CodesPerUnit];
	      
	      if (mutations) {
		Uint end = op + rest; // excluded, R-Coding !!
		if (mutatint[pmutat] <= end) { // mutat hat R codierung !!
		  // insbesondere bei select koennen Mutation uebersprungen
		  // werden. Somit muss zunaechst pmutat gegebenfalls erhoeht
		  // werden:
		  while (pmutat < n_mutat && mutatint[pmutat] <= op) pmutat++;
		  while (pmutat < n_mutat && mutatint[pmutat] <= end) {
		    Uint idx = (mutatint[pmutat] - 1) % CodesPerUnit;
		    //  AND(L1 , BitMaskStart[idx], BitMaskEnd[idx]);
		    //		    XOR(o_code, o_code, L1);
		    o_code ^= BitMaskStart[idx] & BitMaskEnd[idx];
		    pmutat++;
		  }
		  mutations = pmutat < n_mutat;
		}
	      }
	      //   AND(L1, BitMaskStart[start], BitMaskEnd[start + rest - 1]);
	      //   AND(L1, L1, NullEins[o_haplotype]);
	      // AND(o_code, o_code, L1);
	      o_code &= BitMaskStart[start] & BitMaskEnd[start + rest - 1] &
		NullEins[o_haplotype];

	      if (o_haplotype != hapnr) {		
		// if (o_haplotype > hapnr) { SHR32(o_code, o_code, 1); }
		// else SHL32(o_code, o_code, 1);
		if (o_haplotype > hapnr) o_code >>= 1; else o_code <<= 1;
	      }
	      
	      //OR(haplo[bnr], haplo[bnr], o_code);
	      haplo[bnr] |= o_code;
	      op += rest;
	      p += rest; // p is now endpoint + 1 !
	      q += rest;
	    } else {
	      BUG;
	      // bitwise copying
	      if (do_select) {
		if (psel >= len_sel || cur_select != p) {
		  op++;
		  p++;
		  continue;
		}
		cur_select = Selectint[psel++];
	      }
	      Uint H = GetTwoBitHaploCode(o_haplo, op, o_haplotype);
	      if (mutations) {
		// note mutat:R coding; op:C-Coding
		while (pmutat < n_mutat && mutatint[pmutat] <= op) pmutat++;
		if (mutatint[pmutat] == op + 1) {
		  H = 1- H;
		  pmutat++;
		}
		mutations = pmutat < n_mutat;
	      }
	      if (H) PutTwoBitHaploOne(q, hapnr, haplo);
	      op++;
	      p++;
	      q++;
	    }
	  }// for loop through all snps between two recombination points
	} // if space between two recombination points is not empty
	from = till + 1;
      } // for n_recomb
      FREElocal(mutat);
      FREElocal(origin);
    } // hapnr 
  }

  FREElocal(snps_perchrom);
  FREElocal(Select);
  FREElocal(Generation);
  FREElocal(Sex); 
  FREElocal(Nr);
  if (n_protect) UNPROTECT(n_protect);
   
  return NULL;
}

    
SEXP computeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		 SEXP From_SNP, SEXP To_SNP, SEXP Select, SEXP Geno){
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
  assert_MoBPS();
  bool geno = LOGICAL(Geno)[0];
  Long individuals = LENGTH(Generation);
  SEXP Popinfo;
  PROTECT(Popinfo = IcomputeSNPS(population, Generation, Sex, Nr,
			   INTEGER(From_SNP)[0], INTEGER(To_SNP)[0],
				 Select, KT, NULL));
  Uint *popinfo = (Uint*) INTEGER(Popinfo),
    snps = popinfo[CURRENT_SNPS];
  SEXP Ans;
   PROTECT(Ans = createSNPmatrix(snps, individuals, 
				 geno ? TwoBitGeno : TwoBitHaplo, 0,
				 MY_LDABITALIGN, global, utils ));
   //  printf("mopbs hier\n");
  
  unit_t *ans = Align(Ans, opt); 
  ///printf("mopbs hierX\n");
  IcomputeSNPS(population, Generation, Sex, Nr, INTEGER(From_SNP)[0],
	       INTEGER(To_SNP)[0], Select, KT, ans);
  //printf("mopbs hierZX\n");
  if (geno) haplo2geno(Ans, opt);

  //printf("mopbs hier done\n");
  UNPROTECT(2);
  return Ans;
 }



SEXP compute(SEXP popul, SEXP Generation, SEXP Sex, SEXP Nr,
	     SEXP Tau, SEXP Vec, SEXP Betahat,
	     SEXP Select, SEXP Matrix_return) {
  coding_type coding = assert_MoBPS(); 
  KEY_type *KT = KEYT_M();
  option_type *global = &(KT->global);
  utilsoption_type *utils = &(KT->global_utils);
  basic_options *opt = &(utils->basic);
  int cores = GreaterZero(KT->global_utils.basic.cores);
  int
    variant = check_variant(coding, KT->global.tuning.variant, opt->efficient);
  SEXP Popinfo;
  PROTECT(Popinfo = IcomputeSNPS(popul, Generation, Sex, Nr, 1, NA_INTEGER,
				 Select, KT, NULL));

  Uint *popinfo = (Uint*) INTEGER(Popinfo),
    snps = popinfo[CURRENT_SNPS],
    individuals = LENGTH(Generation);
  Long
    ldaHPL = Lda2Bit(snps, MY_LDABITALIGN),
    memInUnits = (Long) individuals * ldaHPL,
    mem = calculateAlignedMem(memInUnits, TwoBitHaplo);
  int *G = (int *) CALLOC(mem, BytesPerUnit);

  unit_t *g = (unit_t*) algn_generalL(G, MY_LDABITALIGN);
  if ((Uint) LENGTH(Vec) != individuals) ERR0("'vec' not of correct length");

  IcomputeSNPS(popul, Generation, Sex, Nr, 1, NA_INTEGER, // R coded indices !!
	       Select, KT, g);
  haplo2geno(g, snps, individuals, ldaHPL, TwoBitHaplo, variant, cores);
	   
  double *A = crossprod((unit_t*) G, // no PROTECT( needed;
			snps, individuals, SnpLDA(snps, coding, MY_LDABITALIGN),
			coding, variant,		
			RowMeans, NoNormalizing, false,
			NULL,
			opt, &(global->tuning));
  
  FREE(G);

  SEXP Ans = IsolveRelMat(individuals, A, // no PROTECT( needed;
			  REAL(Tau), LENGTH(Tau),
			  REAL(Vec),
			  REAL(Betahat), LENGTH(Betahat),
			  2 + (int) LOGICAL(Matrix_return)[0],
			  true, cores);
  FREE(A);
  UNPROTECT(1);
  return Ans;
}

#endif
