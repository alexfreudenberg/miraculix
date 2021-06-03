
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
#define PlainInteger256

#include "intrinsics.h"
#include "miraculix.h"
#include <General_utils.h>

#include <inttypes.h> 
#include "Vector.matrix.h"
#include "AutoMiraculix.h"
#include "options.h"
#include "Haplo.h"
#include "haplogeno.h"
#include "xport_import.h"
//#include "shuffle.h"
#include "align.h"
#include "utils.h"

//#define DO_FLOAT 1  // must be completely re-checked if defined

// #define DO_PARALLEL 1 // done by basic_utils.h

// #define INLINE
#define INLINE static inline



//
#define show(X)
// void show(const char *s) {  if (debugging) PRINTF("%.50s", s); }
				   
void checkMethod(snpcoding method) {
  KEY_type *KT = KEYT();
  if (!is256(method)) 
    ERR1("MoBPS does not work with method '%.20s'\n", SNPCODING_NAMES[method]);

  if (method == TwoBit) {
    if (KT->doshow) {
      PRINTF("MoBPS will run much faster and with a minimum of memory if 'Shuffle' or 'Shuffle256' is used as 'snpcoding'\n");
      KT->doshow = false;
    }
  } else if (method != Shuffle256) {
      if (KT->doshow) {
      PRINTF("MoBPS will run twice as fast if 'Shuffle256' is used as 'snpcoding'\n");
     KT-> doshow = false;
    }
  }
  
}
 
BlockType BitMaskStart[CodesPerBlock], BitMaskEnd[CodesPerBlock];
bool MoBPSNotInit = true;
void InitMoBPS() {
  if (!MoBPSNotInit) BUG;
  MoBPSNotInit = false;
  assert(BytesPerBlock == sizeof(BlockType0));
  BlockType Einsen = 0xFFFFFFFF;
  BitMaskEnd[0] = 0x03;

  for (Uint k=1; k<CodesPerBlock; k++) {
    //    XOR(BitMaskStart[k], BitMaskEnd[k-1], Einsen); // not
    BitMaskStart[k] = BitMaskEnd[k-1] ^ Einsen; // not
    BitMaskEnd[k] = BitMaskEnd[k-1];
    BitMaskEnd[k] |= 0x03 << (BitsPerCode * k);
  }
  BitMaskStart[0] = Einsen;
}


void assert_MoBPS() {
  if (MoBPSNotInit) InitMoBPS();
  option_type *global = &(KEYT()->global);
  checkMethod((snpcoding) global->genetics.method);
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


void INLINE decodeOrigin(Uint x, Uint *ans) {
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
  for (Uint i=0; i<4; ans[i++]++);
  UNPROTECT(1);
  return Ans;
}


SEXP codeOrigins(SEXP M) {
  assert_MoBPS();
  Uint nrow = nrows(M);
  assert(ncols(M) == 4);
  
  SEXP Ans;
  PROTECT(Ans = allocVector(INTSXP, nrow));
 
  Uint *ans = (Uint*) INTEGER(Ans);
  for (Uint i=0; i<nrow; i++) ans[i] = 0;  // unnoetig?!

#define DO1(TYPE, RTYPE)			\
  TYPE *gen = RTYPE(M),				\
    *sex = gen + nrow,				\
    *nr = sex + nrow,				\
    *haplo = nr + nrow					       

#ifdef SCHLATHERS_MACHINE  
#define CONTROL \
    Uint control[4];							\
    decodeOrigin(ans[i], control);					\
    if (control[0]+1 != g || control[1]+1 != s || control[2]+1 != n ||	\
	control[3]+1 != H) {						\
      /*	p rintf("i=%d: generation:%d->%d sex:%d->%d nr:%d->%d haplo:%d->%d\n", i, g, control[0], s, control[1], n, control[2], H, control[3]); */ \
      BUG;								\
    }
#else
  #define CONTROL
#endif  
  
#define DO2								\
  for (Uint i=0; i<nrow; i++) {						\
    Uint g = (Uint) gen[i],						\
      s = (Uint) sex[i],						\
      n = (Uint) nr[i],							\
      H = (Uint) haplo[i];						\
    if (g < 1 || g > MAX_GENE_INPUT|| s < 1 || s > MAX_SEX ||		\
	n < 1 || n > MAX_INDIVIDUALS ||	H < 1 || H > MAX_HAPLO)		\
      ERR1("some value in row %d out of bound", i + 1);			\
    ans[i] = ((((((g - 1) << BITS_SEX) +				\
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
#define SNPS_EQUIDIST 17
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


Uint INLINE GetHaploCode(BlockType0 *code, Uint pos, Uint haplo) {
  Uint elmt = pos % CodesPerBlock; 
  assert(haplo <= 1);
  assert(elmt < CodesPerUnit);
  return (code[pos / CodesPerBlock] >> (2 * elmt + haplo)) & 0x00000001;
}


/*
void INLINE PutHaploCode(Uint value, Uint pos, Uint haplo, BlockType0 *code) {
  Uint elmt = pos % CodesPerBlock;
  BlockType shift = 0x00000001 << (2 * elmt + haplo);
  BlockType *c = code + pos / CodesPerBlock; 
  if (value) *c |= shift;
  else *c &= 0xFFFFFFFF xor shift;
}
*/


void INLINE PutHaploOne(Uint pos, Uint haplo, BlockType0 *code) {
  Uint elmt = pos % CodesPerBlock;
  BlockType *c = code + pos / CodesPerBlock; 
  *c |= 0x00000001 << (2 * elmt + haplo);
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
  assert(STRCMP(CHAR(STRING_ELT(getAttrib(WHERE, R_NamesSymbol), WHAT)), \
		NAME) == 0);
SEXP IcomputeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
	       Rint From_SNP, Rint To_SNP, SEXP Select, Uint *ans){
  
  Uint
    cur_select = NA_INTEGER,
    individuals = length(Generation),
    len_sel = length(Select);
  bool do_select = len_sel > 0;
  if (individuals != (Uint) length(Sex) || individuals != (Uint) length(Nr))
    ERR("'gen', 'gender' and 'nr' do not have same length");
  if (From_SNP <= 0) ERR("'from_p' not valid.");

  CondToInt(do_select, Select);
   
  ASSERT(population, BREEDING, "breeding");
  SEXP breeding = VECTOR_ELT(population, BREEDING);

  ASSERT(population, INFO, "info")
    SEXP info = VECTOR_ELT(population, INFO);
  Uint len_info = length(info);

  static Uint gen_ref[MAX_GENE_INPUT] = { 0 };
  if (GENERATION_REFTABLE >= len_info ||
      length(VECTOR_ELT(info, GENERATION_REFTABLE)) == 0) {
    if (gen_ref[MAX_GENE_INPUT - 1] == 0) { // should never happen except uninit
      for (Uint i=0; i<MAX_GENE_INPUT; i++) gen_ref[i] = i;
      PRINTF("Generation reference table now initialized.\n"); 
    }
  } else {
    ASSERT(info, GENERATION_REFTABLE, "origin.gen");
    Uint *p = (Uint*) INTEGER(VECTOR_ELT(info, GENERATION_REFTABLE)); 
    for (Uint i=0; i<MAX_GENE_INPUT; i++) gen_ref[i] = (Uint) (p[i] - 1);
  }

  ASSERT(info, LENGTH_TOTAL, "length.total");
  Uint n_chrom = length(VECTOR_ELT(info, LENGTH_TOTAL)) - 1;
  if (n_chrom >= MAX_CHROM) ERR1("maximum chromosoms is %d", MAX_CHROM);
  double *length_total = REAL(VECTOR_ELT(info, LENGTH_TOTAL));
  
  ASSERT(info, LENGTH_CHROM, "length");
  double *LengthChrom = REAL(VECTOR_ELT(info, LENGTH_CHROM));
  
  ASSERT(info, NSNPS, "snp");
  Uint n_snps;
  SEXP InfoHaplo = VECTOR_ELT(info, INFO_INFOHAPLO);
  option_type *global = &(KEYT()->global);

  if (TYPEOF(InfoHaplo) != INTSXP || length(InfoHaplo) != INFO_LAST + 1) {
    PROTECT(InfoHaplo = allocVector(INTSXP, INFO_LAST + 1));
    int *infoHaplo = INTEGER(InfoHaplo);
    for (int i=0; i<=INFO_LAST; infoHaplo[i++] = NA_INTEGER);
    infoHaplo[METHOD] = Haplo;
    snpcoding method = check_method((snpcoding) global->genetics.method);
    if (!is256(method)) ERR1("'%.20s'not allowed", SNPCODING_NAMES[1 + method] );
    SET_VECTOR_ELT(info, INFO_INFOHAPLO, InfoHaplo);
    UNPROTECT(1);    
  }
   
  double length_prec_chroms[MAX_CHROM + 1];
  SEXP snps_perchrom = VECTOR_ELT(info, NSNPS);
  ToInt(snps_perchrom);
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
  Rint fromSNP = From_SNP - 1,
    toSNP = To_SNP==NA_INTEGER || To_SNP >= (int) n_snps ? MAXINT : To_SNP -1;
  if (toSNP < fromSNP) ERR("'to_p' smaller than 'from_p'");
  if (fromSNP % CodesPerBlock != 0)
    ERR1("(from_p-1) must be devisible by %d", CodesPerBlock);
  if (toSNP < MAXINT && (toSNP + 1)  % CodesPerBlock != 0)
    ERR1("'to_p' must be devisible by %d or 'NA'", CodesPerBlock);


  Uint cur_nsnps = 0;
  if (do_select) {
    cur_select = Selectint[0] - 1;
    Rint toSNP_P1 = toSNP + 1;
    for (Uint k=0; k<len_sel; k++)
      cur_nsnps +=
	Selectint[k] > (Uint) fromSNP && Selectint[k] <= (Uint) toSNP_P1;
  } else cur_nsnps = 1L + MIN(n_snps - 1L, (Uint) toSNP) - fromSNP;
  INTEGER(InfoHaplo)[CURRENT_SNPS] = cur_nsnps;
  
  if (ans == NULL) return InfoHaplo;
  Uint unitsPerIndiv = GetUPI(cur_nsnps, Haplo);

  ASSERT(info, SNPS_EQUIDIST, "snps.equidistant");
  bool equidist = (bool) LOGICAL(VECTOR_ELT(info, SNPS_EQUIDIST))[0];
  
  ASSERT(info, POSITION, "position");
  SEXP position = VECTOR_ELT(info, POSITION);

  ASSERT(info, SNP_POSITION, "snp.position");
  double *snp_positions = REAL(VECTOR_ELT(info, SNP_POSITION)),
    last_snp_position = snp_positions[n_snps - 1],
    first_snp_position = snp_positions[0];

  //  printf("%d %d %d\n", length(Generation), length(Sex), length(Nr));
  ToInt(Generation); 
  ToInt(Sex);
  ToInt(Nr);
  //  for (int i=0; i<length(Generation); i++)
  //	 printf("%d %d %d\n", Generationint[i], Sexint[i], Nrint[i]);

  BlockType NullEins[2] = { 0x55555555, 0xAAAAAAAA};

  for (Ulong i=0; i<individuals; i++) {
    //    print("i=%d %d\n", i, individuals);
    if (length(breeding) < (int) Generationint[i])
      ERR("value of 'generation' too large.");
    SEXP
      g = VECTOR_ELT(breeding, Generationint[i]-1);
    if (length(g) < (int) Sexint[i]) ERR("'sex' out of range.");
    SEXP 
      s = VECTOR_ELT(g, (int) Sexint[i]-1);
    if (length(s) < (int) Nrint[i]) ERR("value of 'nr' too large.");  
    SEXP
      Indiv = VECTOR_ELT(s, Nrint[i]-1);
    BlockType *haplo = (BlockType0 *) (ans + i * unitsPerIndiv);
    for (Uint hapnr=0; hapnr<=1; hapnr++) {
      if (TYPEOF(VECTOR_ELT(Indiv, RECOMBI1 + hapnr)) != REALSXP)
	ERR("real valued vector for recombination information expected");
      double *recombi = REAL(VECTOR_ELT(Indiv, RECOMBI1 + hapnr));
      Uint
	n_recomb = length(VECTOR_ELT(Indiv, RECOMBI1 + hapnr)) - 1;
      SEXP origin = VECTOR_ELT(Indiv, ORIGIN1 + hapnr);
      ToInt(origin);
      Uint
	n_mutat = length(VECTOR_ELT(Indiv, MUTAT1 + hapnr)),
	// separation of the following 3 allows for selection but also
	// for mutation error such as insertion and deletion
	p = 0, // position what is written if there is no selection
	q = 0, // position what is written actually
	psel = 0, pmutat = 0,
	till = 0,
	from = 0; // current recombination starting point
      bool mutations = n_mutat > 0;
      
      SEXP mutat = mutations ? VECTOR_ELT(Indiv, MUTAT1 + hapnr) : NULL;
      CondToInt(mutations, mutat);
      
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
	    o_info[4],
	    n_protect = 0;
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
	  //	  if (n_snps <= CodesPerBlock*3)//Coding selbst + hinten & vorne Platz
	  //   ERR2("number of snps is %d, but should be more than %d so that the programme can be used.", n_snps, CodesPerBlock * 3);
	  
	  //  printf("%d %d\n", length(coded_haplo), n_snps);
	  //if ((Uint) length(coded_haplo) == n_snps) { // unsinn
	  //	    printf("L %d %d %d H%d %d OG %d %d S%d\n", length(coded_haplo), n_snps,  o_info[ORIGIN_NR], HAPLOTYP1, HAPLOTYP2, o_info[ORIGIN_GENERATION], gen_ref[o_info[ORIGIN_GENERATION]], o_info[ORIGIN_SEX]);
	  //	    PROTECT(coded_haplo = codeHaplo2(VECTOR_ELT(Orig, HAPLOTYP1),//
	  //					     VECTOR_ELT(Orig, HAPLOTYP2)));
	  //SET_VECTOR_ELT(Orig, CODEDHAPLO, coded_haplo);
	  //n_protect++;
	  //}
 
	  BlockType
	    *o_haplo = (BlockType0*) Align256(coded_haplo, ALIGN_CODED);
	 
	  for (Uint op=von; op<=bis; ) {
	    if ( (op % CodesPerBlock) == (q % CodesPerBlock) && !do_select) {
	      // blockwise copying	      
	      Uint 
		bnr = q / CodesPerBlock, 
		start = q % CodesPerBlock,
		o_rest = 1 + bis - op,
		rest = CodesPerBlock - start;
	      if (rest > o_rest) rest = o_rest;
	      BlockType o_code = o_haplo[op / CodesPerBlock];
	      
	      if (mutations) {
		Uint end = op + rest; // excluded, R-Coding !!
		if (mutatint[pmutat] <= end) { // mutat hat R codierung !!
		  // insbesondere bei select koennen Mutation uebersprungen
		  // werden. Somit muss zunaechst pmutat gegebenfalls erhoeht
		  // werden:
		  while (pmutat < n_mutat && mutatint[pmutat] <= op) pmutat++;
		  while (pmutat < n_mutat && mutatint[pmutat] <= end) {
		    Uint idx = (mutatint[pmutat] - 1) % CodesPerBlock;
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
	      Uint H = GetHaploCode((BlockType0 *) o_haplo, op, o_haplotype);
	      if (mutations) {
		// note mutat:R coding; op:C-Coding
		while (pmutat < n_mutat && mutatint[pmutat] <= op) pmutat++;
		if (mutatint[pmutat] == op + 1) {
		  H = 1- H;
		  pmutat++;
		}
		mutations = pmutat < n_mutat;
	      }
	      if (H) PutHaploOne(q, hapnr, haplo);
	      op++;
	      p++;
	      q++;
	    }
	  }// for loop through all snps between two recombination points
	  if (n_protect > 0) UNPROTECT(n_protect);
	} // if space between two recombination points is not empty
	from = till + 1;
      } // for n_recomb
      FREEint(mutat);
      FREEint(origin);
    } // hapnr 
  }

  FREEint(snps_perchrom);
  FREEint(Select);
  FREEint(Generation);
  FREEint(Sex); 
  FREEint(Nr); 

  
  return NULL;
}

    
SEXP computeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		 SEXP From_SNP, SEXP To_SNP, SEXP Select, SEXP Geno){
  assert_MoBPS();
  bool geno = LOGICAL(Geno)[0];
  Uint individuals = length(Generation);
  SEXP Popinfo;
  PROTECT(Popinfo = IcomputeSNPS(population, Generation, Sex, Nr,
			   INTEGER(From_SNP)[0], INTEGER(To_SNP)[0],
				 Select, NULL));
  Uint *popinfo = (Uint*) INTEGER(Popinfo),
    snps = popinfo[CURRENT_SNPS];
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  SEXP Ans;
  PROTECT(Ans = createSNPmatrix(snps, individuals, Haplo));
  
  Uint *ans =  Align256(Ans, ALIGN_COMPUTE); 
  IcomputeSNPS(population, Generation, Sex, Nr, INTEGER(From_SNP)[0],
	       INTEGER(To_SNP)[0], Select, ans);
  if (geno) {
    if (!is256(method)) ERR("unallowed snp coding method found");
      assert(ans == Align256(Ans, ALIGN_COMPUTE));
      haplo2geno(Ans, ans, method);
  }

  UNPROTECT(2);
  return Ans;
 }



SEXP compute(SEXP popul, SEXP Generation, SEXP Sex, SEXP Nr,
	     SEXP Tau, SEXP Vec, SEXP Betahat,
	     SEXP Select, SEXP Matrix_return) {
  assert_MoBPS(); // asserts(is256(method));
  option_type *global = &(KEYT()->global);
  snpcoding method = check_method((snpcoding) global->genetics.method);
  SEXP Popinfo;
  PROTECT(Popinfo = IcomputeSNPS(popul, Generation, Sex, Nr, 1, NA_INTEGER,
				 Select, NULL));

  Uint *popinfo = (Uint*) INTEGER(Popinfo),
    snps = popinfo[CURRENT_SNPS],
    individuals = length(Generation),
    unitsPerIndivHPL = UnitsPerIndiv256(snps);
  Ulong
    memInUnits = (Ulong) individuals * unitsPerIndivHPL,
    mem = calculateAlignedMem(memInUnits, Haplo, 0);
  Uint *G = (Uint *) CALLOC(mem, BytesPerUnit);
  if (G == NULL) ERR1("mem allocation (%lu bytes)", mem * BytesPerUnit);

  Uint *g = (Uint*) algn_general(G, BytesPerBlock256);
  if ((Uint) length(Vec) != individuals) ERR("'vec' not of correct length");

  IcomputeSNPS(popul, Generation, Sex, Nr, 1, NA_INTEGER, // R coded indices !!
	       Select, g);
  haplo2geno(g, snps, individuals, method, unitsPerIndivHPL, g);

  
  double *A = crossprod(G, // no PROTECT( needed;
			  snps, individuals, method, true, true, 0);
  
  FREE(G);

  SEXP Ans = IsolveRelMat(individuals, A, REAL(Tau)[0], // no PROTECT( needed;
			  REAL(Vec),
			  REAL(Betahat)[0],
			  2 + (int) LOGICAL(Matrix_return)[0],
			  true);
  FREE(A);
  UNPROTECT(1);
  return Ans;
}

