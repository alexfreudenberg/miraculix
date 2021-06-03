
#ifndef miraculix_Haplo_H
#define miraculix_Haplo_H 1

//SEXP codeHaplo(SEXP M1, SEXP M2);
Uint BytesPerBlockHaplo();
Uint CodesPerBlockHaplo();
Uint BitsPerCodeHaplo();
void InitGetHaplo(SEXP CM, Uint **code, Uint *units);
Uint GetHaplo(Uint *cm, Uint snp);
//Uint GetHaploX(Uint *cm, Uint snp);
SEXP codeHaplo2(SEXP M1, SEXP M2);
SEXP create_codevectorHaplo(Uint snps, Uint individuals);
SEXP create_codevectorHaplo(Uint snps, Uint individuals, SEXP);
void zeroNthHaplo(SEXP CM, SEXP NN);

typedef Uint (*GetHaploFun)(Uint *cm, Uint snp);

#endif
