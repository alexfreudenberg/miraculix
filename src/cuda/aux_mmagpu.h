
/*
    This file replaces the R-specific types and functions in mmagpu.cu to allow the utilities to be used outside of R.

*/
#ifndef AUX_MMAGPU 
#define AUX_MMAGPU 1

#define CodesPerByte 4L
#define BytesPerUnit 4L
#define Uint unsigned int
#define Ulong unsigned long

#define SEXP int
void PROTECT(SEXP a){}
void UNPROTECT(SEXP a){}

#define PRINTF printf
#define MEMCOPY memcpy
void ERR(const char* s){
    printf("%s",s);
    exit(1);
}
long static inline UnitsPerIndiv256(Ulong snps){
    return (1L + (snps - 1L)/((4L * 8L / 2) * (32L / 4L))) * (32L / 4L);
}  

#endif