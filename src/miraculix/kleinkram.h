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



// by 3.2.2021: xAx:: BLAS lohnt noch nicht
//              A^t A: BLAS lohnt sich ab aA = k x n, k >=8, n > MAXOWN

#ifndef kleinkram_rfutils_h
#define kleinkram_rfutils_h 1

#define SCALAR_VERSION 1 // alternatively 8


typedef char name_type[][MAXCHAR];

void strcopyN(char *dest, const char *src, int n);


#if defined compatibility_to_R_h
#define INT Integer(el, name, 0)
#define LOGI Logical(el, name, 0)
#define NUM Real(el, name, 0)
#define USRLOG UsrBool(el, name, 0)
#define USRLOGRELAXED UsrBoolRelaxed(el, name, 0)
#define CHR Char(el, name)
#define STR(X, N)  strcopyN(X, CHAR(STRING_ELT(el, 0)), N);
#define POS0INT NonNegInteger(el, name) /* better: non-negative */
#define POS0NUM NonNegReal(el, name)
#define NEG0NUM NonPosReal(el, name)
#define POSINT PositiveInteger(el, name) /* better: non-negative */
#define POSNUM PositiveReal(el, name)

usr_bool UsrBool(SEXP p, char *name, Long idx);
usr_bool UsrBoolRelaxed(SEXP p, char *name, Long idx);
SEXP Logic(bool* V, Long n, Long max) ;
SEXP Char(const char **V, Long n, Long max) ;
SEXP String(char V[][MAXCHAR], Long n, Long max);
SEXP Mat(double* V, Long row, Long col, Long max);
SEXP Mat_t(double* V, Long row, Long col, Long max);
SEXP MatInt(int* V, Long row, Long col, Long max) ;
SEXP MatString(char **V, Long row, Long col, Long max);
//SEXP Array3D(int** V, Long depth, Long row, Long col, Long max) ;

SEXP Logic(bool* V, Long n) ;
SEXP Char(const char **V, Long n) ;
SEXP Mat(double* V, Long row, Long col);
SEXP Mat_t(double* V, Long row, Long col);
SEXP MatInt(int* V, Long row, Long col) ;
SEXP MatString(char** V, Long row, Long col);
//SEXP Array3D(int** V, Long depth, Long row, Long col) ;
SEXP String(char *V);
SEXP String(char V[][MAXCHAR], Long n);
SEXP String(int *V, const char * List[], Long n, Long endvalue);

SEXP Num(double* V, Long n, Long max) ;
SEXP Num(double* V, Long n) ;
SEXP Int(int *V, Long n, Long max) ;
SEXP Int(int *V, Long n);
bool Logical(SEXP p, char *name, Long idx);

double Real(SEXP p, char *name, Long idx);
void Real(SEXP el,  char *name, double *vec, Long maxn) ;

int Integer(SEXP p, char *name, Long idx, bool nulltoNA) ;
int Integer(SEXP p, char *name, Long idx);
void Integer(SEXP el, char *name, int *vec, Long maxn) ;
void Integer2(SEXP el, char *name, int *vec) ;
char Char(SEXP el, char *name) ;
int NonNegInteger(SEXP el, char *name) ;
double NonNegReal(SEXP el, char *name) ;
double NonPosReal(SEXP el, char *name) ;
int PositiveInteger(SEXP el, char *name) ;
double PositiveReal(SEXP el, char *name) ;
void String(SEXP el, char *name, char names[][MAXCHAR], Long maxlen);


SEXP TooLarge(int *n, Long l);
SEXP TooSmall();

SEXP ExtendedInteger(double x);
SEXP ExtendedBooleanUsr(usr_bool x);
 

void GetName(SEXP el, char *name, const char * List[], int n,
	     int defaultvalue, int endvalue, int *ans, int maxlen_ans);
int GetName(SEXP el, char *name, const char * List[], int n) ;
int GetName(SEXP el, char *name, const char * List[], int n,
	    int defaultvalue) ;
#endif


double XkCXtl(double *X, double *C, Long nrow, Long dim, Long k, Long l,
	      int cores);
void XCXt(double *X, double *C, double *V, Long nrow, Long dim, int cores);
void AtA(double *a, Long nrow, Long ncol, double *A, int cores);
void AtA(double *a, Long nrow, Long ncol, double *A, int cores, int mode);
void AtAInt(int *a, int *b,
       Long nrow, Long ncol, // a and b have same size
       Long ld, Long ldR,
       Long *C, // result
       int cores,
       Long m, Long n);
void AtAIntBlock3x3(int *a, int *b,
	    Long nrow, Long ncol, // a and b have same size
	    Long ld, Long ldR,
	    Long *C, // result
	    int cores, int scalarVersion,
	    Long m, Long n);
void AtAIntSpeedBlock3x3(int *a, int *b,
	    Long nrow, Long ncol, // a and b have same size
	    Long ld, Long ldR,
	    Long *C, // result
	    int cores,
	    Long m, Long n);
void AtAIntBlock2x2(int *a, int *b,
	    Long nrow, Long ncol, // a and b have same size
	    Long ld, Long ldR,
	    Long *C, // result
	    int cores, int scalarVersion,
	    Long m, Long n);
void AtAInt(int *a, int *b,
       Long nrow, Long ncol, // a and b have same size
       Long ld, Long ldR,
       Long *C, // result
	    int cores, int mode);
void xA(double *x, double*A, Long nrow, Long ncol, double *y, int cores);
void xA_noomp(double *x, double*A, Long nrow, Long ncol, double *y);
void xA(double *x1, double *x2,  double*A, Long nrow, Long ncol, double *y1,
	double *y2);
double xAx(double *x, double*A, Long nrow, int cores);
void Ax(double *A, double*x, Long nrow, Long ncol, double *y, int cores);
void Ax(double *A, double*x1, double*x2, Long nrow, Long ncol, double *y1,
	double *y2);
double xUy(double *x, double *U, double *y, Long dim, int cores);
double xUxz(double *x, double *U, Long dim, double *z, int cores);
double x_UxPz(double *x, double *U, double *z, Long dim, int cores);
double xUx(double *x, double *U, Long dim, int cores);
void matmult(double *A, double *B, double *C, Long l, Long m, Long n,
	     int cores);
void matmulttransposed(double *A, double *B, double *C, Long m, Long l, Long n,
		       int cores);
//void matmulttransposedInt(Long *A, Long *B, Long *c, Long m, Long l, Long n); 
void matmult_2ndtransp(double *A, double *B, double *C, Long m, Long l, Long n,
		       int cores);
void matmult_2ndtransp(double *A, double *B, double *C, Long m, Long l,
		       int cores);
void matmult_tt(double *A, double *B, double *C, Long m, Long l, Long n,
		int cores);
double *matrixmult(double *m1, double *m2, Long dim1, Long dim2, Long dim3,
		   int cores);


#define MULTIPLEMATCHING -2
#define NOMATCHING -1
#define MATCHESINTERNAL -3
#define NOSEPMATCH -4
int Match(char *name, name_type List, int n);
int Match(char *name, const char * List[], int n);
int Match(char *name, const char * List[], int n, char sep, int *Len);
int Match(char *varname, char* name, const char * List[], int n,
	  char sep, int *Len, bool relax);

#define TYPE_INDEP_SCALAR_PROD(A, B, N, ANS) {			\
    Long  k_ =0,				\
    end_ = N - 4;				\
  ANS = 0.0;					\
  for (; k_<end_; k_+=4) {				\
    ANS += A[k_] * B[k_]				\
      + A[k_ + 1] * B[k_ + 1]				\
      + A[k_ + 2] * B[k_ + 2]				\
      + A[k_ + 3] * B[k_ + 3];				\
  }							\
  for (; k_<N; k_++) ANS += A[k_] * B[k_];		\
  }


// unused :
#define FILL_IN(A, N, VALUE) {				\
    Long end_ = N;					\
    for (Long k_=0; k_<end_; (A)[k_++]=VALUE);		\
}

double ownround(double x);

#define Mod(Z, modulus) ((Z) - FLOOR((Z) / (modulus)) * (modulus))
double lonmod(double x, double modulus);

void Minus_int(int *A, Long nrowA, Long ncolA, Long storingNrowA,
	       int *B, Long nrowB, Long ncolB, Long storingNrowB,
	       Long storingNrowR, int *result);

void matchError(char *varname, char* name, int intvalue, 
		const char * List[], int n, int errtype);
#endif
