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

#define RFU_LOCAL 1  

/*

mssngs auf Long, Long

 */

#include "intrinsics.h"
#if ! defined SCHLATHERS_MACHINE
//#error OFF
#endif
#include "def.h" // just to define STAND_ALONE
#if defined SCHLATHERS_MACHINE
//#error ON
#endif
#include "Basic_RandomFieldsUtils.h"
#include "compatibility.lapack.h"
#include "compatibility.C.h"
#include "kleinkram.h"
#include "zzz_RFU.h"
#include "xport_import.h"



// R FU INCLUDE #include  "xport_import_RFU.h"


AVAILABLE_SIMD


#define USE_OWN_ALG(SCALAR_LEN, PARALLEL) true
#define USE_OWN_SCALAR_PROD true

#define SCALAR(A,B,C) Ext_scalarX(A,B,C, SCALAR_AVX)

void strcopyN(char *dest, const char *src, int n) {
  if (n > 1) {
    n--; 
    strncpy(dest, src, n);
  }
  dest[n] = '\0';
}


void AtAInt(int *a, int *b,
       Long nrow, Long ncol, // a and b have same size
       Long ld, Long ldC,
       Long *C, // result
       int VARIABLE_IS_NOT_USED cores,
       Long m, Long n) { // tile sizes
 // C =  A^T %*% B
  Long
    bytes = sizeof(*C) * ncol,
    tileRows = m <= 0 ? nrow : m,
    tileCols = n <= 0 ? ncol : n;

  if (ldC == ncol) MEMSET(C, 0,  bytes * ncol);
  else for (Long i=0; i<ncol; i++) MEMSET(C + ldC * i, 0, bytes);

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif
  for (Long tR=0; tR<ncol; tR += tileCols) {
    Long Rend = MIN(tR + tileCols, ncol);
    for (Long tC=tR; tC<ncol; tC += tileCols) {
      Long Cend = MIN(tC + tileCols, ncol);
      for (Long rowPos = 0; rowPos < nrow; rowPos += tileRows) {
	Long curLen = MIN(tileRows, nrow - rowPos);
	for (Long i=tR; i<Rend; i++) {
	  Long j = a==b && tC == tR ? i : tC;
	  int
	    *A = a + i * ld + rowPos,
	    *B = b + j * ld + rowPos;
	  for ( ; j<Cend; j++, B+=ld) {
	    C[i + ldC * j] += Ext_scalarXint(A, B, curLen, SCALAR_VERSION);
	    //  printf("%ld x %ld : %ld\n", i, j, C[i + ldC * j]);
	  }
	}
      }
    }
  }
  if (a==b)
    for (Long i=0; i<ncol; i++) 
      for (Long j=i; j<ncol; j++) 
	C[j + ldC * i] = C[j * ldC + i];
}



    


void AtAInt(int *a, int *b,
       Long nrow, Long ncol, // a and b have same size
       Long ld, Long ldC,
       Long *C, // result
	    int VARIABLE_IS_NOT_USED cores, int mode) { 
 // C =  A^T %*% B
  // printf(" plain cores = %d\n", cores);
  if (mode < 100) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif 
  for (Long i=0; i<ncol; i++) {
    Long j = a==b ? i : 0;
    int
      *A = a + i * ld,
      *B = b + j * ld;
    for ( ; j<ncol; j++, B+=ld) {
      C[i + j * ldC] = Ext_scalarXint(A, B, nrow, mode);
    }
  }
  if (a==b) 
    for (Long i=0; i<ncol; i++) 
      for (Long j=i; j<ncol; j++) 
	C[j + ldC * i] = C[j * ldC + i];
  return;
  }
  
  if (a != b) ERR0("chosen method only works if the matrices are the same.");
  switch(mode) {
    // for research only
  case 101 : {
    Long repet = 4,
      end = nrow - repet + 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif  
    for (Long i=0; i<ncol; i++) {
      int *A = a + i * ld;
      for (Long j=i; j<ncol; j++) {
	int *B = b + j * ld;	 
	Long  tmp = 0,
	k = 0;
	for ( ; k<end; k+=repet) {	  
	  tmp +=  A[k] * B[k]
	    + A[k + 1] * B[k + 1]
	    + A[k + 2] * B[k + 2]
	    + A[k + 3] * B[k + 3];
	}
	for ( ; k<nrow; k++) tmp += A[k] * B[k];
	C[i * ncol + j] = C[i + ncol * j] = tmp;
      }
    }
  }
    break;

  default : {
  // C =  A^T %*% A
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif  
    for (Long i=0; i<ncol; i++) {
      int *A = a + i * ld;
      for (Long j=i; j<ncol; j++) {
        int *B = b + j * ld;
	Long tmp = 0;
	for (Long k=0; k<nrow; k++)  tmp += A[k] * B[k];
	C[i * ncol + j] = C[i + ncol * j] = tmp;
      }
    }
  }
  

  } // switch

  
} // switch


void AtA(double *a, Long nrow, Long ncol, double *C, int VARIABLE_IS_NOT_USED cores, int mode) {
  switch(mode) {
  case 1 :
  // C =  A^T %*% A
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(dynamic, 20) if (MULTIMINSIZE(ncol)) 
#endif  
    for (Long i=0; i<ncol; i++) {
      double 
	*A = a + i * nrow,
	*B = A;
      for (Long j=i; j<ncol; j++, B+=nrow) {
      C[i * ncol + j] = C[i + ncol * j] = SCALAR(A, B, nrow);
      }
    }
    break;
  default :
    double alpha = 1.0,
      beta = 0.0;
    stopIfNotInt(ncol | nrow);   
    int ncol0 = (int) ncol,
      nrow0 = (int) nrow;
    MEMSET(C, 0, ncol * ncol *sizeof(*C));
    F77dsyrk("U","T", &ncol0, &nrow0,
	      &alpha, a,
	      &nrow0,
	      &beta, C,
	      &ncol0
#ifdef USE_FC_LEN_T	  
		  FCONE FCONE
#endif		  
	     );
    for (Long i=1; i<ncol; i++) {
      for (Long j=0; j<i; j++) {
	C[i + ncol * j] = C[i * ncol + j];
      }
    } // for
  } // switch
}
  
void AtA(double *a, Long nrow, Long ncol, double *C, int cores) {
  AtA(a, nrow, ncol, C, cores, USE_OWN_ALG(nrow, ncol) || nrow * ncol > MAXINT);
}


void xA_noomp(double *x, double*A, Long nrow, Long ncol, double *y) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(*y) * nrow);
  } else {
    for (Long i=0; i<ncol; i++) {
      y[i] = SCALAR(x, A + i * nrow, nrow);
    }
  }
}


void xA(double *x, double*A, Long nrow, Long ncol, double *y, int VARIABLE_IS_NOT_USED cores) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(*y) * nrow);
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) if (MULTIMINSIZE(ncol) && MULTIMINSIZE(nrow))
#endif  
   for (Long i=0; i<ncol; i++) y[i] = SCALAR(x, A + i * nrow, nrow);
  } 
} 

  
void xA(double *x1, double *x2,  double*A, Long nrow, Long ncol, double *y1,
	double *y2) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y1, x1, sizeof(*y1) * nrow);
    MEMCOPY(y2, x2, sizeof(*y2) * nrow);
  } else {
    double *a = A;
    for (Long i=0; i<ncol; i++, a += nrow) {
      y1[i] = SCALAR(x1, a, nrow);
      y2[i] = SCALAR(x2, a, nrow);
    }
  }	
}  
  

double xAx(double *x, double*A, Long nrow, int VARIABLE_IS_NOT_USED cores) {
  if (USE_OWN_SCALAR_PROD || nrow * nrow > MAXINT) {
    double sum = 0.0;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) reduction(+:sum) schedule(static) if (MULTIMINSIZE(nrow) && MULTIMINSIZE(nrow))
#endif  
    for (Long i=0; i<nrow; i++)
      sum += x[i] * SCALAR(x, A + i * nrow, nrow);
    return sum;
  } else {
   double alpha = 1.0,
     beta = 0.0;
    int incOne = 1;
    double *y = (double*)  MALLOC(nrow * sizeof(double));
    // z = x^\top A
    stopIfNotInt(nrow);
    int nrow0 = (int) nrow;
    F77dgemv("T", &nrow0, &nrow0, &alpha, A, &nrow0, x, &incOne, &beta,
	     y, &incOne
#ifdef USE_FC_LEN_T	  
		  FCONE
#endif		  
	     );
#if defined compatibility_to_R_h
    // z^top y
    alpha = F77ddot(&nrow0, x, &incOne, y, &incOne);
#else
    alpha = 0.0;
    for (Long i=0; i<nrow; i++) alpha += x[i] * y[i];
#endif
    FREE(y);
    return alpha;
  }
}

void Ax(double *A, double*x, Long nrow, Long ncol, double *y, int VARIABLE_IS_NOT_USED cores) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y, x, sizeof(*y) * nrow);
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) if (MULTIMINSIZE(ncol) && MULTIMINSIZE(nrow))
    for (Long j=0; j<nrow; j++) {
      double tmp = 0.0;
      Long k = j;
      for (Long i=0; i<ncol; i++, k+=nrow) { 
	tmp += A[k] * x[i];
      }
      y[j] = tmp;
    }
#else
    for (Long i=0; i<nrow; i++) y[i]=0.0;
    for (Long k=0, i=0; i<ncol; i++) { 
      for (Long j=0; j<nrow; j++) {
	y[j] += A[k++] * x[i];
      }
    }
#endif  
  }
}


void Ax(double *A, double*x1, double*x2, Long nrow, Long ncol, double *y1,
	double *y2) {
  if (A == NULL) {
    if (nrow != ncol || nrow <= 0) BUG;
    MEMCOPY(y1, x1, sizeof(*y1) * nrow);
    MEMCOPY(y2, x2, sizeof(*y2) * nrow);
  } else {
    for (Long i=0; i<nrow; i++) y1[i]=y2[i]=0.0;
    for (Long k=0, i=0; i<ncol; i++) { 
      for (Long j=0; j<nrow; j++) {
	y1[j] += A[k] * x1[i];
	y2[j] += A[k++] * x2[i];
      }
    }
  }
}


double XkCXtl(double *X, double *C, Long nrow, Long dim, Long k, Long l,
	      int VARIABLE_IS_NOT_USED cores) {
  // (k-th row of X) * C * (l-th row of X)
  // X is nrow x dim matrix
  // C is dim x dim matrix
  double
    *pX = X + k, 
    *pY = X + l, 
    result = 0.0;
  Long size = nrow * dim;
  
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) reduction(+:result)
#endif
  for (Long j=0; j<size; j+=nrow) {
    double scalar = 0.0;
    Long ci = j * dim;
    for (Long i=0; i<size; i+=nrow) scalar += pX[i] * C[ci++];
    result += scalar * pY[j];
  }
  return result;
}


void XCXt(double *X, double *C, double *V, Long nrow, Long dim /* dim of C */, int VARIABLE_IS_NOT_USED cores) {
  Long size = nrow * dim;
  double  
    *endpX = X + nrow,
    *tmp = (double*) MALLOC(sizeof(double) * size); // tmp = XC
 
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores))  schedule(static)
#endif
  for (double *pX = X; pX < endpX; pX++) {
    double *ptmp = tmp + (pX - X);
    for (Long ci=0, cd=0; cd<size; cd+=nrow) {
      double scalar = 0.0;
      for (Long i=0; i<size; i+=nrow) {
        scalar += pX[i] * C[ci++];
      }
      ptmp[cd] = scalar;
    }
  }

  // V = tmp X^t
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores))  schedule(dynamic, 20)
#endif
  for (Long rv=0; rv<nrow; rv++) {
    for (Long cv=rv; cv<nrow; cv++) {
      double scalar=0.0;
      for (Long i=0; i<size; i+=nrow) {
	scalar += tmp[rv + i] * X[cv + i];
     }
      V[rv + cv * nrow] = V[cv + rv * nrow] = scalar;
    }
  }

  UNCONDFREE(tmp);
}


double xUy(double *x, double *U, double *y, Long dim, int VARIABLE_IS_NOT_USED cores) {
  // U a symmetric matrix given by its upper triangular part
  double xVy = 0.0;
  Long    dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) if (MULTIMINSIZE(dim)) reduction(+:xVy) 
#endif  
  for (Long d=0; d<dim; d++) {
    Long i, 
      j = dim * d;
    double tmp = 0.0;
    for (i=0; i<=d; i++) tmp += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) tmp += x[i] * U[j];
    xVy += tmp * y[d];
  }
  return xVy;
}

/*

  // U a symmetric matrix given by its upper triangular part
  assert(z != NULL);
  Long   dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) if (MULTIMINSIZE(dim))
#endif  
  for (Long d=0; d<dim; d++) {
    double tmp;
    Long i,
      j = dim * d;
    for (tmp = 0.0, i=0; i<=d; i++) tmp += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) tmp += x[i] * U[j];
    if (z!=NULL) z[d] = tmp;
  }
  double xVx;
  SCALAR_PROD(z, x, dim, xVx);
  return xVx;

 */

double xUxz(double *x, double *U, Long dim, double *z, int VARIABLE_IS_NOT_USED cores) {
 double xVx = 0.0;
  Long dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) reduction(+:xVx)
#endif
  for (Long d=0; d<dim; d++) {
    Long i, 
      j = dim * d;
    double tmp = 0.0;
    for (tmp = 0.0, i=0; i<=d; i++) tmp += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) tmp += x[i] * U[j];
    if (z != NULL) z[d] = tmp;
    xVx += tmp * x[d];
  }
  return xVx;
}

double xUx(double *x, double *U, Long dim, int VARIABLE_IS_NOT_USED cores) {
  return xUxz(x, U, dim, NULL, cores);
}

double x_UxPz(double *x, double *U, double *z, Long dim, int VARIABLE_IS_NOT_USED cores) {
// x^t (Ux + z); U dreieckmatrix
  double xVx = 0.0;
  Long    dimM1 = dim - 1;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) reduction(+:xVx)
#endif
  for (Long d=0; d<dim; d++) {
    Long i,
      j = dim * d;
    double tmp = z[d];
    for (i=0; i<=d; i++) tmp += x[i] * U[j++];
    for (j += dimM1; i<dim; i++, j+=dim) tmp += x[i] * U[j];
    xVx += tmp * x[d];
  }
  return xVx;
}



void matmult(double *a, double *b, double *c, Long l, Long m, Long n,
	     int VARIABLE_IS_NOT_USED cores) {
// multiplying an lxm- and an mxn-matrix, saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static)
#endif
   for (Long i=0; i<l; i++) {
     double *A = a + i,
       *C = c + i;
     for (Long j=0; j<n; j++) {
       double tmp = 0.0,
	 *B = b + j * m;
       for (Long k=0; k<m; k++) tmp += A[k*l] * B[k];
       C[j * l] = tmp;
     }
   }
}


double *matrixmult(double *m1, double *m2, Long dim1, Long dim2, Long dim3,
		   int VARIABLE_IS_NOT_USED cores) {
  double *m0 = (double*) MALLOC(sizeof(double) * dim1 * dim3);
  matmult(m1, m2, m0, dim1, dim2, dim3, cores);
  return m0;
}


void Xmatmult(double *A, double *B, double *C, Long l, Long m, Long n,
	      int VARIABLE_IS_NOT_USED cores) {
// multiplying an lxm- and an mxn-matrix, saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static)
#endif
  for (Long i=0; i<l; i++) {
    for (Long jl=i, jm=0, j=0; j<n; j++, jl+=l, jm+=m) {
      double tmp = 0.0;
      Long endfor = jm + m;
      for (Long kl=i, k=jm; k<endfor; k++, kl+=l) tmp += A[kl] * B[k]; 
      C[jl] = tmp;
    }
  }
}

void matmulttransposed(double *A, double *B, double *c, Long m, Long l, Long n,
		       int VARIABLE_IS_NOT_USED cores) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static)
#endif
  for (Long i=0; i<l; i++) {    
    double *C = c + i,
      *Aim = A + i * m;
    for (Long j=0; j<n; j++) C[j * l] = SCALAR(Aim, B + j * m, m);
  }
}




void matmult_2ndtransp(double *a, double *B, double *c, Long l, Long m, Long n,
		       int VARIABLE_IS_NOT_USED cores) {
// multiplying A and t(B) with dim(A)=(l, m) and dim(B)=(n, m),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) if (l * m * n > 1000)
#endif
  for (Long i=0; i<l; i++) {
    double *C = c + i,
      *A = a + i;
    for (Long j=0; j<n; j++) {
       double tmp = 0.0,
	 *Bj = B + j;
       for (Long k=0; k<m; k++) tmp += A[k * l] * Bj[k * n];
       C[j*l] = tmp;
    }
  }
}


void matmult_2ndtransp(double *a, double *B, double *c, Long l, Long m,
		       int VARIABLE_IS_NOT_USED cores) {
// multiplying A and t(B) with dim(A)=(l, m) and dim(B)=(l, m),
// saving result in C
  Long lm = l  * m;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static) if (l * m * l > 1000)
#endif
  for (Long i=0; i<l; i++) {
    double *C = c + i,
      *A = a + i;
    for (Long j=0; j<l; j++) {
       double tmp = 0.0,
	 *Bj = B + j;
       for (Long k=0; k<lm; k+=l) tmp += A[k] * Bj[k];
       C[j*l] = tmp;
    }
  }
}



void Xmatmulttransposed(double *A, double *B, double *C, Long m, Long l, Long n,
			int VARIABLE_IS_NOT_USED cores) {
// multiplying t(A) and B with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static)
#endif
  for (Long i=0; i<l; i++) {
    Long im = i * m;
    for (Long jl=i, jm=0, j=0; j<n; j++, jl+=l, jm+=m) {
      double tmp = 0.0;
      Long endfor = im + m;
      for (Long jmk=jm, k=im; k<endfor; k++) tmp += A[k] * B[jmk++]; 
      C[jl] = tmp;
    }
  }
}



void matmult_tt(double *a, double *B, double *c, Long m, Long l, Long n,
		int VARIABLE_IS_NOT_USED cores) {
// calculating t(A B) with dim(A)=(m,l) and dim(B)=(m,n),
// saving result in C
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static)
#endif
  for (Long i=0; i<l; i++) {
    double *A = a + i,
      *C = c + i * l;
    for (Long j=0; j<n; j++) {
      double tmp = 0.0,
	*Bjm = B + j * m;
      for (Long k=0; k<m; k++) tmp += A[k * l] * Bjm[k];
      C[j] = tmp;
    }
  }
}


void matmulttransposedInt(int *A, int *B, Long *c,
			  Long trueNrow, Long trueNcolA, Long trueNcolB,
			  Long ldA,Long ldB,Long ldC,
			  int VARIABLE_IS_NOT_USED cores) {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(static)
#endif
  for (Long i=0; i<trueNcolA; i++) {    
    Long *C = c + i;
    int *Aim = A + i * ldA;
    for (Long j=0; j<trueNcolB; j++)
      C[j * ldC] = Ext_scalarXint(Aim, B + j * ldB, trueNrow, SCALAR_VERSION );
  }
}



int Match(char *name, name_type List, int n) {
  // == NOMATCHING, -1, if no matching function is found
  // == MULTIPLEMATCHING,-2, if multiple matching fctns are found,  
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  Ulong ln = STRLEN(name);
  int Nr=0;

  while ( Nr < n  && STRNCMP(name, List[Nr], ln)) Nr++;
  if (Nr < n) { 
    if (ln==STRLEN(List[Nr])) // exactmatching -- take first -- changed 1/7/07
      return Nr;
    // a matching function is found. Are there other functions that match?
    int j; 
    bool multiplematching=false;
    j=Nr+1; // if two or more covariance functions have the same name 
    //            the last one is taken 
    while (j<n) {
      while ( (j<n) && STRNCMP(name, List[j], ln)) {j++;}
      if (j<n) {
	if (ln==STRLEN(List[j])) { // exactmatching -- take first 
	  return j;
	}
	else {multiplematching=true;}
      }
      j++;
    }
    if (multiplematching) {return MULTIPLEMATCHING;}
  } else return NOMATCHING;
  return Nr;
}

int Match(char *name, const char * List[], int n, char sep, int *Len) {
  // printf("Matching\n");
   // == -1 if no matching name is found
  // == -2 if multiple matching names are found, without one matching exactly
  int ln,
    Nr=0;

  while ( Nr < n) {
    ln=0;
    while(name[ln] != '\0' && name[ln] == List[Nr][ln]) ln++;
    if (name[ln] == sep) {
      if (Len != NULL) *Len = ln - 1;
      if (List[Nr][ln] == '\0') return Nr;// exact match

      // name abbreviates an element within List:
      // so 1 partial match is found. Are there other elements that match?
      int j; 
      bool multiplematching=false;
      j=Nr+1; // if two or more covariance functions have the same name 
      //            the last one is taken 
      while (j<n) {
	while ( (j<n) && STRNCMP(name, List[j], ln)) j++;
	if (j<n) {
	  size_t tmp = STRLEN(List[j]);
	  stopIfNotInt(tmp);
	  if (ln==(int) tmp) return j;  // exactmatching -- take first 
	  else multiplematching=true;
	}
	j++; // now, only chance not to throw error is to find an exact match
      }
      return multiplematching ? MULTIPLEMATCHING : Nr;
    } else if (name[ln] == '\0') return NOSEPMATCH;
    Nr++;
  }
    
  return NOMATCHING;
}

int Match(char *name, const char * List[], int n) {
  return Match(name, List, n, '\0', NULL);
}


void matchError(char *varname, char* name, int intvalue, 
		const char * List[], int n, int errtype) {
  char msg0[1000];
  switch(errtype) {
  case 0 :
    if (name != NULL)
      SPRINTF(msg0, "'%.50s': unknown value '%.50s'. Possible values are:",
	      varname, name);
    else 
      SPRINTF(msg0,
	      "'%.50s':  value '%d' not in {0,...%d}. Possible values are:",
	      varname, intvalue, n-1);
    int i;
    for (i=0; i<n-1; i++) {
      char msg[1000];
      SPRINTF(msg, "%.900s '%.50s',", msg0, List[i]);    
      STRCPY(msg0, msg);
    }
    RFERROR2("%.900s and '%.50s'.\n", msg0, List[i]);  
    FALLTHROUGH_OK; // non-existent
    
  case 1 :
    RFERROR1("'%.50s': no value given.", varname);
  default :
    RFERROR1("unknown type of matchError: %d\n", errtype);
  }
}


int Match(char *varname, char* name, const char * List[], int n,
	  char sep, int *Len, bool relax) {
  int ans = Match(name, List, n, sep, Len);
  if (ans < 0) {
    if (name[0] == '\0' || (name[0]==' ' && name[1]=='\0')) {
      if (relax) return NOMATCHING;
      matchError(varname, name, ans, List, n, 1);      
    }
    matchError(varname, name, ans, List, n, 0);
  }
  return ans;
}



double ownround(double x) { return TRUNC((x + SIGN(x) * 0.5)); }


double lonmod(double x, double modulus) {  
  double 
    halfmodulus = 0.5 * modulus,
    y = x + modulus + halfmodulus;
  return Mod(y, modulus) - halfmodulus;
}


