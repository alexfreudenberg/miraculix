
/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2015 - 2017 Martin Schlather, Reinhard Furrer, Martin Kroll
 Copyright (C) 2017 - 2020 Martin Schlather
 Copyright (C) 2021 - 2022 Martin Schlather, Alexander Freudenberg
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


#include "Basic_RFUlocal.h"
#include "compatibility.lapack.h"
#include "extern_RFU.h"
#ifdef USEGPU
#include "solve_gpu.h"
#else
SIMD_MISS(solve_61, gpu);
#endif

#if ! defined FCONE // only as RFU has errror that R V 3x is used
  #if defined __GCC__
    #warning FCONE not available
  #endif 
  #define FCONE
#endif
#ifndef LOCAL_ERRORSTRING
#define LOCAL_ERRORSTRING
#endif
#ifdef WHICH_ERRORSTRING
#undef WHICH_ERRORSTRING
#define WHICH_ERRORSTRING pt->err_msg
#endif
#include "kleinkram.h"
#include "zzz_RFU.h"

//#include "xport_import_RFU.h"
extern const char * solve[solveN];
#define SCALAR(A,B,C) scalarX(A,B,C,NR)
#define LINEAR(A,B,C,D) linearX(A,B,C,D, linear_avx)


AVAILABLE_SIMD



const char * InversionNames[nr_InversionMethods] = {
  "cholesky", "svd", "eigen", "sparse",
  "method undefined",
  "qr", "lu", 
  "no method left",
  "GPU-cholesky",
  "R chol implementation",
  "direct formula",
  "diagonal"};
		

#define KAHAN false // O P TIONS.installNrun.kahanCorrection

#define C_MALLOC(INT, WHICH, N, TYPE)	{		\
    INT _N_ = N;							\
    if (pt->n_##WHICH < _N_) {						\
      if (pt->n_##WHICH < 0) BUG;					\
      FREE(pt->WHICH);							\
      pt->n_##WHICH = _N_;						\
      if ((pt->WHICH = (TYPE *) CALLOC(_N_, sizeof(TYPE))) == NULL)	\
	return ERRORMEMORYALLOCATION;					\
    } else {								\
      assert( (_N_ > 0 && pt->WHICH != NULL) || _N_ == 0);		\
      for (INT iii=0; iii<_N_; pt->WHICH[iii++] = 0);			\
    }									\
  }									\
    TYPE VARIABLE_IS_NOT_USED *WHICH = pt->WHICH

#define CMALLOC(WHICH, N, TYPE)  C_MALLOC(int, WHICH, N, TYPE)


/*
//  sqrtPosDef nutzt pt->U fuer das Ergebnis		
#define FREEING_RESULT(WHICH)					\
  assert(int VARIABLE_IS_UNUSED *_i = WHICH);		\
  if (pt->WHICH != NULL && pt->WHICH != result) {	\
    UNCONDFREE(pt->WHICH);				\
    pt->n_##WHICH = 0;					\
  }		
*/				
 	

double Determinant(double *M, int size, bool logarithm) {
  Long sizeP1 = size + 1,
    sizeSq = (Long) size * size;
  if (logarithm) {
    double tmp = 0.0;
    for (Long i=0; i<sizeSq; i+=sizeP1) tmp += LOG(M[i]);
    return tmp;
  } 
  double tmp = 1.0;
  for (Long i=0; i<sizeSq; i+=sizeP1) tmp *= M[i];
  return tmp;
}

double DeterminantLU(double *M, int size, bool logarithm, int *permut) {
  Long sizeP1 = size + 1;
    //    sizeSq = size * size;
  //  for (int i=0;i<sizeSq; i++) {  if (i % size == 0)printf("\n"); printf("%f ", M[i]);  }printf("\n");
  if (logarithm) {
    double tmp = 0.0;
    int vorz= 0;
    for (int i=0; i<size; i++) {
      tmp += LOG(FABS(M[i * sizeP1]));
      vorz += (permut[i] != i+1) + (M[i * sizeP1] < 0);
      //   printf("tmp %f %d %d %d\n", M[i * sizeP1], vorz , i, permut[i]);
    }
    //    printf("RFU det = %f %d\n", tmp, vorz);
    if (vorz % 2)
      RFERROR("calculation of log determinant need positive determinant");
    return tmp;
  } 
  double tmp = 1.0;  
  for (int i=0; i<size; i++) {
    tmp *= M[i * sizeP1];
    if (permut[i] != i+1) tmp = -tmp;
  }
  //  printf("RFU det = %f\n", tmp);
  return tmp;
}


double cumProd(double *D, int size, bool logarithm) {
  if (logarithm) {
    double dummy = 0.0;
    for (int i = 0; i < size; dummy += LOG(D[i++]));
    return dummy;
  }
  double dummy = 1.0;
  for (int i = 0; i < size; dummy *= D[i++]);
  return dummy;
}

#define FREEING(WHICH) assert(int VARIABLE_IS_UNUSED *_i=WHICH); FREE0(x,WHICH)
void solve_DELETE0(solve_storage *x) {
  FREEING(iwork);

  FREEING(pivotsparse);
  FREEING(pivot_idx);
  FREEING(xlnz);
  FREEING(snode);
  FREEING(xsuper);
  FREEING(invp);
  FREEING(cols);
  FREEING(rows);
  FREEING(lindx);
  FREEING(xja);
   
  FREEING(main);
  FREEING(rhs);
  FREEING(w2);
  FREEING(U);
  FREEING(D);
    
  FREEING(w3);
  FREEING(result);
  FREE(x->to_be_deleted);
}

void solve_DELETE(solve_storage **S) {
  solve_storage *x = *S;
  if (x!=NULL) {
    solve_DELETE0(*S);
    UNCONDFREE(*S);
  }
}
void solve_NULL(solve_storage* x) {
  if (x == NULL) return;
  MEMSET(x, 0, sizeof(solve_storage));  
  x->nsuper = x->size = -1;
  x->method = NoInversionMethod;
  for (int i=0; i<SOLVE_METHODS; x->newMethods[i++] = NoInversionMethod);
  x->actual_pivot = PIVOT_UNDEFINED;
}

double static inline det3(double *M, int size) {
  double det;
  switch(size){ // Abfrage nach Groesse der Matrix M + Berechnung der Determinante per Hand
  case 1: det = M[0];
    break;
  case 2: det = M[0] * M[3] - M[1] * M[2];
    break;
  case 3: det = 
      M[0] * (M[4] * M[8] - M[5] * M[7]) 
      - M[1] * (M[3] * M[8] - M[5] * M[6]) 
      + M[2] * (M[3] * M[7] - M[4] * M[6]); // Entwicklung nach 1. Spalte
    break;
  default : BUG;
    break;
  }
  return det;
}

int logdet3(double det, bool posdef, double *logdet, bool logarithm) {
  if (posdef && det < 0) return ERRORFAILED;
  if (logdet != NULL) {
    if (logarithm) {
      if (det <= 0) return ERRORFAILED;
      *logdet = LOG(det);
    } else *logdet = det;
  }
  return NOERROR;
}

int solve3(double *M, int size, bool posdef,
	   double *rhs, int rhs_cols,
	   double *result, double *logdet, bool logarithm,
	   solve_storage *pt	
	   ){
  
  assert(size <= 3);
  if (size <= 0) SERR0("matrix in 'solvePosDef' of non-positive size.");

  double det = det3(M, size);
  if (logdet3(det, posdef, logdet, logarithm) != NOERROR) return ERRORFAILED;

  double detinv = 1.0 / det; // determinant of inverse of M
  
  switch(size){
  case 1 : {// size of matrix == 1
    if (rhs_cols == 0) result[0] = detinv;   
    else for (int i=0; i<rhs_cols; i++) result[i] = rhs[i] * detinv;
  }
    break;
  case 2 : { // size of matrix == 2
    double a = M[0] * detinv,
      d = M[3] * detinv;
    if (rhs_cols == 0) {
      result[0] = d;
      result[1] = -M[1] * detinv;
      result[2] = -M[2] * detinv;
      result[3] = a;
    } else { // rhs_cols != 0
      double *p = rhs, *q = result;
      if (M[1] != 0.0 || M[2] != 0.0) {
	double 
	  b = M[1] * detinv,
	  c = M[2] * detinv;
	for (int i=0; i<rhs_cols; i++, p+=2, q+=2) {
	  double swap = d * p[0] - c * p[1];
	  q[1] = a * p[1] - b * p[0];
	  q[0] = swap;
	}
      } else {
	for (int i=0; i<rhs_cols; i++, p+=2, q+=2) {
	  double swap = d * p[0];
	  q[1] = a * p[1];
	  q[0] = swap;
	}
      }
    }
  }
    break;
  case 3 : {// size of matrix == 3   
    double swap0 = detinv * (M[4] * M[8] - M[5] * M[7]),
      swap1 = detinv * (M[5] * M[6] - M[3] * M[8]),
      swap2 = detinv * (M[3] * M[7] - M[4] * M[6]),
      swap3 = detinv * (M[2] * M[7] - M[1] * M[8]),
      swap4 = detinv * (M[0] * M[8] - M[2] * M[6]),
      swap5 = detinv * (M[1] * M[6] - M[0] * M[7]),
      swap6 = detinv * (M[1] * M[5] - M[2] * M[4]),
      swap7 = detinv * (M[2] * M[3] - M[0] * M[5]),
      swap8 = detinv * (M[0] * M[4] - M[1] * M[3]);
    if(rhs_cols == 0){ // invert matrix
      result[0] = swap0;
      result[1] = swap1;
      result[2] = swap2;
      result[3] = swap3;
      result[4] = swap4;
      result[5] = swap5;
      result[6] = swap6;
      result[7] = swap7;
      result[8] = swap8;
    } else { // solve system given by M and rhs
      double *p = rhs, *q = result;
      for (int i=0; i<rhs_cols; i++, p+=3, q+=3) {
	double swapA = p[0] * swap0 + p[1] * swap3 + p[2] * swap6;
	double swapB = p[0] * swap1 + p[1] * swap4 + p[2] * swap7;
	q[2] = p[0] * swap2 + p[1] * swap5 + p[2] * swap8;
	q[0] = swapA;
	q[1] = swapB;
      }
    }
  }
    break;
  default: BUG;
  }
  
  return NOERROR;
}

int chol3(double *M, int size, double *res, solve_storage *pt){
  // UNBEDINGT in sqrtRHS.cc auch aendern
  assert(size <= 3);
  if (size <= 0) SERR0("matrix in 'solvePosDef' not of positive size.");
  //  if (M[0] < 0) return ERRORFAILED;
  res[0] = SQRT(M[0]);
  if (size == 1) return NOERROR;
  res[1] = 0.0;
  res[size] = res[0] > 0.0 ? M[size] / res[0] : 0.0;
  double dummy = M[size + 1] - res[size] * res[size];
  res[size + 1] = SQRT(MAX(0.0, dummy));
  if (size == 2) return NOERROR;
  res[2] = res[5] = 0.0;
  res[6] = res[0] > 0.0 ? M[6] / res[0] : 0.0;
  res[7] = res[4] > 0.0 ? (M[7] - res[3] * res[6]) / res[4] : 0.0;
  dummy = M[8] - res[6] * res[6] - res[7] * res[7];
  res[8] = SQRT(MAX(0.0, dummy));
  return NOERROR;
}


void Sort(double *RESULT, int size, Long rhs_cols, int *pi, int *rank,
	  double *dummy) {
  if (size > MAXINT) BUG;
  orderingInt(pi, (int) size, 1, rank);
  int i=0;
  Long totalRHS = size * rhs_cols;
  while(i < size && i == rank[i]) i++;
  while (i < size) {
    int stored_i = i,
      read = i;
    double *p_write = NULL,
      *p_read = RESULT + read;
    for (Long  k=0; k<rhs_cols; k++) dummy[k] = p_read[k * size];
    while (true) {
      int write = read;
      p_write = RESULT + write;
      read = rank[write];
      rank[write] = write;
      if (read == stored_i) {
	for (int k=0; k<rhs_cols; k++) p_write[k*size]=dummy[k];
	break;
      }
      p_read = RESULT + read;
      for (Long k=0; k<totalRHS; k+=size) p_write[k] = p_read[k];
    }
    while(i < size && i == rank[i]) i++;
  }
}


// to do: fehlt: chol2solve(chol, x)


void chol2inv(double *MPT, int size, int cores) {
  Long sizeP1 = size + 1,
    sizeSq = (Long) size * size;
  int
    NR = KAHAN ? SCALAR_KAHAN : SCALAR_AVX;

  double *diagonal = (double *) MALLOC(sizeof(double) * (Long) size);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (size > 60) schedule(dynamic, 20)
#endif	    
  for (Long k=0; k<size; k++) {
    double *p_RESULT = MPT + k * sizeP1,
      diagK = diagonal[k] = 1.0 / p_RESULT[0];
    for (Long i=1; i<size - k; i++) {
      double *pM = p_RESULT + i * size;
      p_RESULT[i] = (-diagK * pM[0] - SCALAR(pM + 1, p_RESULT + 1, i -1))
	/ pM[i];
    }
    // i == k
  }
  
  int linear_avx = LINEAR_AVX * avxAvail;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (size > 60) schedule(dynamic, 20)
#endif
  for (Long k=0; k<size; k++) {	      
    double *p_RESULT = MPT + k * size;
    for (Long i=size-1; i>k; i--) {
      double *pM = MPT + i * size,
	r = (p_RESULT[i] /= pM[i]);
      diagonal[k] -= r *pM[k];
      LINEAR(pM + k + 1, -r, i-k-1, p_RESULT + k + 1);
      // for (int j=k+1; j<i; j++) p_RESULT[j] -= r * pM[j]; 
    }
    // i == k
  }
  
  for (Long k=0; k<size; k++) {
    double *pM = MPT + k * size;	         
    pM[k] = diagonal[k] / pM[k];
  }
   
  for (Long i2=0,i=0; i<size; i++, i2+=sizeP1) {
    Long i3 = i2 + 1;
    for (Long j = i2 + size; j<sizeSq; j+=size) 
      MPT[j] = MPT[i3++];
  }
  FREE(diagonal);

}





int doPosDefIntern(double *M0, int size, bool posdef,
		   double *rhs0, Long rhs_cols, double *result,
		   double *logdet, int calculate, solve_storage *Pt,
		   solve_options *sp, int VARIABLE_IS_NOT_USED cores
		   ){
  ASSERT_SOLVE(sp);
  
  // it is ok to have
  // ==15170== Syscall param sched_setaffinity(mask) points to unaddressable byte(s)
  // caused by gcc stuff

  //  printf("doPosDef %ld %d %d rhs=%ld %d %ld %ld calc=%d %ld %ld\n",
  //	 M0, size, posdef, rhs, rhs_cols, result, logdet, calculate, Pt,sp);
  /*
    M0: (in/out) a square matrix (symmetry is not checked) of size x size;
    NOTE THAT THE CONTENTS OF M0 IS DESTROYED IFF NO RHS IS GIVEN
    AND result IS NOT GIVEN.
    In case rhs is not given, the inverse of M0 is returned 
    In case of sqrtonly, M0 is expected to be a positive definite matrix
    posdef (in): whether or not the matrix is positiv (semi)definite --
    to some extend doPosDef can deal with non-positiv definite
    functions

    rhs (in/out) : right hand side of the equality with rhs_cols columns
    NOTE THAT THE CONTENTS OF rhs WILL BE DESTROYED IF rhs IS GIVEN, BUT NOT
    result  
 
    the solution of the equality is returned in rhs
    rhs_cols : number of colums of the matrix on the right hand side
    result (out) : NULL or matrix of the size of the result (inverse matrix or 
    of size of the matrix on the right hand side);see also 'M0' and 'rhs'
    logdet (out): if not NULL the logarithm of the determinant is returned
    pt (in/out) : working space. If NULL, internal working spaces are used.
 
    A non-NULL value gives an advantage only if doPosDef is called
    repeatedly. Then 
    solve_storage *pt = (solve_storage*) m alloc(sizeof(solve_storage);
    solve_NULL(pt);
    prepares pt. Deletion is done only at the very end by
    solve_DELETE(pt);
    In meanwhile the working space is managed by doPosDef;
    Sp (in): parameters. If NULL, the behaviour is as described in
    the R manual for doPosDef.
    The parameters are described in the package 'RandomFields'
    under ?RFoptions
	  
  */

  //  printf("entering\n");
  // http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F01/f01intro.xml#
  
  if (rhs_cols > MAXINT) BUG;
  assert(NA_LOGICAL == INT_MIN && NA_LOGICAL == NA_INTEGER); // nur zur sicherheit, wegen usr_bool
  //          eigentlich sollte usr_bool unabhaengig davon funktionieren
  assert(calculate != DETERMINANT ||
	 (logdet != NULL && result == NULL && rhs0 == NULL));
  assert(calculate != MATRIXSQRT || (rhs0 == NULL && posdef));
  assert((rhs_cols != 0) xor (rhs0 == NULL));

  double *RESULT = result != NULL ? result : rhs_cols > 0 ? rhs0 : M0;
  
  // Pivot_Cholesky:
  //	if (MPT == Morig || (rhs_cols > 0 && rhs == RESULT))
  //	  CERR0("Pivoted cholesky cannot be performed on place! Either you are a programmer or you should contact the maintainer.");

  // !MATRIXSQRT &&  rhs_cols > 0
  // assert(rhs != RESULT);
  la_modes la_mode = OPTIONS.installNrun.la_mode; 
  if (size <= sp->tinysize) {
    if (Pt != NULL) {
      Pt->method = direct_formula;
      Pt->size = size;
    }
    if (calculate == DETERMINANT)
      return logdet3(det3(M0, size), posdef, logdet, sp->det_as_log);
    else if (calculate  == MATRIXSQRT) return chol3(M0, size, RESULT, Pt);
    else return solve3(M0, size, posdef, rhs0,
		       (int) rhs_cols, RESULT, logdet,
		       sp->det_as_log, Pt);
  }

  //  printf("size=%d %d %d sparse=%d\n", size, sp->pivot, PIVOT_AUTO, sp->sparse);
  // ERR0("XXXX");
  // printf("%d %d %d\n",  PIVOT_AUTO, direct_formula,Pt->method);

  assert(SOLVE_METHODS >= 2);
 
  //  printf("A\n");
 
  //  printf("A1\n");
  int  err = NOERROR;
  Long VARIABLE_IS_NOT_USED spam_zaehler = 0; // if spam is not compiled
  int nnzA = 0;
  Long
    sizeSq = (Long) size * size,
    sizeRHS = (Long) size * rhs_cols;
  int sizeP1 = size + 1;
  //  printf("A1dd\n");
#if defined compatibility_to_C_h
   usr_bool sparse = False;
#else
    usr_bool sparse = sp->sparse;
#endif
  double spam_tol = sp->spam_tol;
  //   printf("A2\n");
  bool diag  = false,
    useGPU = la_mode == LA_GPU &&
    (calculate == SOLVE || calculate == DETERMINANT);
 //   printf("A3\n");


  //  printf("A %d %d size=%d %d %d \n",	 sparse,Nan ,size, useGPU, sp->n_spam_min[useGPU]);
  if (sparse == Nan && (sparse = (usr_bool) (size > sp->spam_min_n[useGPU]))) {
    //     printf("AB2\n");
    double mean_diag = 0.0;
    for (Long i=0; i<sizeSq; i += sizeP1) mean_diag += M0[i];
    mean_diag /= (double) size;
    spam_tol *= mean_diag;
    // printf("AC2\n");
    bool random_sample = sizeSq >= sp->spam_sample_n * 3;
    if (random_sample) {
      //  printf("A2E\n");
      double 
	thr = sp->spam_sample_n * (1.0 - sp->spam_min_p[useGPU]);
      Long
	threshold = (Long) (thr + SQRT(thr) * 3),
	notZero = 0;
      for (Long i=0; i<sp->spam_sample_n; i++) {
	//	printf("A2 %d %d\n", i , sp->spam_sample_n);
	if ((notZero += !(FABS(M0[(i * sp->spam_factor) % sizeSq]) <=
			  spam_tol)) >= threshold){
	  sparse = False;
	  break;
	}
      }
      if (PL >= PL_FCTN_DETAILS) {
	PRINTF("random sampling: sparse=%d\n",
	       sparse == Nan ? NA_INTEGER : (int) sparse);
      }
    }
    ///    printf("A2EX %d %d\n", random_sample, sparse == True);
    if (!random_sample || sparse == True) {
      //   printf("AdC2\n");
      Long diag_nnzA = 0;
      //#ifdef DO_PARALLEL
      //#pragma omp parallel for num_threads(cores) schedule(dynamic,10) reduction(+:nnzA,diag_nnzA)
      //#endif
      for (Long i=0; i<size; i++) {
	//	printf("AC2 %d %d\n",i ,size);
	Long end = i * sizeP1;
	Long j;
	// Achtung!!: !(FABS(M0[j++]) <= spam_tol) != FABS(M0[j++]) > spam_tol
	//             bei NA/NaN
	for (j=i * size; j<end; nnzA += !(FABS(M0[j++]) <= spam_tol));
	diag_nnzA += !(FABS(M0[j++]) <= spam_tol);
	end = (i+1L) * size;
	if (!posdef) for (; j<end; nnzA += !(FABS(M0[j++]) <= spam_tol));
      }
      //     printf("AddC2\n");
      diag = (nnzA == 0);
      if (posdef) nnzA *= 2;
      stopIfNotUInt(diag_nnzA);
      nnzA += (int) diag_nnzA;
      sparse = (usr_bool) ((double) nnzA <=
			   (double) sizeSq * (1.0 - sp->spam_min_p[useGPU]));
      spam_zaehler = nnzA + 1;
      //    printf("ddAdC2\n");
     if (PL >= PL_DETAILSUSER) {
	if (diag) { PRINTF("diagonal matrix detected\n"); }
	else if (sparse == True) {
	  PRINTF("sparse matrix detected (%3.2f%% zeros)\n", 
		 100.0 * (1.0 - (double) nnzA / (double) sizeSq));
	} else {
	  PRINTF("full matrix detected (%3.2f%% nonzeros)\n", 
		 100.0 * (double) nnzA / (double) sizeSq); }
      }
    }
  } else {
    //   printf("ttAdC2\n");
    diag = true;
    for (Long i=0; i<size && diag; i++) {
      //      printf("AdC2 %d %d %d\n", i, size, diag);
      Long j,
	end = i * sizeP1;
      for (j=i * size; j<end; j++) {
	//	printf("(%d %d %10g %d)\n", i, j, M0[j], size);
	if (!(FABS(M0[j]) <= spam_tol)) {
	  diag = false;
	  break;
	}
      }
      if (!diag) break;
      j++;
      end = (i+1) * size;
      if (!posdef) {
	for (; j<end; j++) {
	  if (!(FABS(M0[j]) <= spam_tol)) {
	    diag = false;
	    break;
	  }
	}
      }
    }
  }


  solve_storage *pt;
  if (Pt != NULL) {
    pt = Pt; 
  } else {
    pt = (solve_storage*) MALLOC(sizeof(solve_storage));
    solve_NULL(pt);    
  }
  InversionMethod *Meth = pt->newMethods;
  pt->method = NoFurtherInversionMethod;
  pt->size = size;

  
  //  printf("BA\n");
  if (diag) {
    pt->method = Diagonal;
    if (PL>=PL_STRUCTURE) { PRINTF("dealing with diagonal matrix\n"); }
    if (logdet != NULL) {
      *logdet = Determinant(M0, size, sp->det_as_log);
      if (calculate == DETERMINANT) {
	err = NOERROR;
	goto ErrorHandling;
      }
    }
    if (rhs_cols == 0) {
      MEMCOPY(RESULT, M0, sizeSq * sizeof(*RESULT));
      if (calculate == MATRIXSQRT) {
	for (Long i=0; i<sizeSq; i += sizeP1) {
	  RESULT[i] = M0[i] > 0.0 ? SQRT(M0[i]) : 0.0;	
	}
      } else {
	for (Long i=0; i<sizeSq; i += sizeP1) 
	  RESULT[i] = M0[i] <= 0.0 ? 0.0 : 1.0 / M0[i];
      }
    } else {
      CMALLOC(main, size, double);
      for (Long i=0; i<size; i++) {
	Long idx = i * sizeP1;
	main[i] = M0[idx] == 0.0 ? 0.0 : 1.0 / M0[idx];
      }
      Long j;
      for (Long k=j=0; j<rhs_cols; j++)
	for (Long i=0; i<size; i++, k++) RESULT[k] = rhs0[k] * main[i];
    }
    
    err = NOERROR;
    goto ErrorHandling;
  }

  
  // size of matrix at least 4 x 4, and not diagonal
  int to, from;
  from = 0;
  if ((to = sparse == True)) Meth[0] = Sparse;
  while (from < SOLVE_METHODS &&
	 to < SOLVE_METHODS &&
	 sp->Methods[from] != NoFurtherInversionMethod &&
	 sp->Methods[from] != NoInversionMethod) {
    if (sp->Methods[from] == Sparse && sparse == True) from++;
    else Meth[to++] = sp->Methods[from++];
  }

  //  printf("from = %d (%d %d) [%d %d %d] sparse=%d %d\n", from, sp->Methods[0], sp->Methods[1], Meth[0],  Meth[1], Meth[2], sparse == True, Sparse);
  
  if (from == 0) { // user did not give any method
    if (posdef) {
      if (to < SOLVE_METHODS) {
	Meth[to++] = useGPU ? GPUcholesky : Cholesky;
	if (to < SOLVE_METHODS) {
	  Meth[to++] = sp->pivot_mode != PIVOT_NONE && useGPU ? Cholesky : Eigen;
	}
      }
    } else {
      Meth[to++] = LU;
    }
  } else {
    if (useGPU) 
      for (int i=0; i<SOLVE_METHODS; i++)
	if (Meth[i] == Cholesky || Meth[i] == GPUcholesky){
	  Meth[i]=GPUcholesky;
	  break; // so "Cholesky" can be given several times, only
	  // the very first will be changed to GPUcholesky.
	  // if GPU is already given do not change any given Cholesky afterwards
	}
  }

  //  printf("to = %d (%d %d) %d %d\n", to, Meth[0], Meth[1], from, sparse);  assert(Meth[0] != 4);
  
  for (; to<SOLVE_METHODS; Meth[to++]=NoFurtherInversionMethod);//save

  
  
  // cholesky, QR, SVD, Eigen, LU always destroy original matrix M
  int first_not_reading_M0;
  double *MPT, *RHS;
  MPT = M0;// pointer of M matrix, die zerstoert wird
  RHS = rhs0; 
  first_not_reading_M0 = 0;

  if (rhs_cols == 0 && result != NULL) MPT = result;
  else {
    first_not_reading_M0 = sparse + (Meth[sparse] == GPUcholesky);
    // printf("first_not_reading_M0 %d %d\n", first_not_reading_M0, Meth[first_not_reading_M0] != Meth[first_not_reading_M0 - 1] &&
    //	 Meth[first_not_reading_M0] != NoFurtherInversionMethod);
    
    if (first_not_reading_M0 == 0 ||
	(Meth[first_not_reading_M0] != Meth[first_not_reading_M0 - 1] &&
	 Meth[first_not_reading_M0] != NoFurtherInversionMethod)) {
      //printf("rhs_cols +%d\n", rhs_cols);
      if (rhs_cols > 0
	  ||
	  (SOLVE_METHODS > first_not_reading_M0 + 1 &&
	   Meth[first_not_reading_M0 + 1] != Meth[first_not_reading_M0] &&
	   Meth[first_not_reading_M0 + 1] != NoFurtherInversionMethod)
	  ||
	  (Meth[first_not_reading_M0] == SVD && sp->svd_tol > 0.0 &&
	   calculate != SOLVE)
	  ) { // at least two different Methods in the list
	C_MALLOC(Long, main, sizeSq, double);//to pt->main, main local variable
	MPT = pt->main;
	if (rhs_cols > 0) {
	  C_MALLOC(Long, rhs, sizeRHS, double);//to pt->main, main local variab
	  RHS = pt->rhs;
	}
      }
    }
  }

  
  // printf("AFF\n");

  errorstring_type ErrStr;
  STRCPY(ErrStr, "");

  //  printf("Meth=%d %d Chol=%d %d posdef=%d\n", Meth[0], Meth[1], Cholesky, SOLVE_METHODS, posdef);

  //  for (int i=0; i<size; i++) { printf("\n");
  //   for (int j=0; j<size; j++)  printf("%e ", M0[i * size + j]);
  // } printf("\n");

 
  int proposed_pivot;
  proposed_pivot = sp->pivot_mode;
  for (int m=0; m<SOLVE_METHODS && (m==0 || Meth[m] != Meth[m-1]); m++) {
    pt->method = Meth[m];
    if (pt->method == Cholesky && size > OPTIONS.installNrun.LaMaxTakeIntern) {
      pt->method = calculate == DETERMINANT ? LU : Rcholesky;
    }
    if (pt->method < 0) break;
    if (calculate != SOLVE) {
      if (pt->method == NoInversionMethod && m<=sparse) BUG;
      if (pt->method == NoFurtherInversionMethod) break;
      if (PL>=PL_STRUCTURE) { 
	PRINTF("method to calculate the square root : %s\n", 
	       InversionNames[pt->method]);
      }
    } else {
      if (PL>=PL_STRUCTURE) { 
	PRINTF("method to calculate the inverse : %s\n",
	       InversionNames[pt->method]);
      }
    }
     
    if (MPT != M0 && m >= first_not_reading_M0)
      MEMCOPY(MPT, M0,  sizeSq * sizeof(*MPT));

    if (RHS != rhs0) MEMCOPY(RHS, rhs0,  (Long) sizeRHS * sizeof(*RHS));
 
    switch(pt->method) {
    case GPUcholesky :
      if (size > 28000)
	GERR0("due to a bug at NVIDIA the maximum size is currently limited to 28000.");
      
      if (!posdef) CERR0("Cholesky needs positive definite matrix");
#ifdef USEGPU
      if (proposed_pivot > PIVOT_AUTO)
	GERR0("cholesky decomposition on GPU does not allow for pivoting");
      pt->actual_pivot = PIVOT_NONE;
      {
	double LD,
	  *LogDet = logdet == NULL ? &LD : logdet;
	err = // int; see errors_messages.h for values 0...4
	  cholGPU(true, // bool : in : says that values must be copied
		  // so in miraculix, RFU_AlexChol(false, ....)
		  // can be called
		  M0,// in: this matrix is copied by Alex because
		  //               of value 'true' in the first argument,
		  //               so contents never distroyed by Alex
		  size, // in:  size of the matrix
		  rhs0, //in: if NULL the inverse of M is calculated;
		  //  rhs is copied by Alex  because of value 'true'
		  //  in the first argument,, so never distroyed by Alex
		  rhs_cols, // in: number of columns on the right hand side
		  LogDet, // out : logarithm of the determinant of
		  //               the sqare root(!) of the matrix M
		  RESULT); // out: a pointer to the result whether or
	// not rhs is given
	if (err != NOERROR) {
	  if (proposed_pivot == PIVOT_AUTO) proposed_pivot = PIVOT_DO;
	  continue;
	}
	if (logdet != NULL) {
	  *logdet *= 2;
	  if (!sp->det_as_log) *logdet = EXP(*logdet);
	  if (calculate == DETERMINANT) {
	    err = NOERROR;
	    goto ErrorHandling;
	  }
	}
      }
#else
      //     err = NOERROR;
      BUG;
#endif
      break;
    case Rcholesky : {
      //      printf("chol (R)\n");
      if (calculate == SOLVE) {
	if (rhs_cols > MAXINT) BUG;
        int n_rhs = (int) rhs_cols;
	double *mm = NULL;
	CMALLOC(xja, size, int);
	if (rhs_cols == 0) {
	  //	 printf("hier %ld %ld %ld res=%ld %ld %ld\n", RESULT, MPT, M0, result, rhs0, RHS);
	  if (RESULT == MPT){
	    Long bytes = sizeSq * sizeof(*RESULT);
	    mm = (double *) MALLOC(bytes);
	    MEMCOPY(mm, M0, bytes);
	  }
	  
	  MEMSET(RESULT, 0, sizeof(*RESULT) * sizeSq);
	  for (Long i=0; i<sizeSq; i+=sizeP1) RESULT[i] = 1.0;
	  if (size > MAXINT) BUG;
 	  n_rhs = size;
	  
	  /*
	    
	    for (Long i=0; i<size; i++)  {
	    for (Long j=0; j<size; j++) {
	    //printf("%f ", RESULT[i + j * size]); 
	    }
//	    printf("\n");
	    }
//	    printf("\n M\n");
	    for (Long i=0; i<size; i++)  {
	    for (Long j=0; j<size; j++) {
//	    printf("%f ", MPT[i + j * size]);
	    }
//	    printf("\n");
	    }
	    
	  */
	  
	} else {
	  if (RHS != RESULT)
	    MEMCOPY(RESULT, RHS, sizeof(*RESULT) * (Long) size * rhs_cols);
	}
	F77dgesv(&size, &n_rhs, mm == NULL ? MPT : mm, &size, xja, RESULT,
		 &size, &err);
	if (logdet != NULL)
	  *logdet = DeterminantLU(mm == NULL ? MPT : mm, size, sp->det_as_log,
				  xja);
	FREE(mm); // NOTE: return / errrors only afterwards!	
	if (err != NOERROR) GERR0("LU algorithm failed.");
      } else if (calculate == MATRIXSQRT) {
	if (MPT != RESULT) MEMCOPY(RESULT, MPT, sizeSq * sizeof(*RESULT));
        F77dpotrf("U", &size, RESULT, &size, &err
#ifdef USE_FC_LEN_T	  
		  FCONE
#endif		  
		  );
	if (logdet != NULL) {
	  Determinant(RESULT, size, sp->det_as_log);
	  if (sp->det_as_log) *logdet *=2; else *logdet *= *logdet;
	}
	Long sizeM1 = size - 1;
	for (Long i=0; i<sizeM1; i++)
	  MEMSET(RESULT + i * sizeP1 + 1, 0,
		 sizeof(*RESULT) * ((Long) size - i - 1));
	// inversion:
	// F77dpotri("U", &sz, REAL(ans), &sz, &info FCONE);
      } else BUG;
      if (err != NOERROR) {
	CERR1("algorithm failed (err=%d). Consider 'RFoptions(la_mode=\"intern\", pivot=PIVOT_AUTO)'", err);
      }
      break;
    }
    case Cholesky : {
      //      printf("chol\n");
      //
#ifdef DO_PARALLEL
      //   printf("chol (own), %d cores\n",  cores);
#else
      //   printf("chol (own), single core\n");
#endif
#define C_GERR0(X,G) {STRCPY(ErrStr, X); FERR0(X); err = ERRORM; goto G;}
#define C_GERR1(X,Y,G) {SPRINTF(ErrStr,X,Y); FERR0(ErrStr);err = ERRORM; goto G;}
#define C_GERR2(X,Y,Z,G){SPRINTF(ErrStr,X,Y,Z);FERR0(ErrStr);err=ERRORM; goto G;}
#define C_GERR3(X,Y,Z,A,G) {SPRINTF(ErrStr,X,Y,Z,A); FERR0(ErrStr);err = ERRORM; goto G;}
      int NR = KAHAN ? SCALAR_KAHAN : SCALAR_AVX,
	linear_avx = LINEAR_AVX * avxAvail;
      if (size > sp->max_chol) {
	CERR2("Matrix  is too large for Cholesky decomposition. Maximum ist currently a %d x %d matrix. Increase 'max_chol' in 'RFoption' if necessary.",
	      sp->max_chol, sp->max_chol);
      }

      ///      printf("size = %d %d %d\n", size, rhs_cols, size > rhs_cols ? size : rhs_cols);
      
      C_MALLOC(Long, D, size > rhs_cols ? size : rhs_cols, double);
      for (Long i=0; i<size; i++) D[i] = MPT[sizeP1 * i];
      
      pt->actual_pivot = PIVOT_UNDEFINED;
   
      if (proposed_pivot == PIVOT_NONE || proposed_pivot == PIVOT_AUTO) {

	// cmp for instance http://stackoverflow.com/questions/22479258/cholesky-decomposition-with-openmp

	// obere und untere dreiecksmatrix wird gelesen und obere geschrieben
	err = NOERROR;
	pt->actual_pivot = PIVOT_NONE;
	{
	  double *A = MPT;
	  for (Long i=0; i<size; i++, A += size) {
	    double sclr = SCALAR(A, A, i);
	    if (A[i] <= sclr) {
	      if (proposed_pivot == PIVOT_NONE)
		C_GERR2("Got %10e as %d-th eigenvalue.",
			A[i] - sclr, (int) i, ERR_CHOL)
		else C_GERR2("Got %10e as %d-th eigenvalue.",
			     A[i] - sclr, (int) i, Pivot_Cholesky);
	      break;
	    }
	    A[i] = SQRT(A[i] - sclr);

  
	    //	  double invsum = 1.0 / A[i];
	    double sum = A[i];
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (MULTIMINSIZE(size - i)) schedule(dynamic, 8) 
#endif
	    //	    for (double *B=MPT + (i+1) * size; B<endB; B+=size) {
	    for (Long j = i + 1; j < size; j++) {
	      double *B = MPT + j * size;
	      B[i] = (B[i] - SCALAR(A, B, i)) / sum;
	    }
	  }
	}
      
	if (err == NOERROR) {
	  if (calculate == MATRIXSQRT) {
	    int deltaend = size - 1;
	    double *end = MPT + sizeSq;
	    for (double *p=MPT + 1; p<end; p+=sizeP1, deltaend--)
	      FILL_IN(p, deltaend, 0.0);	    
	  } else {
	    if (logdet != NULL) {
	      *logdet = Determinant(MPT, size, sp->det_as_log);
	      if (sp->det_as_log) *logdet *=2; else *logdet *= *logdet;
	      if (calculate == DETERMINANT) {
		err = NOERROR;
		goto ErrorHandling;
	      }
	    }

	    if (rhs_cols == 0) chol2inv(MPT, size, cores);
	    else { // rhs_cols > 0
	      //Long totalRHS = size * rhs_cols;
	      //if (result!=NULL) MEMCOPY(RESULT, rhs, sizeof(*RESULT)*totalRHS);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (rhs_cols > 30) schedule(static)
#endif	      
	      for (Long k=0; k<rhs_cols; k++) {
		double *p_RESULT = RESULT + k * size,
		  *p_rhs = RHS + k * size;
		for (Long i=0; i<size; i++) {
		  double *pM = MPT + i * size;
		  p_RESULT[i] = (p_rhs[i] - SCALAR(pM, p_RESULT, i)) / pM[i];
		}
	      }
	      
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (rhs_cols > 30) schedule(static)
#endif	      
	      for (Long k=0; k<rhs_cols; k++) {	      
		double *p_RESULT = RESULT + k * size;
		for (Long i=size-1; i>=0; i--) {
		  double *pM = MPT + i * size,
		    r = (p_RESULT[i] /= pM[i]);
		  LINEAR(pM, -r, i, p_RESULT);
		  // for (Long j=0; j<i; j++)  p_RESULT[j] -= r * pM[j];
		}
	      }
	    }
	  } // not sqrt only
	} // err == NOERROR
      } else err = ERRORFAILED;
      
    
      Pivot_Cholesky:
      //      printf("Err pivot %d %s %d\n", err, ErrStr, proposed_pivot );
      if (err != NOERROR && proposed_pivot != PIVOT_NONE) {
	if (PL > PL_DETAILS) { PRINTF("trying pivoting\n"); }
	int actual_size = NA_INTEGER;
	// code according to Helmut Harbrecht,Michael Peters,Reinhold Schneider
	// talk: The pivoted Cholesky decomposition and its application to
	//       stochastic PDEs
	// ebenso: untere dreiecksmatrix wird gelesen; obere geschrieben
	if (pt->actual_pivot == PIVOT_NONE) {
	  // wiederherstellung der Diagonalen und der unteren dreiecksmatrix
	  for (Long i=0; i<size; i++) MPT[sizeP1 * i] = D[i];
	}
	int *pi;
	if (proposed_pivot == PIVOT_DO || proposed_pivot == PIVOT_AUTO) {
 	  FREE(pt->pivot_idx); // ALWAYS FREE IT!!! cp Chol(S E X P M)
	  pt->pivot_idx = (int*) MALLOC((Long) size * sizeof(int));
	  pt->n_pivot_idx = size;
	  pt->actual_pivot = PIVOT_DO;
	  for (int i=0; i<size; i++) pt->pivot_idx[i] = i;
	  pt->actual_size = actual_size = size;
	  pi = pt->pivot_idx;
	} else { // PIVOT_IDX 
	  if (sp->n_pivot_idx < size || sp->actual_size > size) {
	    //	    printf("XA, %d %d %d\n", sp->n_pivot_idx , size, sp->actual_size);
	    CERR0("pivot_idx does not have the correct length.\nSet 'RFoption(pivot_idx=, pivot_actual_size=)' to the attributes of a\npivoted Cholesky decomposition.");
	  }
	  actual_size = pt->actual_size = sp->actual_size;
	  if (actual_size > size) BUG;
	  FREE(pt->pivot_idx);
	  Long bytes = (Long) size * sizeof(*(pt->pivot_idx));
	  pt->pivot_idx = (int*) MALLOC(bytes);
	  MEMCOPY(pt->pivot_idx, sp->pivot_idx, bytes);
	  pt->actual_pivot = PIVOT_IDX;
	  pt->n_pivot_idx = sp->n_pivot_idx;
	  pi = sp->pivot_idx;
	}

	err = NOERROR;
	//	printf("hier\n");
	double // *rhs = rhs,
	  //	  *M00 = M0,
	  rel_thres = 0,
	  max_deviation = sp->max_deviation, //  1e-10,
	  max_reldeviation = sp->max_reldeviation; //  1e-10,

	//printf("MTP %d %d %d\n", MPT == M0, rhs_cols,  RHS == RESULT);
	if (MPT == M0 || (rhs_cols > 0 && RHS == RESULT))
	  CERR0("Pivoted cholesky cannot be performed on place! Either you are a programmer or you should contact the maintainer.");
	/* 
	if (MPT == M0) {
	  CMALLOC(main, sizeSq, double);
	  MEMCOPY(main, M0, sizeSq * sizeof(*main));
	  M00 = main;	  
	}
	if (rhs_cols > 0 && rhs == RESULT) {
	  Long totalRHS = (Long) size * rhs_cols;
	  CMALLOC(U, totalRHS, double);
	  MEMCOPY(U, rhs, totalRHS * sizeof(*U));
	  RHS = U;
	  }*/
	for (Long q=0; q<actual_size; q++) {
	  if (pt->actual_pivot == PIVOT_DO) {
	    double 
	      max = RF_NEGINF,
	      deviation = 0.0;
	    Long k,
	      argmax = NA_INTEGER;	    
	    for (k=q; k<size; k++) {
	      double dummy = D[pi[k]];
	      // if (D[pi[k]] < 0)

	      //printf("k=%d %10e %10e\n", k, du1mmy, -1e-15 * size * size);
	      if (dummy < -1e-4 * sp->pivot_relerror * (double) sizeSq){
		C_GERR1("matrix not positive definite or increase 'pivot_relerror' by at least factor %10g.", dummy * -1e4 / (double) sizeSq, ERR_CHOL);
	      }
	      deviation += dummy;
	      if (max < dummy) {
		max = dummy;
		argmax = k;		    
	      }
	    }
	    
	    double dev = rel_thres * max_reldeviation;
	    if (deviation <= max_deviation || (q > 0 && deviation <= dev) ) {
	      if (q > MAXINT) BUG;
	      actual_size = pt->actual_size = (int) q;
	      if (sp->pivot_check != False) {
		double largest = 0;
		for (Long i=q; i<size; i++) {
		  double *mpt = MPT + pi[i] * (Long) size;
		  for (Long j=q; j<=i; j++) {
		    double absm = FABS(mpt[j]);
		    largest = absm > largest ? absm : largest;
		    // if(absm == largest || absm > 5) printf("%10e %d %d; %d\n", absm, i, j, size);
		  }
		}
		if (largest > max_deviation || (q > 0 && largest > dev)) {
		  char msg[500];
		  SPRINTF(msg, "Matrix has a numerically zero, leading value at the %d-th pivoted line, but the largest deviation %10e from zero in the rest of the matrix is greater than the tolerance %10e. %.50s.",
			  (int) q,
			  largest,
			  MAX(max_deviation, dev),
			  sp->pivot_check == True
			  ? "If you are sure that the matrix is semi-definite, set 'RFoptions(pivot_check=NA)' or 'RFoptions(pivot_check=True)'"
			  : "The result can be imprecise");
		  if (sp->pivot_check == True) C_GERR0(msg, ERR_CHOL)
		    else WARN0(msg);
		}
	      }
	      break;
	    }
	    rel_thres += D[pi[q]];
	    int dummy = pi[q];
	    pi[q] = pi[argmax];
	    pi[argmax] = dummy;
	  }
	  
	  
	  Long pqq = pi[q],
	    col_q = pqq * size;

	  
	  if (D[pqq] < 0) {
	    C_GERR1("Negative leading value found at the %d-th pivoted line.",
		    (int) q, ERR_CHOL);
	  }
	  double lqpq = MPT[q + col_q] = SQRT(D[pqq]);	    

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (MULTIMINSIZE(size - q)) schedule(dynamic, 8) 
#endif
	  for (Long i=q+1; i<size; i++) {
	    Long pii = pi[i];
	    Long col_i = pii * size;
	    assert(pii != pqq);
	    double scalar = SCALAR(MPT + col_q, MPT + col_i, q);
	    MPT[q + col_i] = (M0[pqq + col_i] - scalar) / lqpq;
	    D[pii] -=  MPT[q + col_i] *  MPT[q + col_i];
	    // in Harbrecht: lqpq * MPT[q + col_i];
	  }
	  
	  /*	  if (!true) {
		  for (Long k=0; k<size; k++) {  
		  for (Long j=0; j<size; j++)
		  p rintf("%10g ", MPT[size * j + k]);
		  p rintf("\n");
		  } p rintf("\n");
		  }
	  */
	  
	} // for q

	
	if (err == NOERROR) {
	  if (calculate == MATRIXSQRT) {
	    Long i = 0;
	    for ( ; i<actual_size; i++) {
	      FILL_IN(MPT + i + 1 + size * pi[i], size - 1 - i, 0.0);
	    }
	    for ( ; i<size; i++) {
	      FILL_IN(MPT + actual_size + (Long) size*pi[i], size - actual_size, 0.0);
	    }
	    
	  } else { // ! MATRIXSQRT
	    if (logdet != NULL) {
	      Long N = size;
	      if (pt->actual_pivot == PIVOT_DO) {
		N = actual_size;
	      }
	      if (sp->det_as_log) {
		if (N < size && !sp->pivot_partialdet) *logdet = RF_NEGINF;
		else {
		  double logD = 0.0;
		  for (Long i=0; i < N; i++)
		    logD += LOG(MPT[i + pi[i] * (Long) size]);
		  *logdet = logD * 2;
		}
	      } else {
		if (N < size && !sp->pivot_partialdet) *logdet = 0;
		else {
		  double logD = 1.0;
		  for (Long i=0; i < N; i++)
		    logD *= MPT[i + pi[i] * (Long) size];
		  *logdet = logD * logD;
		}
	      }
	      if (calculate == DETERMINANT) {
		err = NOERROR;
		goto ErrorHandling;
	      }
	    }
	    
	    //////////////////////////////////////////////////
	    //////////////////////////////////////////////////

	    if (rhs_cols == 0) {
 	      if (actual_size < size) 
		GERR0("Matrix not definite. Try ")
	      
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (size > 60) schedule(dynamic, 20)
#endif
		  for (Long k=0 ; k<actual_size; k++) {
		    double *p_RESULT = MPT + pi[k] * (Long) size + k,
		      diagK = D[k] = 1.0 / p_RESULT[0];
		    for (Long i=1; i<size - k; i++) {
		      double *pM = MPT + k + pi[k + i] * (Long) size;
		      p_RESULT[i] = (-diagK * pM[0]
				     -SCALAR(pM + 1, p_RESULT + 1, i -1)) / pM[i];
		    }
		  }     
	      
	      
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (size > 60) schedule(dynamic, 20)
#endif	    
	      for (Long k=0; k<size; k++) {	      
		double *p_RESULT = MPT + pi[k] * (Long) size;
		for (Long i=actual_size-1; i>k; i--) {
		  double *pM = MPT + pi[i] * (Long) size,		    
		    r = (p_RESULT[i] /= pM[i]);
		  D[k] -= r * pM[k];
		  LINEAR(pM + k + 1, -r, i-k-1, p_RESULT + k + 1);
		  // for (Long j=k+1; j<i; j++) p_RESULT[j] -= r * pM[j]; 
		}
		// i == k
	      }
	      
	      for (Long k=0; k<size; k++) {	      
		double *p_RESULT = MPT + pi[k] * (Long) size;
		for (Long q=0; q<k; p_RESULT[q++] = RF_INF);
		//	for (Long q=actual_size; q<size; p_RESULT[q++] = RF_NA);
	      }
	      
	      for (Long k=0; k<actual_size; k++) {
		double *pM = MPT + pi[k] * (Long) size;	         
		pM[k] = D[k] / pM[k];
	      }
	      
	      CMALLOC(xja, size, int);
	      Sort(RESULT, size, size, pi, xja, D);
	      for (Long i=0; i<size; i++) {
		for (Long j=i+1; j<size; j++) {
		  Long idx = i + j * size;
		  if (MPT[idx] == RF_INF) MPT[idx] = MPT[j + i * size];
		  else MPT[j + i * size] = MPT[idx];
		}
	      }
	      

	      //////////////////////////////////////////////////
	      //////////////////////////////////////////////////
	      
	    } else { // rhs_cols > 0
	      // if (rhs0 == RESULT) {
	      //	/* crash(); */
	      //	#pragma GCC diagnostic push
	      //#pragma GCC diagnostic ignored "-Wuninitialized"
	      // Long i; PRINTF("%d\n", i);char m[1];m[i] = m[i-9] + 4; if (m[0]) i++; else i--; PRINTF("%s\n", m); // not MEMCOPY
	      //#pragma GCC diagnostic pop
  //  int *x = (int*) MALLOC(1000000);  f ree(x);  f ree(x); x[100] = 100;
	      //    }
	      //  printf("%ld %ld %ld\n", rhs0 ,  RESULT, RHS);
	      //assert(rhs0 != RESULT); 
	      // assert(RHS != RESULT); 
	      double eps = D[0] * sp->pivot_relerror;
			      
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (rhs_cols > 30) schedule(static)
#endif
	      for (Long k=0; k<rhs_cols; k++) {
		double *p_RESULT = RESULT + k * size,
		  *p_rhs = RHS + k * size;
		Long i=0;
		for ( ; i<actual_size; i++) {
		  Long pii = pi[i];
		  double *pM = MPT + pii * size;
		  p_RESULT[i] = (p_rhs[pii] - SCALAR(pM, p_RESULT, i)) / pM[i];
		}		
	 	for ( ; i<size; i++) {
		  Long pii = pi[i];
		  double *pM = MPT + pii * size;
		  p_RESULT[i] = 0.0;
		  if (FABS((p_rhs[pii] - SCALAR(pM, p_RESULT, i))) > eps) {
		    if (Pt == NULL) solve_DELETE(&pt);
#ifdef DO_PARALLEL
		    RFERROR("Equation system not solvable"); // RFERROR necessary.
#else		    
		    GERR1("Equation system not solvable (difference %10e). Try increasing 'pivot_relerror' in 'RFoptions' to get an approximate solution.",
			 p_rhs[pii] - SCALAR(pM, p_RESULT, i));
#endif
		    
		  }
		}
	      }

#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (rhs_cols > 30)
#endif 
	      for (Long k=0; k<rhs_cols; k++) {	      
		double *p_RESULT = RESULT + k * size;
		for (Long i=actual_size - 1; i>=0; i--) {
		  Long pii = pi[i];
		  double *pM = MPT + pii * size,
		    r = (p_RESULT[i] /= pM[i]);
		  LINEAR(pM, -r, i, p_RESULT);
		  // for (Long j=0; j<i; j++) p_RESULT[j] -= r * pM[j];
		}
	      }
	      CMALLOC(xja, size, int);
	      Sort(RESULT, size, rhs_cols, pi, xja, D);
	    }// rhs_cols > 0
	  } // not sqrt only

	} // err == NOERROR	
      }

      ERR_CHOL:
      
      if (err != NOERROR) {	
	if (pt->actual_pivot == PIVOT_NONE)
	  CERR2("Probably matrix not positive definite: %.300s. Consider 'RFoptions(%.80s)'.\n",
		ErrStr,
		calculate != SOLVE || rhs_cols > 0  ?
		"'pivot=PIVOT_AUTO)' or 'RFoptions(pivot=PIVOT_DO"
		: "solve_method=\"eigen\", solve.pseudoinverse=TRUE") // OK

	    
	  else // pt->actual_pivot == PIVOT_DO  or  PIVOT_IDX
	    CERR1("Likely, the matrix is not positive semi definite: %.300s. Consider 'RFoptions(solve_method=\"svn\"\n", ErrStr)
	      }
      
      if (PL >=  PL_DETAILSUSER) {
	PRINTF("Cholesky decomposition successful\n");
      }
      }
      break;
    case QR : {// QR returns transposed of the inverse !! 
      if (rhs_cols > 0 || logdet != NULL || calculate != SOLVE) {
	err = ERRORFAILED;
	continue;
      }

        err = ERRORNOTPROGRAMMEDYET; /// to do: clarify transposed !
      continue;

      CMALLOC(w2, size, double);
      CMALLOC(w3, size, double);

      F77dgeqrf(&size, &size, // QR
		       MPT, &size, // aijmax, &irank, inc, w2, 
		       w3, w2, &size, &err);     
      assert(false); // code zur invertierung fehlt hier!
      
      if (err != NOERROR) {	
	CERR1("'dgeqrf' failed with err=%d.", err);
      }
      if (PL >=  PL_DETAILSUSER) { PRINTF("QR successful\n"); }
      break;
    }
      
    case Eigen : { //  M = U D UT
      int max_eigen = sp->max_svd;
      double eigen2zero = sp->eigen2zero;
      if (size > max_eigen) CERR0("matrix too large for Cholesky or eigen value decomposition. Increase 'max_chol' and 'max_svd' in 'RFoption' if necessary.");
      
      double 
	optimal_work,
	*pt_work = &optimal_work;
      int k=0,
	optimal_intwork,
	*pt_iwork = &optimal_intwork, 
	lw2 = -1,
	lintwork = -1;
      
      C_MALLOC(Long, U, sizeSq, double);
      CMALLOC(D, size, double); 
      CMALLOC(xja, 2 * (Long) size, int);
      CMALLOC(w3, size, double);
 
      for (int i=0; i<=1; i++) {
	double dummy = 0.0,
	  abstol = 0.0;
	int dummy_nr;

	//	printf("%f %f %f %f\n", MPT[0],  MPT[1],  MPT[2],  MPT[3]);
	F77dsyevr("V", "A", "U", &size, // Eigen
			 MPT, &size, &dummy, &dummy, &k, &k, 
			 &abstol,// or DLAMCH
			 &dummy_nr, D, U, &size, 
			 xja, // 2 * size * size of(integer); nonzeros_idx
			 pt_work, &lw2,
			 pt_iwork, &lintwork,
			 &err
#ifdef USE_FC_LEN_T	  
		  FCONE FCONE FCONE
#endif		  
			 );
	//	printf("eigen %d %f %f;  %e %e %e %e\n", err, D[0], D[1], U[0], U[1], U[2], U[3]);
	if (i==1 || err != NOERROR || ISNAN(D[0])) break;
	lw2 = (int) optimal_work;
	lintwork = (int) optimal_intwork;
	CMALLOC(w2, lw2, double);
	CMALLOC(iwork, lintwork, int);
	pt_iwork = iwork;
	pt_work = w2; 
      }

      
      if (err != NOERROR) {
	if (PL>PL_ERRORS) { PRINTF("F77 error code for 'dsyevr' = %d\n", err);}
	CERR1("'dsyevr' failed with err=%d.", err);
	break;
      }
      
      for (Long i=0; i<size; i++)
	if (D[i] < -eigen2zero) {
	  const char
	    *advice[2]={"", " Consider increasing the value of 'eigen2value'."};
	  double min  = D[i];
	  for (Long j=i+1; j<size; j++) if (D[j] < min) min = D[j];
	  GERR3("Negative eigenvalues found (less than -eigen2zero=%10e). Smallest one equals %10e. %.50s", -eigen2zero, min, advice[min > -eigen2zero * 100]);
	  
	} 
      
      if (calculate == MATRIXSQRT) {
	for (Long j=0; j<size; j++) {
	  double dummy;
	  dummy = 0.0;
	  if (D[j] >= eigen2zero) dummy = SQRT(D[j]);
	  for (Long i=0; i<size; i++, k++) RESULT[k] = U[k] * dummy;
	}
      } else {
	// calculate determinant 
	if (logdet != NULL) {
	  *logdet = cumProd(D, size, sp->det_as_log);
	  if (calculate == DETERMINANT) {
	    err = NOERROR;
	    goto ErrorHandling;
	  }
	}

	//	printf("Hiere EIGGEN\n");
	bool pseudoInverse = false;
	for (Long j=0; j<size; j++) {
	  pseudoInverse |= D[j] < eigen2zero;
	  w3[j] = D[j] < eigen2zero ? 0.0 : 1.0 / D[j];
	}
	
	if (rhs_cols > 0) {
	  Long tot = (Long) size * rhs_cols;
	  C_MALLOC(Long, w2, tot, double);	
	  matmulttransposed(U, RHS, w2, size, size, rhs_cols, cores);
	  if (pseudoInverse) {
	    for (k=0; k<tot; ) {
	      for (Long i=0; i<size; i++) {
		double cmp = w2[k];
		w2[k] *= w3[i];
		if (FABS(w2[k] * D[i] - cmp) > eigen2zero) {
		  GERR0("singular matrix problem does not have a solution");
		}
		k++;
	      }
	    }
	  } else {
	    for (k=0; k<tot; )
	      for (Long i=0; i<size; i++) w2[k++] *= w3[i];
	  }
	  matmult(U, w2, RESULT, size, size, rhs_cols,cores);
	} else {
	  if (pseudoInverse && !sp->pseudoinverse)
	    GERR0("Singular matrix: inverse does not exist. Consider 'RFoption(solve.pseudoinverse=TRUE)'");
	  C_MALLOC(Long, w2, sizeSq, double);
	  for (Long kk=0, j=0; j<size; j++) {
	    double dummy = w3[j];
	    for (Long i=0; i<size; i++, kk++) w2[kk] = U[kk] * dummy;
	  }
	  matmult_2ndtransp(w2, U, RESULT, size, size, size); // V * U^T
	}
      }

      if (PL >=  PL_DETAILSUSER) {
	PRINTF("eigen value decomposition successful\n");
      }
      break;
    }
	
    case SVD : {// SVD : M = U D VT
      if (size > sp->max_svd) CERR0("matrix too large for SVD decomposition.");
      int 
	lw2 = -1,
	size8 = size * 8;
      double  optim_lwork,
	eigen2zero = sp->eigen2zero,
	*pt_w2 = &optim_lwork;
      
      C_MALLOC(Long, w3, sizeSq, double);
      C_MALLOC(Long, U, sizeSq, double);
      CMALLOC(D, size, double); 
      CMALLOC(iwork, size8, int);
      CMALLOC(lnz, size, double);
       
      for (Long i=0; i<=1; i++) {
	F77dgesdd("A", &size, &size, // SVD
			 MPT, &size, D, U, &size, w3, &size, 
			 pt_w2, &lw2, iwork, &err
#ifdef USE_FC_LEN_T	  
		  FCONE
#endif		  
		  );
	if (i==1 || err != NOERROR || ISNAN(D[0])) break;
	lw2 = (int) optim_lwork;
	CMALLOC(w2, lw2, double);
	pt_w2 = w2;
      }
      if (err != NOERROR) {
	if (PL>PL_ERRORS) { PRINTF("F77 error code for 'dgesdd' = %d\n", err);}
	CERR1("'dgesdd' failed with err=%d.", err);
	break;
      }

      if (calculate == MATRIXSQRT) {
	double svdtol = sp->svd_tol;
	/* calculate SQRT of covariance matrix */
 	for (Long j=0, k=0; j<size; j++) {
	  double dummy;
	  if (D[j] < -eigen2zero) CERR0("negative eigenvalues found.");
	  dummy = 0.0;
	  if (D[j] >= eigen2zero) dummy = SQRT(D[j]);
	  for (Long i=0; i<size; i++, k++) RESULT[k] = U[k] * dummy;
	}
 
	/* check SVD */
 	if (svdtol > 0.0) {
	  for (Long i=0; i<size; i++) {
	    double *Ui = RESULT + i;
	    for (Long k=i; k<size; k++) {
	      double *Uk = RESULT + k,
		sum = 0.0;
	      for (Long j=0; j<sizeSq; j+=size) {
		sum += Ui[j] * Uk[j];
	      }
	      
	      if (FABS(M0[i * size + k] - sum) > svdtol) {
		if (PL > PL_ERRORS) {
		  PRINTF("difference %10e at (%ld,%ld) between the value (%10e) of the covariance matrix and the square of its root (%10e).\n", 
			 M0[i * size +k] - sum, i, k, M0[i*size+k], sum);
		}
		FERR3("required precision not attained  (%10e > %10e): probably invalid model. See also '%.50s'.", FABS(M0[i * size + k] - sum), svdtol,
		      solve[SOLVE_SVD_TOL]);

		err=ERRORM;
		break;
	      } //else printf("ok (%d,%d) %10g %10g\n", i, k, M0[i*size+k],sum);
	    }
	    if (err != NOERROR) break;		
	  }
	  if (err != NOERROR) break;		
	} // end if svdtol > 0

      } else {
        // calculate determinant 
        if (logdet != NULL) {
	  *logdet = cumProd(D, size, sp->det_as_log);
	  if (calculate == DETERMINANT) {
	    err = NOERROR;
	    goto ErrorHandling;
	  }
	}
 	
  	bool pseudoInverse = false;
	for (Long j=0; j<size; j++) {
	  bool small = FABS(D[j]) < eigen2zero;
	  pseudoInverse |= small;
	  lnz[j] = small ? 0.0 : 1.0 / D[j];
	}
	
	if (rhs_cols > 0) {
	  Long tot = (Long) size * rhs_cols;
	  C_MALLOC(Long, w2, tot, double);	
	  matmulttransposed(U, RHS, w2, size, size, rhs_cols,cores);
	  if (pseudoInverse) {
	    for (Long k=0; k<tot; ) {
	      for (Long i=0; i<size; i++) {
		double cmp = w2[k];
		w2[k] *= lnz[i];
		if (FABS(w2[k] * D[i] - cmp) > eigen2zero) {
		  GERR0("singular matrix problem does not have a solution.");
		}
		k++;
	      }
	    }
	  } else {
	  for (Long k=0; k<tot; )
	    for (Long i=0; i<size; i++) w2[k++] *= lnz[i];
	  }
	  matmulttransposed(w3, w2, RESULT, size, size, rhs_cols,cores);
	} else {
	  // calculate inverse of covariance matrix
	  if (pseudoInverse && !sp->pseudoinverse)
	    GERR0("Singular matrix: inverse does not exist. Consider 'RFoption(solve.pseudoinverse=TRUE)'");
	  for (Long k=0, j=0; j<size; j++) {
	    double dummy = lnz[j];
	    for (Long i=0; i<size; i++) U[k++] *= dummy;
	  }
	  matmult_tt(U, w3, RESULT, size, size, size, cores); // V * U^T
	}
      }

      if (PL >=  PL_DETAILSUSER) { PRINTF("svd successful\n"); }
      break;
    }

    case LU : {// LU
      //printf("LU\n");
      if (calculate == MATRIXSQRT) {
	err = ERRORFAILED;
	continue;
      }
      
      CMALLOC(xja, size, int);		    
      F77dgetrf(&size, &size, MPT, &size, xja, &err);
      if (err != NOERROR) {
	CERR1("'dgetrf' (LU) failed with err=%d.", err);
      }
 
      //printf("LU %d\n", logdet != NULL);

      
      if (logdet != NULL) {
	*logdet =  DeterminantLU(MPT, size, sp->det_as_log, xja);
        if (calculate == DETERMINANT) {
	  err = NOERROR;
	  goto ErrorHandling;
	}
      }

      if (rhs_cols > 0) {
	if (rhs_cols > MAXINT) BUG;
	int rhs_cols0 =  (int) rhs_cols;
	Long totalRHS = (Long) size * rhs_cols;
	if (result != NULL) MEMCOPY(RESULT, RHS, sizeof(*RESULT) * totalRHS);
	F77dgetrs("N", &size, // LU rhs
			 &rhs_cols0, MPT, &size, xja, 
			 RESULT, &size, &err
#ifdef USE_FC_LEN_T		  
		  FCONE
#endif		  
		  );
	if (err != NOERROR) {	
	  CERR1("'dgetrs' (LU) failed with err=%d.", err);
	}
      } else {
	int lw2 = -1;
	double dummy,
	  *p = &dummy;
	for (int i=0; i<=1; i++) { 
	  F77dgetri(&size, MPT, // LU solve
			   &size, xja, p, &lw2, &err);	
	  if (err != NOERROR) break;
	  lw2 = (int) dummy;
	  CMALLOC(w2, lw2, double);
	  p = w2;
	}
      }      
      if (PL >=  PL_DETAILSUSER) { PRINTF("LU decomposition successful\n"); }
      break;
    }
	
    case Sparse : {// sparse matrix
#if defined compatibility_to_C_h
      BUG;
#else
      if (sizeSq > MAXINT) BUG;
      
      Long halfsq = (Long) size * (size + 1) / 2;
      int nnzlindx = -1, 
	doperm = sp->pivotsparse,
	nnzcolindices = 0,
	nnzR = 0,
	cache = 512, // to do: CPU cache size
	nnzcfact[3] = { 5, 1, 5 }, 
	nnzRfact[3] = { 5, 1, 2 };
      double
	cholincrease_nnzcol = 1.25,
	cholincrease_nnzR = 1.25;

      if (!posdef) CERR0("'spam' needs a positive definite matrix.");
      CMALLOC(pivotsparse, size, int);
      if (!doperm) for (int i=0; i<size; i++) pivotsparse[i] = i + 1;

      if (spam_zaehler == 0) { 
	for (Long i=0; i<sizeSq; i++) nnzA += !(FABS(M0[i]) <= spam_tol);
	spam_zaehler = nnzA + 1; // falls nur aus Nullen bestehend
      }

      CMALLOC(xlnz, sizeP1, int);
      CMALLOC(snode, size, int);
      CMALLOC(xsuper, sizeP1, int);
      CMALLOC(iwork, sizeP1, int);
      CMALLOC(invp, size, int);
      CMALLOC(w3, size, double);

      CMALLOC(cols, spam_zaehler, int);
      CMALLOC(rows, sizeP1, int);
   
      int nD = spam_zaehler;
      if (nD < size) nD = size;
      CMALLOC(D, nD, double);
      // prepare spam

      F77spamdnscsr(&size, &size, M0, &size, D,
			   cols, // ja
			   rows, // ia
			   &spam_tol); // create spam object   
      pt->nsuper = 0;
      // calculate spam_cholesky
      err = 4; // to get into the while loop
      while (err == 4 || err == 5) {
	if (nnzcolindices == 0) {
	  double rel = nnzA / (double) size;
	  if (rel < 5) {
	    nnzcolindices = (int) CEIL(nnzA * (1.05 * rel - 3.8));
	    if (nnzcolindices < 1000) nnzcolindices = 1000;
	  } else {
	    nnzcolindices = nnzA;
	  }
	  nnzcolindices *= nnzcfact[doperm];
	  if (nnzcolindices < nnzA) nnzcolindices = nnzA;
	} else if (err == 5) {
	  int tmp = (int) CEIL(nnzcolindices * cholincrease_nnzcol);
	  if (PL > PL_RECURSIVE) {
	    PRINTF("Increased 'nnzcolindices' with 'NgPeyton' method\n(currently set to %d from %d)", tmp, nnzR);
	  }
	  nnzcolindices = tmp;
	}
	if (nnzcolindices < pt->n_lindx) nnzcolindices = pt->n_lindx;
	
	if (nnzR == 0) {
	  double u = FLOOR(.4 * POW(nnzA, 1.2));
	  u = u < 4 * nnzA ? 4 * nnzA : CEIL(u);
	  nnzR = (int) u * nnzRfact[doperm];
	} else if (err == 4) {
	  int tmp = (int) CEIL(nnzR * cholincrease_nnzR);
	  if (PL > PL_RECURSIVE) {
	    PRINTF("Increased 'nnzR' with 'NgPeyton' method\n(currently set to %d from %d)", tmp, nnzR);
	  }
	  nnzR = tmp;
	}
	if (nnzR < pt->n_lnz) nnzR = pt->n_lnz;
	else if (nnzR > halfsq) nnzR = (int) halfsq;	
	
	CMALLOC(lindx, nnzcolindices, int);	
	CMALLOC(lnz, nnzR, double);
	 	
	F77cholstepwise(&size, &nnzA, D, cols, rows, &doperm,
			       invp, pivotsparse, 
			       &nnzlindx, &nnzcolindices, 
			       lindx, // 
			       iwork,// 
			&(pt->nsuper), // length of lindx
			       &nnzR,  // physical length of lindx
			       lnz,   // output:result
			       xlnz,  // cols of lnz "ja"
			       snode,  // supernode membership ??
			       xsuper, // supernode partioning
			&cache, // cache size of the CPU
			&err
			);       
	
	if (err != NOERROR) {
	  CERR1("'cholstepwise' failed with err=%d.", err);
	  break;
	}	 
      } // while
      
      if (err != NOERROR) CERR0("'spam' failed.");
      if (PL >=  PL_DETAILSUSER) { PRINTF("'spam' successful\n"); }
      
      // spam solve
      
      if (calculate == MATRIXSQRT) {
	
	//BUG; // unexpected behaviour in spam
	
	nnzR = xlnz[size] - 1;
	CMALLOC(xja, nnzR, int);
	F77calcja(&size, &(pt->nsuper), pt->xsuper, 
		  pt->lindx, pt->iwork, pt->xlnz, xja);
	for (Long i=0; i<size; invp[i++]--); 
	F77spamcsrdns(&size, pt->lnz, xja, pt->xlnz, RESULT);
	for (Long i=0; i<size; i++) {
	  Long endfor = (i + 1) * size;
	  for (Long j = i * (size + 1) + 1; j<endfor; RESULT[j++]=0.0);
	}
      } else {     
	double *lnz = pt->lnz;
	int *lindx = pt->lindx;
	
	// spam determinant
	if (logdet != NULL) {
	  if (sp->det_as_log) {
	    double tmp = 0.0;
	    for (Long i=0; i<size; i++) tmp += LOG(lnz[xlnz[i] - 1]);	  
	    *logdet = 2.0 * tmp;
	  } else {
	    double tmp = 1.0;
	    for (Long i=0; i<size; i++) tmp *= lnz[xlnz[i] - 1] ;
	    *logdet = tmp * tmp;
	  }
	  if (calculate == DETERMINANT) {
	    err = NOERROR;
	    goto ErrorHandling;
	  }
	}
	
	/*   z = .Fortran("backsolves", m = nrow,
	     nsuper, p, a@colindices,
	     a@colpointers, as.double(a@entries),
	     a@rowpointers, a@invpivot, a@pivot,
	     a@supernodes, vector("double",nrow),
	     sol = vector("double",nrow*p),
	     as.vector(b,"double"),
	     NAOK = .Spam$NAOK,PACKAGE = "spam")$sol
	*/
	int RHS_COLS;
	if (rhs_cols <= 0) { // UNBEDINGT VOR double *RHS;
	  RHS_COLS = size;	
	  FILL_IN(RESULT, sizeSq, 0.0);
	  for (Long i=0; i<sizeSq; i += sizeP1) RESULT[i] = 1.0; 
	} else {
	  if (rhs_cols > MAXINT) BUG;
	  RHS_COLS = (int) rhs_cols;	
	  if (result != NULL) 
	    MEMCOPY(RESULT, RHS, (Long) size * rhs_cols * sizeof(*RESULT));
	}
	
	//printf("nsuper=%d\n", pt->nsuper);
	//	  for (Long ii=0; ii<size; ii++) 
	//printf("%d %d %d %d %10e\n", ii, pt->nsuper, sizeP1, xsuper[ii],
	//	   w3[ii]);
	
	//	  if (false)
	//	  for (Long jsub=0; jsub<=pt->nsuper; jsub++) {
	//	    int fj = xsuper[1 - 1],
	//	      Lj = xsuper[jsub + 1 - 1] -1;
	//	    printf("%d %d %d\n", jsub, fj, Lj);
	//	    for (Long jcol=fj; jcol <= Lj; jcol++) {
	//	      printf("%d,%10g  ", jcol, w3[jcol -  1]);
	//	    }	    
	//	  }
	
	//	  for (Long jcol=1; jcol <= 600; jcol++) {
	//	    w3[jcol - 1] = jcol;
	//   printf("%d,%10g  ", jcol, w3[jcol -  1]);
	//  }	    
	
	
	//	  printf("%ld %ld %d\n", RESULT, rhs, rhs_cols);
	//	  for (Long ii=0; ii<size; ii++) printf("%d %10e\n", ii, RESULT[ii]);
	//	  BUG;
	
	F77backsolves(&size, &(pt->nsuper), &RHS_COLS, 
			     lindx, // colindices
			     iwork, //colpointers
			     lnz, 
			     xlnz, //  rowpointers
			     invp, pivotsparse,
			     xsuper, // supernodes
			     w3, RESULT);	
	if (PL >=  PL_DETAILSUSER) { PRINTF("'spam' successful\n"); }
      }         
      break;
#endif 
   } // Sparse

    case NoInversionMethod: GERR0("no inversion method given.");
    case NoFurtherInversionMethod:
      STRCPY(ErrStr, WHICH_ERRORSTRING);
      GERR1("%.300s (All specified matrix inversion methods have failed.)",
	    ErrStr);
    case direct_formula: case Diagonal: 
      GERR1("strange method appeared:%.200s.", CONTACT);
    default : GERR1("unknown method (%d) in 'RandomFieldsUtils'.", pt->method);
    } // switch

    if (err==NOERROR) break;
  } // for m

 
 ErrorHandling:
  if (Pt == NULL) solve_DELETE(&pt);      
  else Pt->sparse = sparse;
  
  //  if (err != NOERROR) { printf("Err = %s %s %d\n",  ErrStr, WHICH_ERRORSTRING ,err); exit(0);}
  
  return err; //  -method;
}

 


int SolvePosDef(double *M, int size, bool posdef, 
		double *rhs, Long rhs_cols, double *logdet, 
		solve_storage *PT,
		int VARIABLE_IS_NOT_USED cores) {
  if ((rhs == NULL) xor (rhs_cols == 0)) BUG; 
  return doPosDefIntern(M, size, posdef, rhs, rhs_cols,
		  NULL, // result, so result returned in M or rhs
		  logdet, 
		  SOLVE, // calculate
		  PT, // storage
		  &(OPTIONS.solve), cores);
}


int SolvePosDefSp(double *M, int size, bool posdef, 
		  double *rhs, Long rhs_cols,  double *logdet, 
		  solve_storage *PT,  solve_options * sp, int VARIABLE_IS_NOT_USED cores) {  
  if ((rhs == NULL) xor (rhs_cols == 0)) BUG; 
  return doPosDefIntern(M, size, posdef, rhs, rhs_cols, NULL, logdet, SOLVE,
		  PT, sp, cores);
}




int xCinvYdet(double *M, int size, bool posdef, double *X, double *Y,
	      Long cols, double *XCinvY, double *det, bool log,
	      solve_storage *PT, int VARIABLE_IS_NOT_USED cores) {
  // called by randomfields
  int NR = KAHAN ? SCALAR_KAHAN : SCALAR_AVX;
  bool pt = PT != NULL && PT->result != NULL;
  double *result;
  if (pt) result=PT->result;
  else result= (double *) MALLOC(sizeof(*result) * (Long) size * cols);  
  if (result == NULL) return ERRORMEMORYALLOCATION;
  double *res = result;
  solve_options sp;
  MEMCOPY(&sp, &(OPTIONS.solve), sizeof(solve_options));
  sp.det_as_log = log;
  int err =  doPosDefIntern(M,// no PROTECT( needed
		      size, posdef, Y, cols, result, det, SOLVE, PT, &sp,
		      cores);
  for (Long i=0; i<cols; i++, res += size, X += size)
    XCinvY[i] = SCALAR(res, X, size);
  if (!pt) { FREE(result); }
  return err;
}

int xCinvXdet(double *M, int size, double *X, Long X_cols,
	      double *XCinvX, double *det, bool log, solve_storage *PT,
	      int VARIABLE_IS_NOT_USED cores) {
  // called by randomfields
  return xCinvYdet(M, size, true, X, X, X_cols, XCinvX, det, log, PT, cores);
}


double DetPosDefsp(double *M, int size, solve_options *sp, int VARIABLE_IS_NOT_USED cores) {
  double det;
   int err= doPosDefIntern(M, size, true, NULL, 0, NULL, &det,// no PROTECT( needed
		     DETERMINANT, NULL, sp, cores);
  if (err != NOERROR)
    ERR0("error occurred when calculating determinant of a pos def matrix.");
  return det;  
}

double DetPosDef(double *M, int size, int VARIABLE_IS_NOT_USED cores) {
  //never log of det --always small sizes!!
  //not that M will be destroyed!!
  solve_options sp;
  MEMCOPY(&sp, &(OPTIONS.solve), sizeof(solve_options));
  sp.det_as_log = false;
  return DetPosDefsp(M, size, &sp, cores);
}


int InvertMatrix(double *M, int size, int VARIABLE_IS_NOT_USED cores) {
  return doPosDefIntern(M, size, false, NULL, 0, NULL, NULL, SOLVE, NULL, NULL,
		  cores);
}




int chol(double *M, int size, solve_options *sp, int VARIABLE_IS_NOT_USED cores) {
  return doPosDefIntern(M, size, true, NULL, 0, NULL, NULL, MATRIXSQRT, NULL, sp,
		  cores);   
}

int cholesky(double *M, int size, int VARIABLE_IS_NOT_USED cores) {
  solve_options sp;
  MEMCOPY(&sp, &(OPTIONS.solve), sizeof(solve_options));
  sp.Methods[0] = sp.Methods[1] = Cholesky;
  sp.sparse = False; // currently does not work, waiting for Reinhard
  return chol(M, size, &sp, cores);   
}


bool Is_positive_definite(double *C, int dim, int VARIABLE_IS_NOT_USED cores) { // bool not allowed in C
  int err;
  Long bytes = sizeof(*C) * dim * dim;
  double *test;
  test = (double*) MALLOC(bytes);
  MEMCOPY(test, C, bytes);
  err = cholesky(test, dim, cores);
  //  printf("errr = %d\n", err);
  UNCONDFREE(test);
  return(err == NOERROR);
}



/*

  ## extrem wichter check -- folgendes funktioniert bislang bei spam nicht:
  library(RandomFields, lib="~/TMP")
  RFoptions(printlevel = 3, pch="", seed=999, use_spam = TRUE) #//
  z = RFsimulate(RMspheric(), x, max_variab=10000, n=10000, spC=F ALSE)
  C = cov(t(z))
  c = RFcovmatrix(RMspheric(), x) #//
  p rint(summary(as.double(c - C))) ##//
  stopifnot(max(a b s(c-C)) < 0.05)

*/


int SqrtPosDefFree(double *M,  // in out
		   int size,  
		   solve_storage *pt,     // in out
		   solve_options *sp,
		   int VARIABLE_IS_NOT_USED cores
		   ){
  ASSERT_SOLVE(sp);
  int err;
  Long sizeSq = (Long) size * size;
  InversionMethod *Meth = sp->Methods;
  double *res = NULL;
  bool extra_alloc = Meth[0] == NoInversionMethod ||
    Meth[0] == NoFurtherInversionMethod ||
    (Meth[1] != NoInversionMethod && Meth[1] != NoFurtherInversionMethod &&
     Meth[1] != Meth[0]) ||
    (Meth[0] != Cholesky && Meth[0] != Eigen && Meth[0] != SVD);
  assert(pt != NULL);

  if (sp->sparse == True) 
    warning("package 'spam' is currently not used for simulation");
  sp->sparse = False;

  if (extra_alloc) {
    C_MALLOC(Long, result, sizeSq, double);
    res = result;
  } else {
    FREE(pt->result);
    pt->result = M;
    pt->n_result = sizeSq;
  }
  
  // it is ok to have
  // ==15170== Syscall param sched_setaffinity(mask) points to unaddressable byte(s)
  // caused by gcc stuff
  err = doPosDefIntern(M, size, true, NULL, 0, res, NULL, MATRIXSQRT, pt, sp,
		 cores);// no PROTECT( needed
  
  if (extra_alloc) {  


#if defined MSDOS_WINDOWS
    pt->to_be_deleted = M;
#else
      FREE(M);
#endif
  }
  return err;
}


void sqrtRHS_Chol(double *U, int size, double* RHS, Long RHS_size, Long n,
		  double *result,
		  bool pivot, int act_size,
		  int cores,
		  int *pi) {
  //  printf("n=%d,rhss=%d si=%d pivot=%d, act=%d U=%d  RHS=%d %d pi=%d\n",
  //	 n,  RHS_size, size,pivot,  act_size, U!=NULL, RHS!=NULL,   result!=NULL, pi!=NULL );
//  for (Long i=0; i<size; i++) printf("%d ", pi[i]); printf("\n");

  Long nsize = n * (Long) size;
  int NR = KAHAN ? SCALAR_KAHAN : SCALAR_AVX;
  assert(U != NULL);
  if (pivot){
    Long n_act_size = n * (Long) act_size,
      diff = size - act_size,
      n_diff = nsize - n_act_size;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (MULTIMINSIZE(n_act_size)) schedule(dynamic, 8) 
#endif
    for (Long i=0 ; i<n_act_size; i++) {
      Long ii = i % act_size,
	j = i / act_size;
      result[pi[ii] + j * size] = SCALAR(RHS + j * RHS_size,
					 U + pi[ii] * (Long) size, ii + 1);
    }
    pi += act_size;
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (MULTIMINSIZE(n_diff))
#endif
    for (Long i=0; i<n_diff; i++) {
      Long ii = pi[i % diff],
	j = i / diff;
      result[ii + j * size] = SCALAR(RHS + j * RHS_size, U + ii * size,
				     act_size);
    }
  } else {
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) if (MULTIMINSIZE(nsize)) schedule(dynamic, 8) 
#endif
    for (Long i=0; i<nsize; i++) {
      Long ii = i % size,
	j = (i / size) * size;
      result[ii + j] = SCALAR(RHS + j, U + ii * size, ii + 1);
    }
  }
}


int sqrtRHS(solve_storage *pt, double* RHS, double *result, int cores){
  // Multipliziert die Cholesky_Zerlegung auf einen Vektor
  // doppelte Geschwindigkeit der Dreieckmatrix wird hier genutzt
  // (und dass der andere Teil der Matrix nicht wohl definiert ist)
  assert(pt != NULL);
  int size = pt->size;
  switch (pt->method) {
  case Rcholesky : {
    //printf("RchoL\n");
    int incx = 1;
    MEMCOPY(result, RHS, (Long) size * sizeof(*result));
    F77dtrmv("U", "T", "N", &size, pt->result, &size, result, &incx
#ifdef USE_FC_LEN_T	  
		  FCONE FCONE FCONE
#endif		  
	     );
  }
    break;
  case GPUcholesky :
    // Alex: Gegebenenfalls schnelle GPU Version von Deiner Seite
    
  case direct_formula : 
  case Cholesky : {
    //    printf("intern\n");
      bool pivot = (pt->actual_pivot == PIVOT_DO ||
		  pt->actual_pivot == PIVOT_IDX) &&
      pt->method != direct_formula;
    if (pivot && pt->n_pivot_idx != size) BUG;
    
    sqrtRHS_Chol(pt->result, size, RHS, size, 1, result, pivot,
		 pivot ? pt->actual_size : NA_INTEGER, cores,
		 pt->pivot_idx);
    return NOERROR;
  }
  case SVD : case Eigen : {  
    double *U = pt->result;
    assert(U != NULL);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(cores) schedule(static) if (MULTIMINSIZE(size))
#endif
    for (Long i=0; i<size; i++){
      double dummy = 0.0;
      Long k = i;
      for (Long j=0; j<size; j++, k+=size) dummy += U[k] * RHS[j];
      result[i] = (double) dummy; 
    }
  }
    break;

  case Sparse : {
    BUG; // SEE ALSO solve, calculate, tmp_delete !!
    int one = 1;
    assert(pt->D != NULL);
    F77amuxmat(&size, &size, &one, RHS, pt->D, pt->lnz, 
		      pt->xja, pt->xlnz);
    for (Long i=0; i<size; i++) result[i] = pt->D[pt->invp[i]];
  }
    break;

  case Diagonal : {  
    Long  i, j,
      sizeP1 = size + 1;
    double *D = pt->result;
    assert(D != NULL);
    for (i=j=0; j<size; j++, i+=sizeP1) result[j] = RHS[j] * D[i];
  }
    break;
    
  default :
    // printf("pt->method %d\n", pt->method);
    BUG;
  }
  
  return NOERROR;
}
