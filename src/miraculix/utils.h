
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



#ifndef rfutils_utils_H
#define rfutils_utils_H 1


void bes_k_simd (double *xv, double alpha, int sx, double *yv);
void set_num_threads();

void colMaxsI(double *M, Long r, Long c, double *ans, int cores);
void colMaxsIint(int *M, Long r, Long c, int *ans, int cores);
void rowProdI(double *M, Long r, Long c, double *ans);
void rowMeansI(void *M, int sexp, Long r, Long c, double *weight, double *ans);
void dbinormI(double *x, double *y, Long nrow, double *sigma, double *ans);
void dotXVI(double *M, Long r, Long c, double *V, double *Ans);
#if defined rfutils_options_H
void AtAInt(int *a, int *b,
	    Long nrow, Long ncol, // a and b have same size
	    Long ld, Long ldC,
	    Long *C, // result
	    int VARIABLE_IS_NOT_USED cores,
	    solve_options *options);


#endif

  
void sqrtRHS_Chol(double *U, int size, double* RHS, Long RHS_size,
  Long n, double *result, bool pivot, int act_size, int cores,
  int *pi);

#endif
