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



#ifndef miraculix_matrix_H
#define miraculix_matrix_H 1

// #include "MX.h"

SEXP IsolveRelMat(Long individuals, double *Atau, double *tau, int ntau,
		  double *Vec, double* beta, int nbeta,
		  int returns, bool destroy,
		  int cores);
  
void matrixvector012I(double * m, int r, int c,  SEXP vector, double *a);
void vector012matrixI(SEXP vector, double * m, int r, int c,  double *a);

void VectorRelMatrix(SEXP SxI, SEXP SxI2, double *V, Long repetV, 
		     int cmp, 
		     option_type *global, utilsoption_type *utils,
		     double *ans);
void VectorCrossMatrix(SEXP SxI, double *V, Long repetV, 
		      option_type *global, utilsoption_type *utils,
		       double *ans);

void genoVector_raw_LongDouble(SEXP SxI, LongDouble *V, Long repetV, Long ldV, 
			   option_type *global, utilsoption_type *utils,
			   LongDouble *ans, Long ldAns);
void vectorGeno_raw_LongDouble(SEXP SxI, LongDouble *V, Long repetV, Long ldV, 
			   option_type *global, utilsoption_type *utils,
			   LongDouble *ans, Long ldAns);

void genoVector_means_double(SEXP SxI, double *V, Long repetV, Long ldV, 
			     option_type *global, utilsoption_type *utils,
			     double *ans, Long ldAns);
void genoVector_means_double(SEXP SxI, double *V, Long repetV,  
			     option_type *global, utilsoption_type *utils,
			     double *ans);

void vectorGeno_means_double(SEXP SxI, double *V, Long repetV, Long ldV, 
			     option_type *global, utilsoption_type *utils,
			     double *ans, Long ldAns);
void vectorGeno_means_double(SEXP SxI, double *V, Long repetV, 
			     option_type *global, utilsoption_type *utils,
			     double *ans);

gV_vG_header(Ulong,);
gV_vG_header(LongDouble,);
gV_vG_header(double,); 



void genoVector_means_LongDouble(SEXP SxI, LongDouble *V, Long repetV,Long ldV, 
			     option_type *global, utilsoption_type *utils,
			     LongDouble *ans, Long ldAns);
void genoVector_means_LongDouble(SEXP SxI, LongDouble *V, Long repetV,  
			     option_type *global, utilsoption_type *utils,
			     LongDouble *ans);

void vectorGeno_means_LongDouble(SEXP SxI, LongDouble *V, Long repetV,Long ldV, 
			     option_type *global, utilsoption_type *utils,
			     LongDouble *ans, Long ldAns);
void vectorGeno_means_LongDouble(SEXP SxI, LongDouble *V, Long repetV, 
			     option_type *global, utilsoption_type *utils,
			     LongDouble *ans);


#endif

