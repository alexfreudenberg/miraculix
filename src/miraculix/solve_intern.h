
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


#ifndef RFutils_solve
#define RFutils_solve 1


int doPosDefIntern(double *M0, int size, bool posdef,
  double *rhs0, Long rhs_cols, double *result,
  double *logdet, int calculate, solve_storage *Pt,
  solve_options *sp, int cores);


void sqrtRHS_Chol(double *U, int size, double* RHS, Long RHS_size,
  Long n, double *result, bool pivot, int act_size, int cores,
  int *pi);

#endif
