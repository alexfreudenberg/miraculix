
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


#ifndef RFutils_public_H
#define RFutils_public_H 1

#ifdef __cplusplus
extern "C" {
#endif 
  SEXP scalarR(SEXP x, SEXP y, SEXP mode);
  SEXP struve(SEXP X, SEXP Nu, SEXP Factor_Sign, SEXP Expscaled);
  SEXP besselk_simd(SEXP X, SEXP Nu);
  SEXP I0ML0(SEXP X);
  SEXP gaussr(SEXP X, SEXP Derivative); 
  SEXP WMr(SEXP X, SEXP Nu, SEXP Derivative, SEXP Factor);
  SEXP logWM2r(SEXP X, SEXP Nu1, SEXP Nu2, SEXP Factor);

  SEXP SolvePosDefR(SEXP M, SEXP rhs, SEXP logdet);
  SEXP Chol(SEXP M);
  
  SEXP RFoptions(SEXP options);

  void loadoptionsRFU();
  void detachoptionsRFU();
  
  SEXP sortX(SEXP Data, SEXP From, SEXP To, SEXP NAlast);
  SEXP orderX(SEXP Data, SEXP From, SEXP To, SEXP NAlast);

  SEXP colMaxs(SEXP M);
  SEXP rowMeansX(SEXP M, SEXP Factor);
  SEXP rowProd(SEXP M);
  SEXP chol2mv(SEXP Chol, SEXP N);
  SEXP tcholRHS(SEXP C, SEXP RHS);
  SEXP DivByRow(SEXP M, SEXP V);
  SEXP quadratic(SEXP x, SEXP A);
  SEXP dbinorm(SEXP X, SEXP Sigma);
  SEXP dotXV(SEXP M, SEXP V);
  //  void Ordering(double *d, int *len, int *dim, int *pos);
  SEXP crossprodX(SEXP X, SEXP Y, SEXP mode);

  SEXP DebugCall();
  SEXP getPackagesToBeInstalled(SEXP Force);
  SEXP isGPUavailable();
  SEXP isNEONavailable();
  SEXP isX86_64();
  void setCPUs(int *n);
  void recompilationNeeded(int *n);
  SEXP SIMDmessages(SEXP pkgs);
  SEXP debuggingLevel();
  SEXP gpu_info(SEXP DEVICES);
  SEXP instruction_set(SEXP which,  SEXP pkgs, SEXP used);

  SEXP testStrassen(SEXP A, SEXP B, SEXP truenrowA, SEXP truenrowB);
  SEXP Update_utilsoption();
    
#ifdef __cplusplus
}
#endif



#endif 
 
