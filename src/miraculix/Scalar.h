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



#ifndef miraculix_ScalarInt_H
#define miraculix_ScalarInt_H 1

#define SCALAR_INT_8 0

#define SCALAR_INT_16 1
#define SCALAR_INT_AVX2 2

#if defined AVX2
//
#define SCALAR_INT_DEFAULT SCALAR_INT_AVX2
//#define SCALAR_INT_DEFAULT SCALAR_INT_16
#else
#define SCALAR_INT_DEFAULT SCALAR_INT_8
#endif

#define SCALARINT(x, y, len) scalarInt(x, y, len, SCALAR_INT_DEFAULT);
Long scalarInt(int *x, int *y, int len, int n);
Long scalarInt8by8(int *x, int *y, int len);
Long scalarIntAVX2(int * V1, int * V2, int N);

void matmulttransposedInt(int *A, int *B, double *c, int m, int l, int n,
			   int cores);
void crossprod_Int(int *x, int *y, int nrow, int ncol, int len,
		   int VARIABLE_IS_NOT_USED cores, int *ans);

#endif
