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

#ifndef basic_rfutils_local_h
#define basic_rfutils_local_h 1


#include "def_rfu.h"
#include "intrinsics.h"
#include "Basic_RandomFieldsUtils.h"

#define F77dgeqrf F77call(dgeqrf)
#define F77dsyevr F77call(dsyevr)
#define F77dgetrf F77call(dgetrf)
#define F77dgetrs F77call(dgetrs)
#define F77dgetri F77call(dgetri)
#define F77dgesv F77call(dgesv)
#define F77dpotrf F77call(dpotrf)
#define F77dtrmv F77call(dtrmv)

#if defined compatibility_to_C_h
#define F77U
#else
#define F77U
#endif

#define F77nameUU(X,Y) F77name(X##Y)
#define F77callUU(X,Y) F77call(X##Y)
#define F77nameU(X,Y) F77nameUU(X,Y)
#define F77callU(X,Y) F77callUU(X,Y)

F77nameU(spamdnscsr,F77U)(int *nrow, int* ncol, double* dns, int* ndns, double* a, int* ja, int* ia, double* eps);//
#define F77spamdnscsr F77callU(spamdnscsr,F77U)

F77nameU(cholstepwise,F77U)(int*, int*, double* , int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*  , int*, int*, int*, int*, int*);
#define F77cholstepwise F77callU(cholstepwise,F77U)

F77nameU(calcja,F77U)(int*, int*, int*, int*, int*, int*, int*);
#define F77calcja F77callU(calcja,F77U)

F77nameU(spamcsrdns,F77U)(int*,  double *, int *, int*, double*  ); // ok
#define F77spamcsrdns F77callU(spamcsrdns,F77U)

F77nameU(backsolves,F77U)(int*, int*, int*, int*, int*, double*  , int*, int*, int*, int*, double*  , double*  );
#define F77backsolves F77callU(backsolves,F77U)

F77nameU(amuxmat,F77U)(int*, int*, int*, double*  , double*  , double*  , int*, int*);
#define F77amuxmat F77callU(amuxmat,F77U)

#endif
