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



#ifndef miraculix_gpusolve
#define miraculix_gpusolve 1

//#include "MX.h"

void gpuSolve(double *matrix, Long individuals, double *result);
void gpu_relmat_custom(Uint*, double*, Uint, Uint);
void gpu_relmat_cublas(Uint*, double*, Uint, Uint);

#define PADDIM 4L
#define THREADS_PER_BLOCK 1024 //2048 / 32

#endif
