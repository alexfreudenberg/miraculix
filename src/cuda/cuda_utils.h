/*
 Authors 
 Alexander Freudenberg, alexander.freudenberg@stads.de

 Copyright (C) 2022-2023 Alexander Freudenberg

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

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cusolverSp.h>
#include <cusparse.h>


void debug_info(const char *s, ...);

int checkCuda();
int checkDevMemory(size_t required_mem);
int switchDevice();

int checkError(const char *func, int line, cudaError_t err);
int checkError(const char *func, int line, cublasStatus_t err);
int checkError(const char *func, int line, cusparseStatus_t err);
