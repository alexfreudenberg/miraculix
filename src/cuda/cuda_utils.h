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


/**
 * @brief Retrieve the environment variable "PRINT_LEVEL", convert it to an integer and return it.
 * 
 * @return The integer value of "PRINT_LEVEL".
 */
int get_print_level();

/**
 * @brief Print debugging information, if PRINT_LEVEL is greater than 0.
 *
 * @param s The debug message to print.
 * @param ... Additional arguments to format the debug message.
 */
void debug_info(const char *s, ...);

/**
 * @brief Print information about the compile time of the library.
 *
 * @param message The compile information to print.
 */
void print_compile_info(const char *message);


/**
 * @brief Check if a compatible CUDA driver and runtime is installed.
 * 
 * @return 0 if compatible driver and runtime is installed, otherwise returns an error code.
 */
int checkCuda();

/**
 * @brief Check if the GPU device has enough RAM available.
 *
 * @param required_mem The amount of RAM required, in bytes.
 * @return 0 if there is enough RAM available, otherwise returns an error code.
 */
int checkDevMemory(size_t required_mem);

/**
 * @brief Switch the current context to the device requested in the CUDA_DEVICE environment variable.
 * 
 * @return 0 if the context was successfully switched, otherwise returns an error code.
 */
int switchDevice();



int checkError(const char *func, int line, cudaError_t err);
int checkError(const char *func, int line, cublasStatus_t err);
int checkError(const char *func, int line, cusparseStatus_t err);
int checkError(const char *func, int line, cusolverStatus_t err);
