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


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <inttypes.h>
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

#include "cuda_utils.h"

#define STR(x) XSTR(x)
#define XSTR(x) #x

int get_print_level(){
  char *print_level_env = getenv("PRINT_LEVEL");
  if (print_level_env != NULL) {
    return atoi(print_level_env);
  }
  else {
    return 0;
  }
}

void debug_info(const char *s, ...) {
  if (get_print_level() > 0) {
    va_list argptr;
    va_start(argptr, s);
    printf("\033[36m\t ");
    vprintf(s, argptr);
    printf(" \033[37m\n");
    va_end(argptr);
  }
}

void print_compile_info(const char *message){
  if (get_print_level() >= 0) {
    printf("------------");
    printf("------------");
    printf("------------");
    printf("------------");
    printf("------------\n");
    printf("\tmiraculix - %s\n", message);
#if defined COMMIT_ID
    printf("Compiled on %s %s, git commit %s\n", __DATE__, __TIME__, COMMIT_ID);
#endif
    printf("------------");
    printf("------------");
    printf("------------");
    printf("------------");
    printf("------------\n");
  }
}

int checkError(const char *func, int line, cudaError_t err) {
  if (err != cudaSuccess) {
    printf("Internal error in function %s at line %d: %s\n", func, line,
           cudaGetErrorString(err));
    return 1;
  }
  return 0;
}

int checkError(const char *func, int line, cublasStatus_t err) {
  if (err != CUBLAS_STATUS_SUCCESS) {
    printf("Internal error in CUBLAS function %s at line %d: %s\n", func, line,
           cublasGetStatusString(err));
    return 1;
  }
  return 0;
}

int checkError(const char *func, int line, cusparseStatus_t err) {
  if (err != CUSPARSE_STATUS_SUCCESS) {
    printf("Error in call to cuSPARSE in function %s at line %d: %s\n", func, line,
           cusparseGetErrorString(err));
    return 1;
  }
  return 0;
}

int checkError(const char *func, int line, cusolverStatus_t err) {
  if (err != CUSOLVER_STATUS_SUCCESS) {
    printf("Error in call to cuSOLVER in function %s at line %d: %d\n", func, line,
           err);
    return 1;
  }
  return 0;
}


int checkCuda(){
  //
  // Driver check
  // The following section checks if a compatible driver and runtime is
  // installed
  //
  cudaError_t err;
  int driverVersion = 0, runtimeVersion = 0;

  err = cudaDriverGetVersion(&driverVersion);
  if (checkError(__func__, __LINE__, err) != 0)
    return 1;
  if (driverVersion == 0) { // Check if there's a CUDA driver on the system
    printf("No CUDA driver detected.");
    return 1;
  }

  // Check if the CUDA runtime is compatible 
  cudaRuntimeGetVersion(&runtimeVersion);
  int driverMajor = driverVersion / 1000,
      driverMinor = (driverVersion - driverMajor * 1000) / 10,
      runtimeMajor = runtimeVersion / 1000,
      runtimeMinor = (runtimeVersion - runtimeMajor * 1000) / 10;

  if (driverMajor <
      10) { // CUDA drivers below 10 don't have the required feature set
    printf("Your CUDA driver is of version %d.%d. CUDA driver versions below "
           "10.0 are not supported.\n",
           driverMajor, driverMinor);
    return 1;
  }
  if (driverVersion < runtimeVersion) { // Mismatches between runtime and driver
                                        // might lead to problems
    printf("This software has been compiled with CUDA version %d.%d, but "
           "your driver is of version %d.%d. If you run into errors, consider "
           "uprading your CUDA driver.\n",
           runtimeMajor, runtimeMinor, driverMajor, driverMinor);
  }

 return 0;
}

int checkDevMemory(size_t required_mem){
  //
  // Device memory check
  // Checks if the device has enough RAM available

  cudaError_t err;
  if (checkCuda() != 0)
      return 1;

  size_t free_mem = 0, total_mem = 0;
  err = cudaMemGetInfo(&free_mem, &total_mem);
  if (checkError(__func__, __LINE__, err) != 0)
      return 1;

  if (free_mem < required_mem) {
      printf("CUDA Error: Not enough memory available. \nRequired %zu GB, free "
            "%zu GB, total on device %zu GB \n ",
            required_mem / size_t(pow(1024, 3)), free_mem / size_t(pow(1024, 3)),
            total_mem / size_t(pow(1024, 3)));
      return 1;
  }

  return 0;
}

int switchDevice(){
  //
  // Select GPU device
  // The following section switches the current context to the requested device
  // 

  cudaError_t err;
  int  device            = 0;
  char *requested_device = getenv("CUDA_DEVICE");
  int  device_count      = 0;
  bool verbose           = get_print_level() >= 0;
  cudaDeviceProp prop;

  if (checkCuda() != 0){
    return 1;
  }
  cudaGetDeviceCount(&device_count);

  if (requested_device != NULL) {
      device = atoi(requested_device);
      if (verbose) {
        printf("Environment variable CUDA_DEVICE is set to %s.\n",
               requested_device);
      }
  } else {
      if (verbose) {
        printf("Environment variable CUDA_DEVICE is not set.\n");
      }
  }

  // Check if the requested device is available
  if (device >= device_count){
    printf("Device not available.\n");
    return -1;
  }

  if (verbose) {
      cudaGetDeviceProperties(&prop, device);
      printf("Using device %s (device no %d).\n", prop.name, device);
  }

  // Switch to device
  err = cudaSetDevice(device);
  if (checkError(__func__, __LINE__, err) != 0)
    return -1;

  return device;
}

int switchDevice(int device){
  //
  // Select GPU device
  // The following section switches the current context to the requested device
  // 

  cudaError_t err;
  bool verbose = get_print_level() >= 0;
  int  device_count      = 0;
  cudaDeviceProp prop;

  if (checkCuda() != 0) {
        return 1;
  }

  cudaGetDeviceCount(&device_count);

  // Check if the requested device is available
  if (device >= device_count){
    printf("Device not available.\n");
    return -1;
  }

  if (verbose) {
      cudaGetDeviceProperties(&prop, device);
      printf("Using device %s (device no %d).\n", prop.name, device);
  }
  
  // Switch to device
  err = cudaSetDevice(device);
  if (checkError(__func__, __LINE__, err) != 0)
    return -1;

  return device;
}


