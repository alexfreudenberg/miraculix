

#include <chrono>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <inttypes.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <unistd.h>
#include "dgemm_compressed_cuda.h"

int main() {
  int device_count;
  cudaGetDeviceCount(&device_count);
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasDestroy(handle);
  return 0;
}