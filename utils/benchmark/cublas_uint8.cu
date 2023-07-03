/*
 Alexander Freudenberg, alexander.freudenberg@stads.de

 Copyright (C) 2020-2023 Alexander Freudenberg

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

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "cuda_utils.h"


extern "C" {
int cublas_uint8_gemm(unsigned char *snp_matrix, int snps, int indiv,
                       double *ans) {
    /*
    xx
    */
    cublasStatus_t status = CUBLAS_STATUS_SUCCESS;
    cudaError_t err = cudaSuccess;
    cublasHandle_t handle;

    cublasGemmAlgo_t    algo         = CUBLAS_GEMM_DEFAULT;
    cudaDataType_t      input_type   = CUDA_R_8I;
    cudaDataType_t      output_type  = CUDA_R_32F;
    cublasComputeType_t compute_type = CUBLAS_COMPUTE_32F;

    void  *d_A = NULL,
          *d_B = NULL;
    float *d_C = NULL,
          *h_C = NULL;

    int nrowA = indiv,
        ncolA = snps;

    const float alpha = 1.0,
          beta        = 0.0;
    
    size_t size_of_input  = sizeof(uint8_t) * snps * indiv;
    size_t size_of_output = sizeof(float) * indiv * indiv;
    // Create handle
    cublasCreate(&handle);

    // Allocate memory
    err = cudaMalloc(&d_A, size_of_input);  
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);;
    err = cudaMalloc(&d_B, size_of_input);  
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);;
    err = cudaMalloc((void**)&d_C, size_of_output);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    err = cudaMallocHost((void**)&h_C, size_of_output);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);

    // Copy data to device
    err = cudaMemcpy(d_A, snp_matrix, size_of_input, cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    err = cudaMemcpy(d_B, snp_matrix, size_of_input, cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    
    // Calculate GEMM
    status =
        cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_T, nrowA, nrowA, ncolA,
                     &alpha, d_A, input_type, nrowA, d_B, input_type, nrowA,
                     &beta, d_C, output_type, nrowA, compute_type, algo);
    cudaDeviceSynchronize();
    if (checkError(__func__, __LINE__, status) != 0)
        return (1);

    // Copy data back to host
    err = cudaMemcpy(h_C, d_C, size_of_output, cudaMemcpyDeviceToHost);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);

    // Cast to double
    for(int i = 0; i < indiv * indiv; i++){
        ans[i] = (double) (h_C[i]);
    }

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFreeHost(h_C);
    cublasDestroy(handle);
    return 0;
}

}