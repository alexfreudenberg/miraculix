
/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de


 Copyright (C) 2020 -- 2021   Alexander Freudenberg, Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


/*

fuer mich: 
tar zcvf x.tgz miraculix/miraculix RandomFieldsUtils/RandomFieldsUtils miraculix/test.R miraculix/makefile RandomFieldsUtils/test.R RandomFields/RandomFields RandomFieldsUtils/RandomFields  

 */




//////////////////////////////////////////////////
// DO NOT MOVE OR DELETE INCLUDES OR CHANGE ORDER
// very nasty compile errors caused by redefinitions

#include <string>
#include <unistd.h>
#include <chrono>

#include <cuda_runtime_api.h> 
#include <cuda_runtime.h> 

#include <cublasLt.h> 
#include <cublas_v2.h> 

#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "intrinsics_specific.h"


#include "mmagpu.h"


//////////////////////////////////////////////////
// DO NOT MOVE OR DELETE INCLUDES OR CHANGE ORDER
// very nasty compile errors caused by redefinitions

#include <cublasLt.h> 
#include <cublas_v2.h> 





#define BitsPerCode 2 
#define MY_VARIANT VARIANT_GPU
#define MY_LDABITALIGN MY_LDABITALIGN_2BIT

#define n_streams 10L
#define COMPRESSION_GPU 2L
#define MEMORY_FACTOR 4L
#define PADDIM 4L

ASSERT_SIMD(mma_61, gpu);


void crossprod_PlainGPU(Uint* M, Uint snps, Uint individuals,
			Uint VARIABLE_IS_NOT_USED lda, int cores,
			double* A ){

  KEY_type *KT = KEYT();
  utilsoption_type *global_utils = &(KT->global_utils);
  int *devices = global_utils->installNrun.gpu_devices;
  int N = global_utils->installNrun.Ngpu_devices;
  assert(N <= MAX_GPU_DEVICES);
  int maxStreams = global_utils->installNrun.maxStreams;

  
    /*
        Calculates the crossproduct of M with cublas and stores the result in A.
        Input:
            M: non-encoded matrix of dimension snps x indiv (k x n) storing genomic information
            A: pointer to result matrix
            snps: Number of snps
            individuals: number of individuals
        Output:
            A: matrix containing the type-casted result of M^T * M
        
        Note: cublas is fortran based and therefore assumes M is column-major. Therefore to calculate
            A we instruct cublasgemmex to calculate M * M^T and adjust its parameters.
            Furthermore, gemmex requires the matrix M to have a row number that is a multiple of four
            Therefore this function implements a zero-padding to add extra rows
    */
    
    //Define auxiliary variables as needed for gemmex
        Uint n = individuals;
        Uint m = individuals;
        Uint k = snps;
    
    //Auxiliary padding variables for padding
        Uint k_pad_diff = (PADDIM - k % PADDIM) % PADDIM;
        Uint k_pad = k + k_pad_diff;
        Uint dim = m * k_pad;
    
    //Start timing copy and calculation time
    #ifdef DEBUG
        std::chrono::time_point<std::chrono::high_resolution_clock> timer_start;
        std::chrono::time_point<std::chrono::high_resolution_clock> timer_stop;
        timer_start = std::chrono::high_resolution_clock::now();
    #endif
    //Declare cublas variables and allocate memory
        cublasHandle_t handle;
        cublasCreate(&handle);
        int8_t *d_M, *h_M;
        int32_t *d_C, *h_C;
        int32_t alpha = 1.f;
        int32_t beta = 0.f;
        cudaMalloc(&d_M, sizeof(int8_t) * dim);
        cudaMalloc(&d_C, sizeof(int32_t) * n * m );
        cudaMallocHost((void **)&h_M, sizeof(int8_t) * dim);
        cudaMallocHost((void **)&h_C, sizeof(int32_t) * n * m);
    
    
    
    //Type-cast matrix M to int8 and store the result in page-locked memory
    //Zero-pad matrix to get a row number that is a multiple of four
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(GreaterZero(cores))   
    #endif
        for(int i = 0; i < n; i++){
            for(int j = 0; j < k_pad; j++){
            h_M[j + i * k_pad] = (int8_t) (j< k ?  M[j + i * k] : 0 );
            }
        }
    
    
    //Copy int8 matrix to device
    cudaMemcpy(d_M, h_M, sizeof(int8_t) * dim, cudaMemcpyHostToDevice);  

    //Calculate the crossproduct and check for errros
        cublasStatus_t stat = cublasGemmEx(handle,
            CUBLAS_OP_T,
            CUBLAS_OP_N,
            n,
            m,
            k_pad,
            &alpha,
            d_M,
            CUDA_R_8I,
            k_pad, // I have no idea why this doesnt need to be individuals, same below
            d_M,
            CUDA_R_8I,
            k_pad,
            &beta,
            d_C,
            CUDA_R_32I,
            n,
            CUDA_R_32I, //CUBLAS_COMPUTE_32I,
            CUBLAS_GEMM_DEFAULT
            );
        
        if(stat) PRINTF("GemmEx failed.");
        cudaDeviceSynchronize();
    
    
    //copy result back to host
        cudaMemcpy(h_C, d_C, sizeof(int32_t) * n * m, cudaMemcpyDeviceToHost);
    
    //Convert result to double and store it in output matrix A
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(GreaterZero(cores))   
    #endif
        for (int i = 0; i < n * m; i++) A[i] = (double) h_C[i];
    
    //Free memory 
        cublasDestroy(handle);
        cudaFree(d_M);
        cudaFree(d_C);
        cudaFreeHost(h_C);
        cudaFreeHost(h_M);
    
    //Stop timer
    #ifdef DEBUG
        timer_stop = std::chrono::high_resolution_clock::now();
        PRINTF("Time: %.3f s\n", ((float) std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count())/1000000.0 );
    #endif
}


/*
This is a possible alternative implementation in cublasLt. Might be faster?
    cublasLtCreate(&handle);

    cublasLtMatmulDescCreate(&operationDesc, CUBLAS_COMPUTE_32I, CUDA_R_32I);
    cublasLtMatmulDescSetAttribute(operationDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb));
    cublasLtMatmulDescSetAttribute(operationDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb));

    cublasLtMatrixLayoutCreate(&Adesc, CUDA_R_8I, snps, in dividuals, snps);
    cublasLtMatrixLayoutCreate(&Bdesc, CUDA_R_8I, in dividuals, snps, snps);
    cublasLtMatrixLayoutCreate(&Cdesc, CUDA_R_32I, individuals, individuals, individuals);


    cublasLtMatmulPreferenceCreate(&preference);
    cublasLtMatmulPreferenceInit(preference);

//Hopefully not needed:
    cublasLtMatmulPreferenceSetAttribute(preference, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspaceSize, sizeof(workspaceSize));

    cublasLtMatmulAlgoGetHeuristic(handle, operationDesc, Adesc, Bdesc, Cdesc, Cdesc, preference, 1, &heuristicResult, &returnedResults);

    if (returnedResults == 0) {
        ERR1("Status %.50s", CUBLAS_STATUS_NOT_SUPPORTED);
    }

    cublasLtMatmul(handle,
        operationDesc,
        d_alpha,  // 1
        d_M,
        Adesc,
        d_M,
        Bdesc,
        d_beta,   // 0
        d_C,
        Cdesc,
        d_C,
        Cdesc,
        &heuristicResult.algo,
        workspace,
        workspaceSize,
        0);

    cublasLtMatmulDescDestroy(operationDesc);
    cublasLtDestroy(handle);
*/
