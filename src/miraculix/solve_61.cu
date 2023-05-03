
/* 
 Authors
 Martin Schlather, martin.schlather@uni-mannheim.de

 (library for simulation of random fields)

 Copyright (C) 2021 -- 2021 Alexander Freudenberg

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


#define checkCudaErrors(err) {					      \
    cudaError_t e=cudaGetLastError();                                 \
    if(e!=cudaSuccess) {                                              \
      PRINTF("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
      exit(0);								\
    }									\
  }

#include <cusolverDn.h>
//#include <helper_cuda.h>
#include <cuda_runtime_api.h>
#include <cublasLt.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <chrono>

#include "Basic_miraculix.h"
#include "xport_import.h"
#include "options.h"
#include "solve_gpu.h"
//#include "MXinfo.h"
//#include "xport_import.h"

//Two small kernels for printing data on the GPU, shoul be called with 1 block and 1 thread
/*__global__ void print_kernel(int32_t* d_C, int in dividuals, int snps) {
  for (int i = 0; i < in dividuals * snps;i++) p rintf("% " PRId32, d_C[i] );
}

__global__ void print_kernel(int8_t* d_C, int in dividuals, int snps) {
    for (int i = 0; i < in dividuals * snps;p rintf( (++i % individuals) ? "" : "\n")   ) p rintf("% " PRId8, d_C[i] );
}
*/

void gpuSolve(double *matrix, Uint individuals, double *vector){
/*
    This function solves the problem
        A x = b
    on an available GPU and writes the solution to the original memory
    Input: 
        matrix: pointer to rowwise allocated matrix A
        individuals: number of individuals in matrix, i.e. dimension
        vector: pointer to vector b
    Ouput:
        vector: contains solution x after the function has been called
*/

//declare/define process variables


  
    int bufferSize = 0;
    int *info = NULL;
    int h_info = 0;
    double *buffer = NULL;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;

//declare device variables
    Uint ind_sq = individuals * individuals;
    double *d_matrix = NULL;
    double *d_b = NULL; 

//initialize handle and stream, calculate buffer size needed for cholesky
    checkCudaErrors(cusolverDnCreate(&handle));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));

    checkCudaErrors(cusolverDnSetStream(handle, stream));

    checkCudaErrors(cusolverDnDpotrf_bufferSize(handle, uplo, individuals, matrix,
        individuals, &bufferSize));
    checkCudaErrors(cudaMalloc(&info, sizeof(int)));
    checkCudaErrors(cudaMalloc(&buffer, sizeof(double) * bufferSize));
//allocate memory on device  
    checkCudaErrors(cudaMalloc((void**)&d_matrix, sizeof(double)*ind_sq));
    checkCudaErrors(cudaMalloc((void **)&d_b, sizeof(double) * individuals));
    checkCudaErrors(cudaMemset(info, 0, sizeof(int)));

//coppy data to device
    checkCudaErrors(cudaMemcpy(d_matrix, matrix, sizeof(double)*ind_sq, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_b, vector, sizeof(double) * individuals, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());

//write cholesky factorization to device copy of A
    cusolverDnDpotrf(handle, uplo, individuals,
            d_matrix, individuals, buffer, bufferSize, info);
            
    //Synchronize is necessary, otherwise error code "info" returns nonsense 
    checkCudaErrors(cudaDeviceSynchronize());

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("%s\n", cudaGetErrorString(err));

//check for errors
    checkCudaErrors(
        cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaDeviceSynchronize());

    PRINTF("Code %i\n", h_info);
    if (0 != h_info) {
        ERR0("Error: Cholesky factorization failed\n");
    }
//calculate x = A\b
    checkCudaErrors(
        cusolverDnDpotrs(handle, uplo, individuals, 1, 
            d_matrix, individuals, d_b,
             individuals, info));

    checkCudaErrors(cudaDeviceSynchronize());
    

//copy  solution from device to vector on host
    cudaMemcpy(vector, d_b, sizeof(double) * individuals, cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Memcpy: %s\n", cudaGetErrorString(err));

            
//free allocated memory
    checkCudaErrors(cudaFree(info));
    checkCudaErrors(cudaFree(buffer));
    checkCudaErrors(cudaFree(d_matrix));
    checkCudaErrors(cudaFree(d_b));
    checkCudaErrors(cusolverDnDestroy(handle));
    checkCudaErrors(cudaStreamDestroy(stream));
};

__global__ void scalar_int_naive(int8_t* d_M, double* d_A, Uint n, Uint k, Uint n_threads, size_t pitch){
    int thread = blockIdx.x * blockDim.x + threadIdx.x;
    int i = -1, j = 0;
    Uint sum = 0;

  //  printf("Pitch: %d, Necessary: %d\n", (int) pitch, k * (int) sizeof(int8_t));
  //  for (int i = 0; i< pitch * n; i++ )  printf( (i%pitch)? "%d ": "\n%d ", d_M[i]);
    if(thread < n_threads){
        //This calculates the position the thread is supposed to compute in the result matrix A
        for (; thread >= 0; j = thread + (++i), thread -= (n- i));
        //Now i,j denotes the index of A_{i,j} 
     /*   for (int m = 0; m < k; m++) sum += (Uint) (d_M[ i * k + m] * d_M[ j * k + m]); 
        d_A[ j * n + i] = d_A[ i * n + j] = (double) sum;
*/
    }
}
void gpu_relmat_custom(Uint* M, double* A, Uint snps, Uint individuals){
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
    Uint k = snps;
    Uint dim = n * k;
    cudaError_t err;
//Start timing copy and calculation time
#ifdef DEBUG
    std::chrono::time_point<std::chrono::high_resolution_clock> timer_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> timer_stop;
    timer_start = std::chrono::high_resolution_clock::now();
#endif
// allocate memory
    double *d_A;
    int8_t *h_M, *d_M;
    int32_t *h_C;
    size_t pitch_M, pitch_A;
    cudaMallocPitch((void**)&d_M, &pitch_M, (size_t)( sizeof(int8_t) * k), (size_t) n);
    cudaMallocPitch((void**)&d_A, &pitch_A,(size_t) (sizeof(double) * n),(size_t) n );
    cudaMallocHost((void **)&h_M, sizeof(int8_t) * dim);
    cudaMallocHost((void **)&h_C, sizeof(int32_t) * n * n);
    checkCudaErrors(err);

    KEY_type *KT = KEYT();
    int cores = KT->global_utils.basic.cores;

    
//Type-cast matrix M to int8 and store the result in page-locked memory
//Zero-pad matrix to get a row number that is a multiple of four
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores))   
#endif
    for (int i = 0; i < n; i++){
        for (int j = 0; j < k; j++){
        h_M[j + i * k] = (int8_t) (j< k ?  M[j + i * k] : 0 );
        }
    }


//Copy int8 matrix to device
    size_t M_width = (size_t)( sizeof(int8_t) * k);
    PRINTF("Copying: Pitch %d and h_M %d\n", (int) pitch_M * k, sizeof(int8_t) * dim);
    cudaMemcpy2D(d_M, pitch_M, h_M, M_width, M_width, n, cudaMemcpyHostToDevice);  
    checkCudaErrors(err);
//Calculate the crossproduct and check for errros
    Uint n_threads = n* (n+1)/2 ;
    scalar_int_naive <<< 1,1 >>> (d_M, d_A, n, k, n_threads, pitch_M);

//    scalar_int_naive <<< (1 + ( n_threads - 1)/ THREADS_PER_BLOCK), THREADS_PER_BLOCK >>> (d_M, d_A, n, k, n_threads);
    // PRINTF("GemmEx failed.");
    cudaDeviceSynchronize();


//copy result back to host
    cudaMemcpy(A, d_A, sizeof(double) * n * n, cudaMemcpyDeviceToHost);

//Free memory 
    cudaFree(d_M);
    cudaFree(d_A);
    cudaFreeHost(h_C);
    cudaFreeHost(h_M);
    cudaDeviceSynchronize();

//Stop timer
#ifdef DEBUG
    timer_stop = std::chrono::high_resolution_clock::now();
    PRINTF("Time: %.3f s\n", ((float) std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count())/1000000.0 );
#endif
} 



void gpu_relmat_cublas(Uint* M, double* A, Uint snps, Uint individuals){
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
    
    
	KEY_type *KT = KEYT();
  int cores = KT->global_utils.basic.cores;

    //Type-cast matrix M to int8 and store the result in page-locked memory
    //Zero-pad matrix to get a row number that is a multiple of four
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(GreaterZero(cores))   
    #endif
        for (int i = 0; i < n; i++){
            for (int j = 0; j < k_pad; j++){
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
