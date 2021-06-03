
/*
 Authors 
 Alexander Freudenberg, afreuden@mail.uni-mannheim.de


 Copyright (C) 2021  Alexander Freudenberg

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, writne to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/



#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cuda.h>
#include <curand.h>
#include <cublasLt.h>
#include <cublas_v2.h>
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include "rgpu.h"
#include <omp.h>
#include <chrono>

extern "C" {

#define MIN(A, B)( A<B ? A : B)
#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)


void err_check(const char* string){
    // Short inline error checking function
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s %s\n", string, cudaGetErrorString(err)); 
}


static void finalizer(SEXP ptr){
    // Finalizer function to clean up device pointers once R's gabarge collector finds them
    cudaError_t err =  cudaFree(R_ExternalPtrAddr(ptr));
    printf("Free status %d\n", err);
    err_check("Freeing error: ");
    printf("Finalized\n");
}

__global__ static void print_kernel(double* ptr){
    printf("Value on GPU is %f \n", *ptr);
}

void* copy_to_device(double* matrix, int size){
    /* Copies data to the device and returns the device pointer */
    void* d_matrix =NULL;
    cudaMalloc((void**)&d_matrix, sizeof(double) *size);  
    cudaMemcpy(d_matrix, matrix, sizeof(double) * size, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    err_check("Copying failed");

    return(d_matrix);
}

SEXP cu_matmat(SEXP mat_A, SEXP mat_B, int64_t tilesize){
    /* 
        Function to calculate the matrix product with cublas. A,B need to be checked for conforming dimensions first
    */
    
    // Initialize
    clock_t start_total = clock();
    int64_t nrow_A = nrows(mat_A),
        ncol_A = ncols(mat_A),
        nrow_B = nrows(mat_B),
        ncol_B = ncols(mat_B);
    // if(ncol_A != ncol_B){
    //      printf("Non confirming matrix dimensions (%d %d), (%d %d)",nrow_A, ncol_A, nrow_B, ncol_B);
    //      return(R_NilValue);
    // }
    tilesize = min(min(tilesize, nrow_A), nrow_B);
    double* A=REAL(mat_A);
    double* B=REAL(mat_B);

    //Allocate output matrix
    SEXP mat_C=PROTECT(allocMatrix(REALSXP, ncol_A, ncol_B));
    double* C=REAL(mat_C);
    double *d_A=NULL, *d_B=NULL, *d_C=NULL; 
    double *h_A = NULL, *h_B = NULL, *h_C=NULL;

    // Get free memory
    size_t free_mem;
    cudaMemGetInfo(&free_mem, NULL);
    err_check("Memory could not be queried");

    // Set number of streams
    int n_streams = min( 16, (int) (free_mem / (  (nrow_A +  nrow_B + tilesize) * tilesize * sizeof(double)) ) );
    n_streams = min(n_streams, (int) ((nrow_A-1)/tilesize +1) );

    // Allocate memory
    cudaMallocHost((void**)&h_C, n_streams * tilesize * tilesize * sizeof(double)); 
    err_check("Memory host could not be allocated");

    cudaMalloc((void**)& d_A, n_streams * tilesize * nrow_A *sizeof(double));
    err_check("Memory A could not be allocated");

    cudaMalloc((void**)& d_B, n_streams * tilesize * nrow_B * sizeof(double));
    err_check("Memory B could not be allocated");

    cudaMalloc((void**)& d_C, n_streams * tilesize * tilesize *sizeof(double));
    err_check("Memory C could not be allocated");
    // printf("Alloc time total %.3f\n",(double)(clock() - start_total) / CLOCKS_PER_SEC );

    // Iterate over tiles of matrix multiplications
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(n_streams) schedule(dynamic)   
    #endif
    for(int64_t i =0; i < ncol_A; i+=tilesize ){
        // Create stream
        int threadidx = omp_get_thread_num();
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(& stream);
        if(err)printf("Thread %d", threadidx);
        err_check("Stream could not be created");
        cublasHandle_t handle;
        cublasCreate(&handle);
        err_check("Handle could not be created");
        int64_t tile_n = min(tilesize, ncol_A -i);

        // Copy data to the device
        cudaMemcpyAsync((void*) ( d_A + threadidx * tilesize * nrow_A), 
        A + i * nrow_A, tile_n * nrow_A * sizeof(double), cudaMemcpyHostToDevice, stream);
        err_check("Memory A couldn't be copied:");

        double alpha = 1;
        double beta = 0;
        for(int64_t j = 0; j< ncol_B; j+=tilesize){
            //Copy data
            int64_t tile_m = min(tilesize, ncol_B -j);   
            cudaMemcpyAsync((void*) (d_B + threadidx * tilesize * nrow_B), 
            B + j * nrow_B, tile_m * nrow_B* sizeof(double), cudaMemcpyHostToDevice, stream);
            err_check("Memory B couldn't be copied:");

            cudaStreamSynchronize(stream);

            //Calculate the crossproduct and check for errros
            cublasStatus_t stat = cublasDgemm(handle,
                CUBLAS_OP_T,
                CUBLAS_OP_N,
                tile_n,
                tile_m,
                nrow_A,
                &alpha,
                d_A + threadidx * tilesize * nrow_A,
                nrow_A, 
                d_B + threadidx * tilesize * nrow_B,
                nrow_B,
                &beta,
                d_C + threadidx * tilesize * tilesize,
                tile_n
                );
            cudaStreamSynchronize(stream);
            clock_t start = clock();

            if(stat!= CUBLAS_STATUS_SUCCESS) printf("Gemmex failed with error code %d\n", stat);
            err_check("Matmul error:");
            
            // Copy output back
            cudaMemcpyAsync((void*) (h_C + threadidx * tilesize * tilesize), d_C + threadidx * tilesize * tilesize, sizeof(double) * tile_m * tile_n, cudaMemcpyDeviceToHost, stream);
            cudaStreamSynchronize(stream);            
            err_check("Memory C couldn't be copied:");
            
            // Write output to matrix
            double *sub_C = C + j * ncol_A + i;
            #ifdef DO_PARALLEL
            #pragma omp parallel for num_threads(n_streams)
            #endif
            for(int64_t Ci =0; Ci <  tile_m; Ci++ ){
            //    printf(" %d",Ci);
                memcpy(sub_C + Ci * ncol_A, h_C + threadidx * tilesize * tilesize + Ci * tile_n, sizeof(double) * tile_n  );
            }
        }
        // Free memory
        cublasDestroy(handle);
        err_check("Handle could not be destroyed");
        cudaStreamDestroy(stream);
        err_check("Stream could not be destroyed");

    }
    cudaFree(d_A);cudaFree(d_B);cudaFree(d_C);
    cudaFreeHost(h_C);
    err_check("Pointers could not be freed");
    
    return(mat_C);
}

void cu_matadd(void* d_A, void* d_B, int ncol_A, int nrow_A){
    /*
        Calculates the sum of two matrices on the GPU. Useful for large matrices.
    */

    // Initialize
    cublasHandle_t handle;
    cublasCreate(&handle);
    double alpha = 1;
    double beta = 1;

     //Calculate the sum and check for errros
     cublasStatus_t stat = cublasDgeam(handle,
        CUBLAS_OP_N,
        CUBLAS_OP_N,
        nrow_A,
        ncol_A,
        &alpha,
        (double*) d_A,
        ncol_A, 
        &beta,
        (double*) d_B,
        ncol_A,
        (double*) d_A,
        nrow_A
        );
        cudaDeviceSynchronize();
        if(stat!= CUBLAS_STATUS_SUCCESS) printf("Gemmex failed with error code %d\n", stat);
        print_kernel <<<1,1>>> ((double*) d_A + ncol_A * nrow_A -1);
        
        err_check("Matmul error:");
    cublasDestroy(handle);

}


SEXP cu_matmul(SEXP mat_A, SEXP mat_B, SEXP operation){
/*
    Matrix multiplication on a CUDA device
    This function is supposed accepts matrices A,B and calculates its product. The operation input allows for different operations, currently only matrix multiplication is implemented.
    
*/
    enum gemm_ops {MATMAT, MATVEC, MATADD};
    int tilesize = 2048;
    int protects = 0;
    SEXP ptr, ncol_return, nrow_return, ans;
    int op = INTEGER(operation)[0];
   
    switch(op){
        case MATMAT:
            ans = cu_matmat(mat_A, mat_B, tilesize);
            protects +=1; 
            break;
        // Possible extensions for large matrix vector multiplications: Not yet implemented because fast enough on CPU
        // case MATVEC: cu_matvec(d_A, d_B); break;
        // case MATADD: 
        //     cu_matadd(d_A, d_B, ncol_A, nrow_A); 
        //     d_C = d_A;
        //     cudaFree(d_B);
        //     break;
        default: printf("Unknown operation"); return(R_NilValue);
    }

    UNPROTECT(protects);
    return(ans);
}

SEXP access_pointer(SEXP input_list){
    /*
        Function to copy data from a GPU device pointer back to main memory and return it as an R list
    */

    //Initialize 
    SEXP ans, ptr;
    int protects = 0;
    int nrow = INTEGER(VECTOR_ELT(input_list,1))[0];
    int ncol = INTEGER(VECTOR_ELT(input_list,2))[0];
    ptr = VECTOR_ELT(input_list,0);

    // Allocate R data matrix
    PROTECT(ptr);
    protects +=1;
    ans = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    protects +=1;
    double *h_ptr = REAL(ans);

    // Access address of external pointer
    void* cuda_ptr = R_ExternalPtrAddr(ptr);
    printf("Pointer is %ld\n", (long) cuda_ptr);

    // Copy content
    cudaMemcpy(h_ptr, cuda_ptr, sizeof(double) * ncol * nrow, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    
    //Error checking
    err_check("Copying from device to host failed\n ");
    printf("Copying finished");
    UNPROTECT(protects);
    printf("Unprotec");
    return(ans);
}


SEXP CUDA_pointer(){
    SEXP ans, ptr;
    void *d_matrix = NULL;
    cudaMalloc((void**)&d_matrix, sizeof(double) * 1);
    double value = 1.0;

    cudaMemcpy(d_matrix, &value, sizeof(double) * 1, cudaMemcpyHostToDevice);
    print_kernel <<<1,1>>> ((double*)d_matrix);
    cudaDeviceSynchronize();


    ans = PROTECT(allocVector(INTSXP,1));
    ptr = R_MakeExternalPtr(d_matrix, R_NilValue,R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizerEx(ptr, finalizer, TRUE);
    printf("Pointer is %ld\n", (long) d_matrix);
    printf("Pointer is %ld\n", (long) R_ExternalPtrAddr(ptr));


    INTEGER(ans)[0]= 1;

    UNPROTECT(2);
    return(ptr);
}




}
