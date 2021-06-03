/* 
    This file contains helper functions to speed-up SNP matrix simulations. 
    It can be dynamically loaded in an R script with dyn.load
    Compilation:
// nvcc -Xcompiler -fpic -I /usr/share/R/include --shared -lcurand -lcublas sample*cu -o sample_snpmatrix.so && Rscript RFutils_ex.R 
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


extern "C" {

#define MIN(A, B)( A<B ? A : B)
#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)

int samplesnps(int* dim, int* seed, int* M)
{
    curandGenerator_t gen;
    unsigned int *devData, *hostData;
    int n = *dim;
    int tilesize = 2048 * 2048;
    for(int i = 0; i< n; i+=tilesize){
        /* Allocate n floats on host */
        cudaMallocHost((void**)&hostData, tilesize * sizeof(unsigned int));
        /* Allocate n floats on device */
        CUDA_CALL(cudaMalloc((void **)&devData, sizeof(unsigned int) * tilesize));

        /* Create pseudo-random number generator */
        CURAND_CALL(curandCreateGenerator(&gen, 
                    CURAND_RNG_PSEUDO_DEFAULT));
        
        /* Set seed */
        CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, 
            *seed));
        CURAND_CALL(curandSetGeneratorOrdering(gen, 
            CURAND_ORDERING_PSEUDO_SEEDED));
        /* Generate n floats on device */
        CURAND_CALL(curandGenerate(gen, devData, tilesize));
        cudaDeviceSynchronize();
        /* Copy device memory to host */
        CUDA_CALL(cudaMemcpy(hostData, devData, tilesize * sizeof(unsigned int),
            cudaMemcpyDeviceToHost));
        for(int j = 0; j < MIN(tilesize, n - i); j++) M[i + j] = hostData[j]%3;
    };
    /* Cleanup */
    CURAND_CALL(curandDestroyGenerator(gen));
    CUDA_CALL(cudaFree(devData));
    cudaFree(hostData);    
    return EXIT_SUCCESS;
}

void check_gpu_available(){
    int count = 0;
    cudaGetDeviceCount(&count);
    if(count==0) ERR("No CUDA device found");
}

void err_check(const char* string){
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s %s\n", string, cudaGetErrorString(err)); 
}


static void finalizer(SEXP ptr){
    cudaError_t err =  cudaFree(R_ExternalPtrAddr(ptr));
    printf("Free status %d\n", err);
    err_check("Freeing error: ");
    printf("Finalized\n");
}

__global__ static void print_kernel(double* ptr){
    printf("Value on GPU is %f \n", *ptr);
}

void* copy_to_device(double* matrix, int size){
    /* Copies data to the device */
    //printf("Pointer is %ld\n", (long) d_matrix);
    void* d_matrix =NULL;
    cudaMalloc((void**)&d_matrix, sizeof(double) *size);  
    cudaMemcpy(d_matrix, matrix, sizeof(double) * size, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    err_check("Copying failed");

    return(d_matrix);
}