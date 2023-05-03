
/* 
 Authors
 Alexander Freudenberg, afreuden@mail.uni-mannheim.de

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

#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cuda_runtime_api.h>
#include <cublasLt.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <chrono>
#include <vector>
#include "parallel_simd.h"
#ifdef TIME_AVAILABLE
#include <time.h>
#endif

#include "Basic_utils_local.h"
#include "errors_messages.h"
#include "RandomFieldsUtils.h"
#include "solve_gpu.h"
#include "options.h"
#include "xport_import.h"




ASSERT_SIMD(solve_61, gpu);

__global__ void logdet_kernel(double *d_matrix, Uint *d_size, double *d_logdet){
    __shared__ double logdet_loc;
    __shared__ double submatrix[THREADS_PER_BLOCK];
    logdet_loc = 0.0;
    *d_logdet = 0.0;
    int idx = blockDim.x * blockIdx.x + threadIdx.x,
        thread = threadIdx.x;
    if(idx < *d_size){
    //     if(THREADS_PER_BLOCK<=thread && PL >= PL_RECURSIVE)
	//   PRINTF("Size %d, access %d",THREADS_PER_BLOCK,thread );
        submatrix[thread] = d_matrix[idx * (*d_size +1)];
    }
    __syncthreads();
    atomicAdd(&logdet_loc, idx >= *d_size ? 0 : (LOG(submatrix[thread])));

    __syncthreads();
    if (threadIdx.x == 0) {atomicAdd(d_logdet, logdet_loc);
    };
};

int cholGPU(bool copy, double *matrix, Uint input_size, double *B,
	    Uint rhs_cols,
     double *LogDet, double *RESULT){
    /*
        This function solves the problem
            A x = b
        on   an available GPU and writes the solution to the original memory
        Input: 
            matrix: pointer to rowwise allocated matrix A
            individuals: number of individuals in matrix, i.e. dimension
            vector: pointer to vector b
        Ouput:
            vector: contains solution x after the function has been called
    */

  KEY_type *KT = KEYT();
  installNrun_options *iNr = &(KT->global_utils.installNrun);
  int *devices = iNr->gpu_devices;
  int N = iNr->Ngpu_devices;
  assert(iNr->Ngpu_devices <= MAX_GPU_DEVICES);
  int maxStreams = iNr->maxStreams;

#ifdef TIME_AVAILABLE  
    clock_t start = clock();
#endif    
    //declare/define process variables
    Ulong size = (Ulong) input_size;
    int bufferSize = 0;
    int *info = NULL;
    int h_info = 0;
    double *buffer = NULL;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;

    //declare device variables
    double *d_matrix = NULL;
    double *d_B = NULL;
    double *d_logdet = NULL;
    Uint *d_size = NULL;

    //initialize handle and stream, calculate buffer size needed for cholesky
    cusolverDnCreate(&handle);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cusolverDnSetStream(handle, stream);

    cusolverDnDpotrf_bufferSize(handle, uplo, size, matrix,
        size, &bufferSize);
    cudaDeviceSynchronize();
    //PRINTF("Buffersize: %f", ((float) bufferSize)/1073741824.0);
    cudaMalloc(&info, sizeof(int));
    cudaMalloc(&buffer, sizeof(double) * bufferSize);
    //allocate memory on device  
    cudaMalloc((void**)&d_matrix, sizeof(double) * size * size);
    cudaMalloc((void**)&d_B, sizeof(double) * size * rhs_cols);
    cudaMemset(info, 0, sizeof(int));

    // if (PL > PL_RECURSIVE)
    //   PRINTF("Size of alloc %ld",  sizeof(double) * size * size);
    //copy data to device
    cudaMemcpy(d_matrix, matrix, sizeof(double) * size * size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, sizeof(double) * size * rhs_cols, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    //write cholesky factorization to device copy of A
    cusolverDnDpotrf(handle, uplo, size,
            d_matrix, size, buffer, bufferSize, info);
            
    //Synchronize is necessary, otherwise error code "info" returns nonsense 
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("%s\n", cudaGetErrorString(err));

    //check for errors
    cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    if (0 != h_info) {
        if(h_info >0)PRINTF("Error: Cholesky factorization failed at minor %d \n", h_info);
        if(h_info <0)PRINTF("Error: Wrong parameter in cholesky factorization at %d entry\n", h_info);
        err = cudaDeviceReset();
        if(err != cudaSuccess)PRINTF("Device reset not successful");
        return(1);
    }
    //calculate x = A\b
    cusolverDnDpotrs(handle, uplo, size, rhs_cols, 
            d_matrix, size, d_B,
             size, info);

    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Potrs: %s\n", cudaGetErrorString(err));
    
    if(LogDet != NULL){
        cudaMalloc((void**)&d_logdet, sizeof(double));
        cudaMalloc((void**)&d_size, sizeof(Uint));
        cudaMemcpy(d_size, &size, sizeof(Uint), cudaMemcpyHostToDevice);
        logdet_kernel <<< (size - 1)/THREADS_PER_BLOCK +1 ,THREADS_PER_BLOCK>>> (d_matrix, d_size, d_logdet);
        cudaDeviceSynchronize();
        cudaMemcpy(LogDet, d_logdet, sizeof(double), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaFree(d_size);
        cudaFree(d_logdet);
    }
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Err at Logdet: %s\n", cudaGetErrorString(err));

     //*LogDet = 1.0;
    //copy  solution from device to vector on host
    cudaMemcpy(RESULT, d_B, sizeof(double) * size * rhs_cols, cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if (err != cudaSuccess) PRINTF("Memcpy: %s\n", cudaGetErrorString(err));
    //free allocated memory
    cudaFree(info);
    cudaFree(buffer);
    cudaFree(d_matrix);
    cudaFree(d_B);
    cusolverDnDestroy(handle);
    cudaStreamDestroy(stream);
#ifdef TIME_AVAILABLE
    // if (PL >= PL_RECURSIVE)
    // PRINTF("Time: %.3f", (double)(clock() - start) / CLOCKS_PER_SEC);
#endif    
    return 0;
};





// void mgpuSolve(double *matrix, Uint individuals, double *vector){
//     /*
//         This function solves the problem
//             A x = b
//         on an MULTIPLE GPUs and writes the solution to the original memory of b
//         Input: 
//             matrix: pointer to rowwise allocated matrix A
//             individuals: number of individuals in matrix, i.e. dimension
//             vector: pointer to vector b
//         Ouput:
//             vector: contains solution x after the function has been called
//     */

//     // Define auxiliary variables
//     cusolverMgHandle_t handle = NULL;
//     const int max_devices = 8; // Maximum number of devices to be used
//     int nbGpus = 0;
//     std::vector<int> deviceList;
//     const int N = individuals, lda = N; // Dimension of matrix
//     const int IA  = 1;
//     const int JA  = 1;
//     const int T_A = 256; //Tile size
//     const int IB  = 1;
//     const int JB  = 1;
//     const int T_B = 1000, ldb = N; 
//     int info = 0;
//     int64_t lwork_potrf = 0, lwork_potrs = 0, lwork = 0 ;

//     cudaLibMgMatrixDesc_t descrA, descrB;
//     cudaLibMgGrid_t grid;
//     double **array_d_A = NULL;
//     double **array_d_B = NULL;
//     double **array_d_work = NULL;

//     // Create handles and select devices
//     cusolverStatus_t status = cusolverMgCreate(&handle);
//     if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Handle couldn't be created");
    
//     cudaError_t cudaStat = cudaGetDeviceCount( &nbGpus );
//     nbGpus = (nbGpus < max_devices)? nbGpus : max_devices;
//     cudaDeviceProp prop;
//     cudaGetDeviceProperties(&prop, 0);
//     int cc_major = prop.major, cc_minor = prop.minor;
//     for(int i = 0; i< nbGpus; i++){
//         cudaDeviceProp prop;
//         cudaGetDeviceProperties(&prop, i);
//         if(prop.major == cc_major & prop.minor == cc_minor)
//                 deviceList.push_back(i);
//     }
//     nbGpus = deviceList.size();
//     status = cusolverMgDeviceSelect(
//         handle,
//         nbGpus,
//         &deviceList[0]);
//     if(CUSOLVER_STATUS_SUCCESS != status) PRINTF("Devices couldn't be selected.");

//     // Enable peer access for selected devices
//     for(int i = 0; i < nbGpus; i++){
//         cudaSetDevice(deviceList[i]);
//         for(int j = 0; j< nbGpus; j++){
//             if(i == j)continue;
//             cudaStat = cudaDeviceEnablePeerAccess(deviceList[j],0);
//             if(cudaStat != cudaSuccess)PRINTF("Device %d can't access device %d.",deviceList[i],deviceList[j]);
//             PRINTF("Access enabled for devices (%d,%d)",deviceList[i],deviceList[j]);
//         }
//     }
//     // Create device grid for vectors A, B
//     status = cusolverMgCreateDeviceGrid(&grid, 1, nbGpus, &deviceList[0], CUDALIBMG_GRID_MAPPING_COL_MAJOR );
//     if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Grid couldn't be created.");

//     // Creeate matrix descriptions
//     status = cusolverMgCreateMatrixDesc(
//         &descrA,
//         N,   /* nubmer of rows of (global) A */
//         N,   /* number of columns of (global) A */
//         N,   /* number or rows in a tile */
//         T_A, /* number of columns in a tile */
//         CUDA_R_64F,
//         grid );
//     if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Matrix descriptions couldn't be created.");
//     status = cusolverMgCreateMatrixDesc(
//         &descrB,
//         N,    /* nubmer of rows of (global) B */
//         1, /* number of columns of (global) B */
//         N,    /* number or rows in a tile */
//         T_B,  /* number of columns in a tile */
//         CUDA_R_64F,
//         grid );
//     if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Matrix description B couldn't be created.");


//     // Allocate arrays of device pointers which point at the memory allocated on each device
//     array_d_A = (double**) MALLOC (sizeof(double*) * nbGpus );
//     array_d_B = (double**)MALLOC(sizeof(double*)*nbGpus);
//     array_d_work = (double**)MALLOC(sizeof(double*)*nbGpus);
//     MEMSET(array_d_work, 0, sizeof(void*)*nbGpus);

//     // Calculate block size on device
//     const int A_num_blks = ( N + T_A - 1) / T_A;
//     const int B_num_blks = ( N + T_B - 1) / T_B;
//     const int A_blks_per_device = (A_num_blks + nbGpus-1)/nbGpus;
//     const int B_blks_per_device = (B_num_blks + nbGpus-1)/nbGpus;

//     // Allocate memory on each device
//     for( int p = 0 ; p < nbGpus ; p++){
//         cudaSetDevice(deviceList[p]);
//         cudaStat = cudaMalloc( &(array_d_A[p]), sizeof(double)*lda*T_A*A_blks_per_device );
//         if(cudaSuccess != cudaStat)PRINTF("Memory for matrix A couldn't be allocated on device %d.",deviceList[p]);
//         cudaStat = cudaMalloc( &(array_d_B[p]), sizeof(double)*ldb*T_B*B_blks_per_device );
//         if(cudaSuccess != cudaStat)PRINTF("Memory for matrix B couldn't be allocated on device %d.",deviceList[p]);
//     }

//     // Copy arrays A and B to device
//     for( int k = 0 ; k < A_num_blks ; k++){
//     /* k = ibx * nbGpus + p */
//         const int p   = (k % nbGpus);
//         const int ibx = (k / nbGpus);
//         double *h_Ak = matrix + (size_t)lda*T_A*k;
//         double *d_Ak = array_d_A[p] + (size_t)lda*T_A*ibx;
//         const int width = MIN( T_A, (N - T_A*k) );
//         cudaStat = cudaMemcpy(d_Ak, h_Ak, sizeof(double)*lda*width, cudaMemcpyHostToDevice);
//         if(cudaSuccess != cudaStat)PRINTF("Matrix A couldn't be copied at block (%d, %d).", p,ibx);
//     }
//     for( int k = 0 ; k < B_num_blks ; k++){
//     /* k = ibx * nbGpus + p */
//         const int p   = (k % nbGpus);
//         const int ibx = (k / nbGpus);
//         double *h_Bk = vector + (size_t) T_B*k;
//         double *d_Bk = array_d_B[p] + (size_t) T_B*ibx;
//         cudaStat = cudaMemcpy(d_Bk, h_Bk, sizeof(double)*T_B, cudaMemcpyHostToDevice);
//         if(cudaSuccess != cudaStat)PRINTF("Matrix B couldn't be copied at block (%d, %d).", p,ibx);
//     }

//     // Calculate buffersizes necessary for potrf and potrs
//     cudaDeviceSynchronize();
//     status = cusolverMgPotrf_bufferSize(
//         handle,
// 		CUBLAS_FILL_MODE_LOWER,
//         N,
//         (void**)array_d_A,
//         IA, /* base-1 */
//         JA, /* base-1 */
//         descrA,
//         CUDA_R_64F,
//         &lwork_potrf);
//     if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Buffer size potrf couldn't  be calculated");    
//     cudaDeviceSynchronize();
//     status = cusolverMgPotrs_bufferSize(
//         handle,
// 		CUBLAS_FILL_MODE_LOWER,
//         N,
//         1, /* number of columns of B */
//         (void**)array_d_A,
//         IA,
//         JA,
//         descrA,
//         (void**)array_d_B,
//         IB,
//         JB,
//         descrB,
//         CUDA_R_64F,
//         &lwork_potrs);
//     cudaDeviceSynchronize();
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess) PRINTF("Buffersize calculation: %s\n", cudaGetErrorString(err));
//     if(CUSOLVER_STATUS_SUCCESS != status)PRINTF("Buffer size potrs couldn't  be calculated");    

//     lwork = (lwork_potrf > lwork_potrs)? lwork_potrf : lwork_potrs;

//     // Allocate workspace size
//     for(int idx = 0 ; idx < nbGpus ; idx++){
//         int deviceId = deviceList[idx];
//         cudaSetDevice( deviceId );
//         void *d_workspace = NULL;
//         cudaStat = cudaMalloc(&d_workspace, sizeof(double)*lwork);
//         if( cudaSuccess != cudaStat )PRINTF("Workspace couldn't be allocated.");
//         ((void**)array_d_work )[idx] = d_workspace;
//     }

//     // Calculate potrf to workspace
//     status = cusolverMgPotrf(
//         handle,
// 		CUBLAS_FILL_MODE_LOWER,
//         N,   
//         (void**)array_d_A,
//         IA,
//         JA,
//         descrA,
//         CUDA_R_64F,
//         (void**)array_d_work,
//         lwork,
//         &info  /* host */
//     );
//     cudaDeviceSynchronize();
//     if(CUSOLVER_STATUS_SUCCESS != status) PRINTF("Potrf couldn't be calculated");
//     if(info != 0)PRINTF("Info code %d", info);
//     // Calculate potrs to B
//     status = cusolverMgPotrs(
//         handle,
// 		CUBLAS_FILL_MODE_LOWER,
//         N,
//         1, /* number of columns of B */
//         (void**)array_d_A,
//         IA,
//         JA,
//         descrA,
//         (void**)array_d_B,
//         IB,
//         JB,
//         descrB,
//         CUDA_R_64F,
//         (void**)array_d_work,
//         lwork,
//         &info  /* host */
//     );
//     cudaDeviceSynchronize();
//     if(CUSOLVER_STATUS_SUCCESS != status) PRINTF("Potrs couldn't be calculated");
//     if(info != 0)PRINTF("Info code %d", info);

//     // Copy solution B back to host
//     for( int k = 0 ; k < B_num_blks ; k++){
//     /* k = ibx * nbGpus + p */
//         const int p   = (k % nbGpus);
//         const int ibx = (k / nbGpus);
//         double *h_Bk = vector + (size_t) T_B*k;
//         double *d_Bk = array_d_B[p] + (size_t) T_B*ibx;
//         cudaStat = cudaMemcpy(h_Bk, d_Bk, sizeof(double)*T_B, cudaMemcpyDeviceToHost);
//         if(cudaSuccess != cudaStat)PRINTF("Matrix B couldn't be copied at block (%d, %d).", p,ibx);
//     }

//     // Free memory on device and host
//     for(int i = 0; i< nbGpus; i++){
//         cudaSetDevice(deviceList[i]);
//         cudaDeviceReset();
//     }
//     FREE(array_d_A); FREE(array_d_B); FREE(array_d_work);
// }
