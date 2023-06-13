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

// Include order is important!!!
// Namespace conflicts if order is changed
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#include "cuda_utils.h"
#include "mmagpuIntern.h"

const int default_tile_size = 2048;

int gpuCrossprodIntern(unsigned char *snp_matrix, int snps,
                               int indiv, double *ans, bool is_plink_format) {   
    /*
    xxx

    */
                              
    // Get number of threads
    cudaError_t err         = cudaSuccess,
                private_err = cudaSuccess;
    cudaStream_t stream;

    clock_t start, setup, calculation, end;
    start = clock();
    // Input data
    unsigned char *d_Z_block1, *d_Z_block2;
    // Buffer for output
    int32_t *d_M, *h_M;

    bool verbose = get_print_level() >= 0;
    const char *env_num_threads = getenv("OMP_NUM_THREADS");
    int num_threads = 4;
    if (env_num_threads != NULL) {
        num_threads = atoi(env_num_threads);
    }


    const long n_bytes_per_indiv =
        (snps - 1) / 4 + 1; // number of columns of Z if individuals
                             // are zero-padded to be a multiple of 4
    const long n_bytes_per_indiv_padded =
        (1 + (n_bytes_per_indiv - 1) / 32) * 32;
    // number of columns of Z if individuals
    // are zero-padded to be a multiple of 32 bytes
    const long n_snps_per_byte = 8 / 2;
    const int n_snps_per_u4b = 4 / 2;

    // sanity checks
    // limit Tilesize to individuals
    const char *env_tile_size = getenv("TILE_SIZE");
    int mem_tile_size = default_tile_size;
    if (env_tile_size != NULL) {
        mem_tile_size = atoi(env_tile_size);
    }
    
    if (verbose) {
        printf("Using %d OMP threads and memory tile size of %d.\n", num_threads, mem_tile_size);
    }

    debug_info("Function arguments: snps %d, indiv %d, bytes_per_indiv %d, bytes_per_indiv_padded %d, k %d.\n", snps, indiv, n_bytes_per_indiv, n_bytes_per_indiv_padded, n_bytes_per_indiv_padded * n_snps_per_u4b);
    
    mem_tile_size = min(indiv, mem_tile_size);

    if (checkCuda() != 0) {
        return 1;
    }
    int device = switchDevice();
    if (device == -1) {
        return 1;
    }

    // Calculate total memory requirements
    size_t required_mem = num_threads * (2 * n_bytes_per_indiv_padded * mem_tile_size +
                          mem_tile_size * mem_tile_size * sizeof(unsigned int));
    if (checkDevMemory(required_mem) != 0) {
        return 1;
    }

    int size_of_input = n_bytes_per_indiv_padded * mem_tile_size;
    int size_of_output = sizeof(int) * mem_tile_size * mem_tile_size;
    // Initialization of buffers: Calculate num_threads of tile matrix
    // multiplications in parallel and allocate the corresponding amount of
    // memory
    err = cudaMalloc((void **)&d_Z_block1, num_threads * size_of_input);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    err = cudaMalloc((void **)&d_Z_block2, num_threads * size_of_input);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    err = cudaMalloc((void **)&d_M, num_threads * size_of_output);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    err = cudaMallocHost((void **)&h_M, num_threads * size_of_output);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);

    // initialization of cutlass gemm operators
    using ColumnMajor = cutlass::layout::ColumnMajor;
    using RowMajor = cutlass::layout::RowMajor;
    using TensorOp = cutlass::arch::OpClassTensorOp;
    using ElementA_ = cutlass::uint4b_t;
    using LayoutA_ = RowMajor;
    using ElementB_ = cutlass::uint4b_t;
    using LayoutB_ = ColumnMajor;
    using ElementC_ = int;
    using LayoutC_ = RowMajor;
    using ElementAccumulator_ = ElementC_;
    using OperatorClass_ = TensorOp;
    using ArchTag_ = cutlass::arch::Sm75;
    using ThreadblockShape_ =
        typename cutlass::gemm::device::DefaultGemmConfiguration<
            OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
            ElementAccumulator_>::ThreadblockShape;
    using WarpShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::WarpShape;
    using InstructionShape_ =
        typename cutlass::gemm::device::DefaultGemmConfiguration<
            OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
            ElementAccumulator_>::InstructionShape;
    using EpilogueOutputOp_ =
        typename cutlass::gemm::device::DefaultGemmConfiguration<
            OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
            ElementAccumulator_>::EpilogueOutputOp;
    using ThreadblockSwizzle_ =
        typename cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>;
    const int Stages = cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::kStages;
    const int AlignmentA = cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::kAlignmentA;
    const int AlignmentB = cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::kAlignmentB;
    const bool SplitKSerial = false;
    using Operator_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
        OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
        ElementAccumulator_>::Operator;
    const bool IsBetaZero = false;
    using CutlassGemm = cutlass::gemm::device::Gemm<
        ElementA_, LayoutA_, ElementB_, LayoutB_, ElementC_, LayoutC_,
        ElementAccumulator_, OperatorClass_, ArchTag_, ThreadblockShape_,
        WarpShape_, InstructionShape_, EpilogueOutputOp_, ThreadblockSwizzle_,
        Stages, AlignmentA, AlignmentB, SplitKSerial, cutlass::arch::CustomOp,
        IsBetaZero>;

    // Define a CUTLASS GEMM type
    CutlassGemm gemm_operator;

    if (checkError(__func__, __LINE__, cudaGetLastError()) != 0)
        return (1);

    err = cudaMemcpyToSymbol(d_conversion_table, PLINK_CONVERSION_TABLE, 256 * sizeof(unsigned char), 0, cudaMemcpyHostToDevice);
    if (checkError(__func__, __LINE__, err) != 0)
        return (1);
    replace_functor functor(d_conversion_table);


    setup = clock();
    debug_info("Crossproduct: Time for initializing: %.3fs",
                (double)(setup - start) / CLOCKS_PER_SEC);

    // Main loop
    // Calculates matrix multiplications in parallel: Each thread in this loop sends
    // its data to a different stream on the device. The threads calculate
    // concurrently and send the output back to main memory. Memory copies are
    // asynchronous to take full advantage of the memory bandwidth.
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(num_threads) private(private_err,stream) shared(err) schedule(dynamic)
    #endif
    for (long i = 0; i < indiv; i += mem_tile_size) {
        if(err != cudaSuccess){
            continue;
        }

        int threadidx = omp_get_thread_num();

        private_err = cudaStreamCreate(&stream);
        if (checkError(__func__, __LINE__, private_err) != 0) {
            printf("Thread = %d, i = %d\n", threadidx, i);
            err = private_err;
            continue;
        }

        cudaStreamSynchronize(stream);

        cutlass::uint4b_t *d_tile1 =
            (cutlass::uint4b_t *)(d_Z_block1 +
                                  threadidx * mem_tile_size * n_bytes_per_indiv_padded);
        cutlass::uint4b_t *d_tile2 =
            (cutlass::uint4b_t *)(d_Z_block2 +
                                  threadidx * mem_tile_size * n_bytes_per_indiv_padded);

        unsigned char *x = snp_matrix + i * n_bytes_per_indiv;

        int rows_remaining     = indiv - i;
        int x_tile_size        = min(mem_tile_size, rows_remaining);
        int x_tile_size_padded = ((x_tile_size-1)/4 +1 ) * 4; // padded tile sizes are required as CUTLASS needs dimensions that are a multiple of 4

        private_err = cudaMemcpy2DAsync(
            d_tile1, n_bytes_per_indiv_padded, x, n_bytes_per_indiv,
            n_bytes_per_indiv, x_tile_size, cudaMemcpyHostToDevice,
            stream);
        

        cudaStreamSynchronize(stream);
        if (checkError(__func__, __LINE__, private_err) != 0) {
            printf("Thread = %d, i = %d\n", threadidx, i);
            err = private_err;
            continue;
        }
        if(!is_plink_format){
            device_convert_plink_2bit((unsigned char*) d_tile1, n_bytes_per_indiv_padded * x_tile_size, functor);
        }
        private_err = cudaGetLastError();
        if (checkError(__func__, __LINE__, private_err) != 0) {
            printf("Thread = %d, i = %d\n", threadidx, i);
            err = private_err;
            continue;
        }
        
        for (long j = i; j < indiv; j += mem_tile_size) {
            if (err != cudaSuccess) {
              continue;
            }

            unsigned char *y = snp_matrix + j * n_bytes_per_indiv;

            int columns_remaining = indiv - j;
            int y_tile_size = min(mem_tile_size, columns_remaining);
            int y_tile_size_padded = ((y_tile_size-1)/4 +1 ) * 4; // padded tile sizes are required as CUTLASS needs dimensions that are a multiple of 4
                        
            private_err = cudaMemcpy2DAsync(
                d_tile2, n_bytes_per_indiv_padded, y, n_bytes_per_indiv,
                n_bytes_per_indiv, y_tile_size, cudaMemcpyHostToDevice,
                stream);

            cudaStreamSynchronize(stream);
            if (checkError(__func__, __LINE__, private_err) != 0) {
                printf("Thread = %d, i = %d, j = %d\n", threadidx, i, j);
                err = private_err;
                continue;
            }

            if(!is_plink_format){
                device_convert_plink_2bit((unsigned char*) d_tile2, n_bytes_per_indiv_padded * y_tile_size, functor);
            }
            private_err = cudaGetLastError();
            if (checkError(__func__, __LINE__, private_err) != 0) {
                printf("Thread = %d, i = %d\n", threadidx, i);
                err = private_err;
                continue;
            }

            // initialize gemm arguments
            CutlassGemm::Arguments args(
                {int(x_tile_size_padded), int(y_tile_size_padded),
                int(n_bytes_per_indiv_padded * n_snps_per_u4b)},
                {d_tile1, int(n_bytes_per_indiv_padded * n_snps_per_u4b)},
                {d_tile2, int(n_bytes_per_indiv_padded * n_snps_per_u4b)},
                {d_M + threadidx * mem_tile_size * mem_tile_size, int(y_tile_size_padded)},
                {d_M + threadidx * mem_tile_size * mem_tile_size, int(y_tile_size_padded)},
                {1, 0});
            cudaStreamSynchronize(stream);

            // compute Multiplication
            cutlass::Status status;
        #pragma omp critical
            status = gemm_operator(args, nullptr, stream);

            cudaStreamSynchronize(stream);
            if (checkError(__func__, __LINE__, (cudaError_t) status) != 0 ) {
                printf("Thread = %d, i = %d, j = %d, status = %d\n", threadidx, i, j, status);
                err = (cudaError_t) status;
                continue;
            }

            private_err = cudaGetLastError();
            if (checkError(__func__, __LINE__, private_err) != 0) {
                printf("Thread = %d, i = %d, j = %d\n", threadidx, i, j);
                err = private_err;
                continue;
            }

            // Copy results back to host
            private_err =
                cudaMemcpyAsync(h_M + threadidx * mem_tile_size * mem_tile_size,
                                d_M + threadidx * mem_tile_size * mem_tile_size,
                                mem_tile_size * mem_tile_size * sizeof(int),
                                cudaMemcpyDeviceToHost, stream);

            cudaStreamSynchronize(stream);
            if (checkError(__func__, __LINE__, private_err) != 0) {
                printf("Thread = %d, i = %d, j = %d\n", threadidx, i, j);
                err = private_err;
                continue;
            }

            for (long d1 = 0; d1 < x_tile_size; d1++) {
                for (long d2 = 0; d2 < y_tile_size; d2++) {
                    // Get result
                    int Mij = *(h_M + threadidx * mem_tile_size * mem_tile_size +
                                d2 + d1 * y_tile_size_padded);

                    // Write results to M[i,j] and M[j,i]  
                    ans[(i + d1) + (j + d2) * indiv] = (double)Mij;
                    ans[(j + d2) + (i + d1) * indiv] = (double)Mij;
                }
            }
        }

        cudaStreamSynchronize(stream);
        private_err = cudaGetLastError();
        if (checkError(__func__, __LINE__, private_err) != 0) {
            printf("Thread = %d, i = %d\n", threadidx, i);
            err = private_err;
            continue;
        }

        cudaStreamDestroy(stream);
    }

    calculation = clock();
    debug_info("Crossproduct: Time for calculation: %.3fs",
                (double)(calculation - setup) / CLOCKS_PER_SEC);

    // Free memory
    cudaFree(d_Z_block1);
    cudaFree(d_Z_block2);
    cudaFree(d_M);
    cudaFreeHost(h_M);

    end = clock();
    debug_info("Crossproduct: Total time: %.3fs",
                (double)(end - start) / CLOCKS_PER_SEC);
    return 0;
}

void device_convert_plink_2bit(unsigned char *d_block, unsigned int extent, replace_functor functor){
  thrust::device_ptr<unsigned char> thrust_ptr(d_block);
//   thrust::transform(thrust_ptr, thrust_ptr + extent, [&](unsigned char value){ return PLINK_CONVERSION_TABLE[value];});
  thrust::transform(thrust_ptr, thrust_ptr + extent, thrust_ptr, functor);
}

extern "C" {

int crossprod_mmagpu(unsigned char *snp_matrix, int snps, int indiv,
                      double *ans, bool is_plink_format) {
    return gpuCrossprodIntern(snp_matrix, snps, indiv, ans, is_plink_format);
}

}
