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
                               int indiv, double *ans) {   
    /*
    xxx

    */
                              
    // Get number of threads
    cudaError_t err         = cudaSuccess,
                private_err = cudaSuccess;
    cudaStream_t stream;

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
    if (verbose) {
        printf("Using %d OMP threads.\n", num_threads);
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
    debug_info("Using tile size of %d.\n", mem_tile_size);

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

    // Main loop
    // Calculates matrix multiplications in parallel: Each thread in this loop sends
    // its data to a different stream on the device. The threads calculate
    // concurrently and send the output back to main memory. Memory copies are
    // asynchronous to take full advantage of the memory bandwidth.
    #ifdef DO_PARALLEL
    // #pragma omp parallel for num_threads(num_threads) private(private_err,stream) shared(err) schedule(dynamic)
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

        int rows_remaining = indiv - i;
        int x_tile_size = min(mem_tile_size, rows_remaining);

        // private_err = cudaMemcpyAsync(d_tile1, x, x_tile_size * n_bytes_per_indiv,
        //                     cudaMemcpyHostToDevice, stream);
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

        for (long j = i; j < indiv; j += mem_tile_size) {
            if (err != cudaSuccess) {
              continue;
            }

            unsigned char *y = snp_matrix + j * n_bytes_per_indiv;

            int columns_remaining = indiv - j;
            int y_tile_size = min(mem_tile_size, columns_remaining);
            
            debug_info("i: %d, j: %d, snps_per_byte: %d, bytes_per_snp: %d, columns_remaining: %d, x_tile_size: %d,  y_tile_size: %d, mem_tile_size: %d", i, j, n_snps_per_byte, n_bytes_per_indiv, columns_remaining, x_tile_size, y_tile_size, mem_tile_size);
            
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
            private_err = cudaGetLastError();
            if (checkError(__func__, __LINE__, private_err) != 0) {
                printf("Thread = %d, i = %d\n", threadidx, i);
                err = private_err;
                continue;
            }

            // initialize gemm arguments
            CutlassGemm::Arguments args(
                {int(x_tile_size), int(y_tile_size),
                int(n_bytes_per_indiv_padded * n_snps_per_u4b)},
                {d_tile1, int(n_bytes_per_indiv_padded * n_snps_per_u4b)},
                {d_tile2, int(n_bytes_per_indiv_padded * n_snps_per_u4b)},
                {d_M + threadidx * mem_tile_size * mem_tile_size, int(y_tile_size)},
                {d_M + threadidx * mem_tile_size * mem_tile_size, int(y_tile_size)},
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
            private_err = cudaMemcpyAsync(h_M + threadidx * mem_tile_size * mem_tile_size,
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
                              d2 + d1 * y_tile_size);
                  double *ans0 = ans + (i + d1), *ans1 = ans + (i + d1) * indiv;

                  ans0[(j + d2) * indiv] = (double)Mij;
                  ans1[j + d2] = (double)Mij;
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

    // Free memory
    cudaFree(d_Z_block1);
    cudaFree(d_Z_block2);
    cudaFree(d_M);
    cudaFreeHost(h_M);

    return 0;
}




void err_check(const char* string){
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s %s\n", string, cudaGetErrorString(err)); 
}
static void gpuCrossprodIntern_legacy(unsigned int *CGM, size_t snps, size_t individuals, double *ans)
{
#define CodesPerByte 4L//8L / 2L
#define BitsPerCode 2
#define MEMORY_FACTOR 4L
#define n_streams 10L
#define COMPRESSION_GPU 2L
#define PADDIM 4L
    
    using data_type = cutlass::uint4b_t;

    size_t TileSize = 2048;
    // force multiples of 32 byte
    const size_t BytesPerRowPadded =
        (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / 32) * 32;
    const size_t BytesPerRow =  (1 + (snps - 1) / 4L) * 4L / sizeof(unsigned int);

    printf("snps %d, individuals %d, bytesperrow %d, intsperrow %d, intermediate %d\n", snps, individuals, BytesPerRowPadded, BytesPerRow, BytesPerRow);

    // sanity checks
    // limit Tilesize to individuals
    TileSize = TileSize > individuals ? individuals : TileSize;

    // get total GPU memory
    size_t free_mem;
    cudaMemGetInfo(&free_mem, nullptr);

    // calculates total memory requirements
    size_t req_mem = 2 * BytesPerRowPadded * TileSize * CodesPerByte +
                     TileSize * TileSize * sizeof(unsigned int);
    if (req_mem > free_mem)
        printf("Not enough memory available.");

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    printf("Using device %s\n", prop.name);

    // Input data
    data_type *d_x, *d_y;
    // Buffer for output
    int32_t *d_val;
    // Buffer for copying back results from device to host
    int32_t *h_val;

    const int size_of_input =
        BytesPerRowPadded * TileSize * CodesPerByte / MEMORY_FACTOR;
    const int size_of_output = sizeof(int32_t) * TileSize * TileSize;
    // Initialization of buffers: We calculate n_streams of tile matrix
    // multiplications in parallel and allocate the corresponding amount of
    // memory
    cudaMalloc((void **)&d_x, n_streams * size_of_input);
    cudaMalloc((void **)&d_y, n_streams * size_of_input);
    cudaMalloc((void **)&d_val, n_streams * size_of_output);
    cudaMallocHost((void **)&h_val, n_streams * size_of_output);
    err_check("Memory allocation: ");

    // initialization of cutlass gemm operators
    using ColumnMajor = cutlass::layout::ColumnMajor;
    using RowMajor = cutlass::layout::RowMajor;
    using TensorOp = cutlass::arch::OpClassTensorOp;

    /* Links:
    https://github.com/NVIDIA/cutlass/blob/master/media/docs/functionality.md
    https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s8t_s8n_s32n_tensor_op_s32_sm75.cu
        */
    using ElementA_ = data_type;
    using LayoutA_ = RowMajor;
    using ElementB_ = data_type;
    using LayoutB_ = ColumnMajor;
    using ElementC_ = int32_t;
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
        ElementA_, // Data-type of A matrix
        LayoutA_,  // Layout of A matrix
        ElementB_, // Data-type of B matrix
        LayoutB_,  // Layout of B matrix
        ElementC_, // Data-type of C matrix
        LayoutC_,  // Layout of C matrix
        ElementAccumulator_, OperatorClass_, ArchTag_, ThreadblockShape_,
        WarpShape_, InstructionShape_, EpilogueOutputOp_, ThreadblockSwizzle_,
        Stages, AlignmentA, AlignmentB, SplitKSerial, cutlass::arch::CustomOp,
        IsBetaZero>;

    // Define a CUTLASS GEMM type
    CutlassGemm gemm_operator;

// Main loop
// Calculates matrix multiplications in parallel: Each thread in this loop sends its data to a different stream on the device. The threads calculate concurrently and send the output back to main memory. Memory copies are asynchronous to take full advantage of the memory bandwidth.
#ifdef DO_PARALLEL
// #pragma omp parallel for num_threads(n_streams < 1 + (individuals - 1) / TileSize ? n_streams : 1 + (individuals - 1) / TileSize) schedule(dynamic)
#endif
  for (int64_t i = 0; i < individuals; i += TileSize)
  {
    int threadidx = omp_get_thread_num();
    cudaStream_t stream;
    cudaError_t err = cudaStreamCreate(&stream);
    cudaStreamSynchronize(stream);

    if (err != cudaSuccess)
      printf("Stream couldn't be created");

    // Pointer to the first element of current rows
    unsigned int *x = (CGM + i * BytesPerRow);
    cutlass::uint4b_t *x_dev =
        (cutlass::uint4b_t *)(d_x + threadidx * TileSize * BytesPerRowPadded);
    cutlass::uint4b_t *y_dev =
        (cutlass::uint4b_t *)(d_y + threadidx * TileSize * BytesPerRowPadded);

    // Number of rows in matrix
    size_t const rows_left = individuals - i;
    // Size x of current tile
    size_t const x_tile_size = TileSize < rows_left ? TileSize : rows_left;

    cudaMemcpy2DAsync(x_dev, BytesPerRowPadded, x, 1 + (snps - 1) / CodesPerByte,
                      1 + (snps - 1) / CodesPerByte, x_tile_size,
                      cudaMemcpyHostToDevice, stream);
    cudaStreamSynchronize(stream);
    err_check("Copy 1:");

    // Inner loop
    for (int64_t j = i; j < individuals; j += TileSize)
    {

      // Same as above with y
      size_t const columns_left = individuals - j;
      size_t const y_tile_size = TileSize < columns_left ? TileSize : columns_left;
      unsigned int *y = (CGM + j * BytesPerRow);

    // cudaMemcpyAsync(y_dev, y, y_tile_size * BytesPerRowPadded, cudaMemcpyHostToDevice, stream);
      cudaMemcpy2DAsync(y_dev, BytesPerRowPadded, y, 1 + (snps - 1) / CodesPerByte,
                        1 + (snps - 1) / CodesPerByte, x_tile_size,
                        cudaMemcpyHostToDevice, stream);

      err_check("Copy 2:");
      cudaStreamSynchronize(stream);

      // initialize gemm arguments
      CutlassGemm::Arguments args(
          {int(x_tile_size), int(y_tile_size), int(BytesPerRowPadded * CodesPerByte / COMPRESSION_GPU)},
          {x_dev, int(BytesPerRowPadded * CodesPerByte / COMPRESSION_GPU)},
          {y_dev, int(BytesPerRowPadded * CodesPerByte / COMPRESSION_GPU)},
          {d_val + threadidx * TileSize * TileSize, int(y_tile_size)},
          {d_val + threadidx * TileSize * TileSize, int(y_tile_size)},
          {1, 0});
      cudaStreamSynchronize(stream);

      cutlass::Status status;
      // This needs to be critical, because cutlass isn't thread safe!!! Undefined behaviour and wrong results otherwise
#pragma omp critical
      status = gemm_operator(args, nullptr, stream);
      cudaStreamSynchronize(stream);
      if (status != cutlass::Status::kSuccess)
        printf("GEMM operation failed\n");
      err_check("GEMM operation:");

      // Copy results back to host
      cudaMemcpyAsync(h_val + threadidx * TileSize * TileSize, d_val + threadidx * TileSize * TileSize, TileSize * TileSize * sizeof(int32_t), cudaMemcpyDeviceToHost, stream);
      err_check("Copying back:");
      cudaStreamSynchronize(stream);

// Loop over tile and store values in output matrix
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(n_streams) schedule(static)
#endif
      for (int64_t di = 0; di < x_tile_size; ++di)
      {
        for (int64_t dj = 0; dj < y_tile_size; ++dj)
        {
          // Get result
          const auto Mij = *(h_val + threadidx * TileSize * TileSize + dj + di * y_tile_size);
          // Create pointers to the output matrix (because it is symmetric we use two pointers)
          double *ans0 = ans + (i + di),
                 *ans1 = ans + (i + di) * individuals;
          // Store result in ouput matrix
          ans0[(j + dj) * individuals] = ans1[j + dj] = (double)Mij;
        }
      }
    }
    cudaStreamDestroy(stream);
  }

  // Free memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_val);
  cudaFreeHost(h_val);
}


extern "C" {

void crossprod_mmagpu(unsigned char *snp_matrix, int snps, int indiv,
                      double *ans) {
    gpuCrossprodIntern(snp_matrix, snps, indiv, ans);
}

}
