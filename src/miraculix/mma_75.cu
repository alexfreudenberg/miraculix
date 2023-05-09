!
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

#define BitsPerCode 2
#define MY_VARIANT VARIANT_GPU
#define MY_LDABITALIGN MY_LDABITALIGN_2BIT

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

// CUTLASS includes
#include "cutlass/gemm/device/gemm.h"
#include "cutlass/numeric_types.h"
#include "cutlass/arch/mma.h"
#include "cutlass/layout/matrix.h"
#include "cutlass/numeric_types.h"
#include "cutlass/arch/wmma.h"
#include "cutlass/array.h"
#include "cutlass/numeric_types.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/arch/arch.h"
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm.h"

#include "mmagpu_cutlass.h"
#include "mmagpu.h"

#ifdef DO_PARALLEL
#include "omp.h"

#define n_streams 10L
#define COMPRESSION_GPU 2L
#define MEMORY_FACTOR 4L
#define PADDIM 4L

using data_type = cutlass::uint4b_t;

void err_check(const char *string)
{
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
    printf("%s %s\n", string, cudaGetErrorString(err));// Rprint
}

__global__ static void print_kernel(int32_t *ptr)
{
#ifdef SCHLATHERS_MACHINE
  printf("Value on GPU is %d \n", *ptr); // Rprint
#endif
}

// Test for right compute capability
void check_cc(int device, int compute_capability)
{
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device);
  int device_compute_capability = deviceProp.major * 10 + deviceProp.minor;

  if (device_compute_capability <  compute_capability)
    ERR1("Targeted function is not suitable for devices of compute capability", device_compute_capability);
}

void check_7_5(){
  check_cc(0,75);
}

static void gpuCrossprodIntern(unsigned int *CGM, size_t snps,
			       size_t individuals,
                               double *ans, size_t TileSize)
{

  KEY_type *KT = KEYT();
  utilsoption_type *global_utils = &(KT->global_utils);
  int *devices = global_utils->installNrun.gpu_devices;
  int N = global_utils->installNrun.Ngpu_devices;
  assert(N <= MAX_GPU_DEVICES);
  int maxStreams = global_utils->installNrun.maxStreams;
  int CORES = global_utils->basic.cores;
  
  // force multiples of 32 byte
  const size_t BytesPerRow = (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / 32) * 32;
  const size_t IntsPerRow = 1 + (BytesPerRow - 1) / sizeof(unsigned int);

  // sanity checks
  // limit Tilesize to individuals
  TileSize = TileSize > individuals ? individuals : TileSize;

  // Determine calculation device among devices
  cudaDeviceProp prop;
  MEMSET(&prop, 0, sizeof(cudaDeviceProp));
  int device = -1;

  // calculates total memory requirements
  const int size_of_input = BytesPerRow * TileSize * CodesPerByte / MEMORY_FACTOR;
  const int size_of_output = sizeof(int32_t) * TileSize * TileSize;
  size_t req_mem = n_streams * ( 2 * size_of_input + size_of_output);
  size_t free_mem;

  for(int m = 0; m < N; m++){
    cudaGetDeviceProperties(&prop, devices[m]);
    if(prop.major * 10 + prop.minor < 75) continue;
    
    cudaSetDevice(devices[m]);
    cudaMemGetInfo(&free_mem, nullptr);
    if (req_mem < free_mem){
      device = devices[m];
      break;
    }
  }
  if(device < 0)
    ERR0("No suitable device with enough memory available");

  PRINTF("Using device %s\n", prop.name);

  // Input data
  data_type *d_x;
  data_type *d_y;
  // Buffer for output
  int32_t *d_val;
  // Buffer for copying back results from device to host
  int32_t *h_val;

  // Initialization of buffers: We calculate n_streams of tile matrix multiplications in parallel and allocate the corresponding amount of memory
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
  using ThreadblockShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
      OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
      ElementAccumulator_>::ThreadblockShape;
  using WarpShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
      OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
      ElementAccumulator_>::WarpShape;
  using InstructionShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
      OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
      ElementAccumulator_>::InstructionShape;
  using EpilogueOutputOp_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
      OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
      ElementAccumulator_>::EpilogueOutputOp;
  using ThreadblockSwizzle_ = typename cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>;
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
  using Operator_ = typename cutlass::gemm::device::DefaultGemmConfiguration<OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_, ElementAccumulator_>::Operator;
  const bool IsBetaZero = false;

  /* Links:
  https://github.com/NVIDIA/cutlass/blob/master/media/docs/functionality.md
  https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s8t_s8n_s32n_tensor_op_s32_sm75.cu
  */

  using CutlassGemm = cutlass::gemm::device::Gemm<
      ElementA_, // Data-type of A matrix
      LayoutA_,  // Layout of A matrix
      ElementB_, // Data-type of B matrix
      LayoutB_,  // Layout of B matrix
      ElementC_, // Data-type of C matrix
      LayoutC_,  // Layout of C matrix
      ElementAccumulator_,
      OperatorClass_,
      ArchTag_,
      ThreadblockShape_,
      WarpShape_,
      InstructionShape_,
      EpilogueOutputOp_,
      ThreadblockSwizzle_,
      Stages,
      AlignmentA,
      AlignmentB,
      SplitKSerial,
      cutlass::arch::CustomOp,
      IsBetaZero>;

  // Define a CUTLASS GEMM type
  CutlassGemm gemm_operator;

// Main loop
// Calculates matrix multiplications in parallel: Each thread in this loop sends its data to a different stream on the device. The threads calculate concurrently and send the output back to main memory. Memory copies are asynchronous to take full advantage of the memory bandwidth.
    int cores = n_streams < 1 + (individuals - 1) / TileSize ? n_streams
      : 1 + (individuals - 1) / TileSize
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(cores)) schedule(dynamic)
#endif
  for (int64_t i = 0; i < individuals; i += TileSize)
  {
    int threadidx = omp_get_thread_num();
    cudaStream_t stream;
    cudaError_t err = cudaStreamCreate(&stream);
    cudaStreamSynchronize(stream);

    if (err != cudaSuccess)
      ERR0("Stream couldn't be created");

    // Pointer to the first element of current rows
    unsigned int *x = (CGM + i * IntsPerRow);
    data_type *x_dev = d_x + threadidx * TileSize * BytesPerRow;
    data_type *y_dev = d_y + threadidx * TileSize * BytesPerRow;

    // Number of rows in matrix
    size_t const rows_left = individuals - i;
    // Size x of current tile
    size_t const x_tile_size = TileSize < rows_left ? TileSize : rows_left;

    cudaMemcpyAsync(x_dev, x, x_tile_size * BytesPerRow, cudaMemcpyHostToDevice, stream);

    cudaStreamSynchronize(stream);
    err_check("Copy 1:");

    // Inner loop
    for (int64_t j = i; j < individuals; j += TileSize)
    {

      // Same as above with y
      size_t const columns_left = individuals - j;
      size_t const y_tile_size = TileSize < columns_left ? TileSize : columns_left;
      unsigned int *y = (CGM + j * IntsPerRow);

      cudaMemcpyAsync(y_dev, y, y_tile_size * BytesPerRow, cudaMemcpyHostToDevice, stream);
      err_check("Copy 2:");
      cudaStreamSynchronize(stream);

      // initialize gemm arguments
      CutlassGemm::Arguments args(
          {int(x_tile_size), int(y_tile_size), int(BytesPerRow * CodesPerByte / COMPRESSION_GPU)},
          {x_dev, int(BytesPerRow * CodesPerByte / COMPRESSION_GPU)},
          {y_dev, int(BytesPerRow * CodesPerByte / COMPRESSION_GPU)},
          {d_val + threadidx * TileSize * TileSize, int(y_tile_size)},
          {d_val + threadidx * TileSize * TileSize, int(y_tile_size)},
          {1, 0});
      cudaStreamSynchronize(stream);

      // compute Multiplication
      cutlass::Status status = gemm_operator(args, nullptr, stream);
      cudaStreamSynchronize(stream);
      err_check("Calculation:");

      // Copy results back to host
      cudaMemcpyAsync(h_val + threadidx * TileSize * TileSize, d_val + threadidx * TileSize * TileSize, TileSize * TileSize * sizeof(int32_t), cudaMemcpyDeviceToHost, stream);
      err_check("Copying back:");

      cudaStreamSynchronize(stream);

      if (*(h_val + threadidx * TileSize * TileSize) == 0)
      {
#if defined SCHLATHERS_MACHINE
        PRINTF("Computation failed at thread %d, (%d,%d)\n", threadidx, i, j);
        print_kernel<<<1, 1>>>((int32_t *)d_val + threadidx * TileSize * TileSize);
#endif

        j -= TileSize;
        continue;
      }
      err_check("Copy back:");

// Loop over tile and store values in output matrix
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(GreaterZero(n_streams)) schedule(static) // not cores)
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
        // PRINTF("\n");
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

static void crossprodIntern(Uint *CM, Uint snps, Uint individuals, Uint lda,
                            double *ans)
{
  // tilse_size needs to stay the same: for smaller values we experience undocumented calculation failures on the device
  const size_t tilesize = 2048;

  // Initialize host pointers and copy input data cuda managed memory
  Uint *h_CM;
  const size_t BytesPerIndiv = lda * SizeOfInt;
  cudaMallocHost((void **)&h_CM, individuals * BytesPerIndiv);
  MEMCOPY(h_CM, CM, individuals * BytesPerIndiv);

  gpuCrossprodIntern(h_CM, snps, individuals, ans, tilesize);
  cudaFreeHost(h_CM);
}


void crossprod_mmagpu(Uint *CGM, Uint snps, Uint individuals, Uint lda,
                      double *ans)
{

  crossprodIntern(CGM, snps, individuals, lda, ans);
}

#endif
