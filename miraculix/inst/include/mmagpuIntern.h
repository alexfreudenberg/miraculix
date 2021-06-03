
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2019 -- 2019  Martin Schlather

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


// Alex Diese Datei macht nur Sinn, wenn wir mehrere Versionen von mmapgu.cu
// haben (werden), die sich nur durch Vorab-#define's unterscheiden. Ansonsten gehoert der Inhalt in die Datei mmagpu.cu
// Vgl. shuffleIntern.h und shuffle128.cc, shuffle256.cc

#ifndef miraculix_mmagpuIntern_H
#define miraculix_mmagpuIntern_H 1

#include <inttypes.h>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

// templates are to controll amount of bytes used for computation
// ! IMPORTANT T must be unsigned !

#ifndef PACKEDGPU_CUH
#define PACKEDGPU_CUH


// CUTLASS includes
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/system/cuda/memory.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/discard_iterator.h>

#include <chrono>
#include <cstring>

// Defines cutlass::gemm::device::Gemm, the generic Gemm computation template class.
#include "cutlass/gemm/device/gemm.h"
#include "cutlass/numeric_types.h"



// CUTLASS includes
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
#include "cutlass/util/host_tensor.h"
#include "cutlass/util/reference/device/gemm.h"
#include "cutlass/util/reference/host/tensor_compare.h"
#include "cutlass/util/reference/host/tensor_copy.h"
#include "cutlass/util/reference/host/tensor_fill.h"
#include "cutlass/util/tensor_view_io.h"

#define CUTLASS_ARCH_MMA_SM75_ENABLED 1

namespace cutlass{
namespace arch{
    struct CustomOp;
    template <>
struct Mma<
  gemm::GemmShape<8,8,32>,
  32,
  uint4b_t,
  layout::RowMajor,
  uint4b_t,
  layout::ColumnMajor,
  int,
  layout::RowMajor,
  CustomOp> {

  using Shape = gemm::GemmShape<8,8,32>;

  using ElementA = uint4b_t;
  using LayoutA = layout::RowMajor;
  using FragmentA = Array<uint4b_t, 8>;

  using ElementB = uint4b_t;
  using LayoutB = layout::ColumnMajor;
  using FragmentB = Array<uint4b_t, 8>;

  using ElementC = int;
  using LayoutC = layout::RowMajor;
  using FragmentC = Array<int, 2>;

  using Operator = CustomOp;
  using ArchTag = arch::Sm75;

  /// Computes multiply-add
  CUTLASS_HOST_DEVICE
  void operator()(
    FragmentC &d,
    FragmentA const &a,
    FragmentB const &b,
    FragmentC const &c
  ) const {

#if defined(CUTLASS_ARCH_MMA_SM75_ENABLED)

  unsigned const & A = reinterpret_cast<unsigned const &>(a);
  unsigned const & B = reinterpret_cast<unsigned const &>(b);

  int const *C = reinterpret_cast<int const *>(&c);
  int *D = reinterpret_cast<int *>(&d);

  asm volatile(
    "{\n\t"
    ".reg .u32 u1, u2;\n\t"
    "and.b32 u1, %2, %6;\n\t"
    "and.b32 u2, %3, %6;\n\t"
    "bar.warp.sync 0xffffffff;\n\t"
    "mma.sync.aligned.m8n8k32.row.col.s32.u4.u4.s32 {%0,%1}, {u1}, {u2}, {%4,%5};\n\t"
    "}\n"
      : "=r"(D[0]), "=r"(D[1])
      : "r"(A), "r"(B), "r"(C[0]), "r"(C[1]), "r"(0x33333333));
  
  asm volatile(
    "{\n\t"
    ".reg .u32 u1, u2;\n\t"
    ".reg .s32 s1, s2;\n\t"
    "shr.b32 u1, %2, 2;\n\t"
    "shr.b32 u2, %3, 2;\n\t"
    "and.b32 u1, u1, %6;\n\t"
    "and.b32 u2, u2, %6;\n\t"
    "mma.sync.aligned.m8n8k32.row.col.s32.u4.u4.s32 {%0,%1}, {u1}, {u2}, {%4,%5};\n\t"
    "}\n"
      : "=r"(D[0]), "=r"(D[1])
      : "r"(A), "r"(B), "r"(D[0]), "r"(D[1]), "r"(0x33333333));
#else
    assert(0);
#endif
  }
};
}
}

#include <cuda_runtime_api.h>
#include <cuda_runtime.h>

//#include  "options.h"
//#include "xport_import.h"
//#include "align.h"
//#include "haplogeno.h"
//#include "Haplo.h"
//#include <inttypes.h>
using data_type = cutlass::uint4b_t;

// Define data type on device
typedef std::conditional< std::is_same<data_type,uint8_t>::value,
                            unsigned int,
                            unsigned short>::type dev_data_type;
// For uint4b only half the memory is needed, for 2bit only a fourth is neede

// const int memory_constant = std::is_same<data_type,uint8_t>::value ? 1 : (use2bit ? 4 : 2); // Alex: diese Zeile ist mir unklar. Isnbesondere sehe ich keine Wirkung von memory_constant=1 ausser in der Speicherberechnung. Irgendwie muss doch dies auch in die Rechnungen eingehen? S. unten wo es verwendet wird
#define n_streams 10L


// ------------------------------------------------------------
// Operator that trasforms the 2-bit inout data
// into 8-bit output data and also copies the data
// from host to device
template<typename T> struct unpack_data{};

template<>
struct unpack_data<uint8_t>
{
    __host__ __device__ unsigned int operator()(const unsigned char& i) const
    {
        unsigned int result = 0;
        // unpack 2-bit into 8-bit values
        result |= ((i >> 6) & 0x03) << 0;
        result |= ((i >> 4) & 0x03) << 8;
        result |= ((i >> 2) & 0x03) << 16;
        result |= ((i >> 0) & 0x03) << 24;

        return result;
    };
};

template<>
struct unpack_data<cutlass::uint4b_t>
{
    __host__ __device__ unsigned short operator()(const unsigned char& i) const
    {
     unsigned short result = 0;

       // unpack 2-bit into 4-bit values of a 16bit integer
       result |= ((i >> 6) & 0x03) << 0;
       result |= ((i >> 4) & 0x03) << 4;
       result |= ((i >> 2) & 0x03) << 8;
       result |= ((i >> 0) & 0x03) << 12;

       return result;
    };
};
__global__ static void print_kernel(unsigned int* y, int length) {
    for(int8_t i = 0; i < length;i++) PRINTF("%d ", y[i] );
};


// ------------------------------------------------------------
// Actual crossproduct function
static void gpuCrossprodIntern(unsigned int* CGM, size_t snps, size_t individuals,
    double* ans, size_t TileSize) {

  // Martin an Alex: Keine globalen  Variablen, insb keine Optionen.
  // Gegebenenfalls muss der aktuelle Wert auf der Datenstruktur
  // gespeichert werden.
  
  
  option_type *global = &(KEYT()->global);
  int  memory_constant = global->genetics.mmamethod; // Werte 4, 2,1 gemaess oben
bool use2bit = global->genetics.mmamethod == 4;
int compressed = (use2bit ? 2 : 1); // Alex: besserer Namen?; darf/ muss bei memory_constant ==1 auch use2bit gewaehlt werden koennen? Dann brauchen wir 4 verschiedene Werte fuer global->genetics.mmamethod oder gar 2 Variablen.
    
    
    // timers
    std::chrono::time_point<std::chrono::high_resolution_clock> timer_start0, timer_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> timer_stop;
    
    timer_start0 = std::chrono::high_resolution_clock::now();
    // helpful values
    //  const size_t ValuesPerByte = 8L / 2L; Alex: sollte CodesPerBytes sein

#ifdef GPU_STANDALONE
    // standalone doesn't need alignment to specific width
    const size_t BytesPerRow = (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / sizeof(unsigned int)) * sizeof(unsigned int);
#else
    // force multiples of 32 byte
    const size_t BytesPerRow = (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / 32) * 32;
#endif
    const size_t IntsPerRow = 1 + (BytesPerRow - 1) / sizeof(unsigned int);

    // sanity checks
    // limit Tilesize to individuals
    TileSize = TileSize > individuals ? individuals : TileSize;

    // get total GPU memory
    size_t free_mem;
    cudaMemGetInfo(&free_mem, nullptr);

    const double memoryFactor = 0.05;
    const size_t totalGPUMemory = free_mem * (1 - memoryFactor);

    // calculates total memory requirements
    size_t neededGPUMemory = 2 * BytesPerRow * TileSize * CodesPerByte + TileSize * TileSize * sizeof(unsigned int);

    // adjust TileSize according to available Memory
    while (neededGPUMemory > totalGPUMemory)
    {
        PRINTF("Tilesize reduced\n");
        // reduce tilesize by 10%
        TileSize *= 0.9;

        // recalculate memory requirements
        neededGPUMemory = 2 * BytesPerRow * TileSize * CodesPerByte + TileSize * TileSize * sizeof(unsigned int);

        if (TileSize == 0)
            ERR("Not enough GPU memory available!");
    }
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    
#ifdef DEBUG
    PRINTF("Using device %s\n", prop.name);
    std::cout << std::endl << "Debug Outputs:" << std::endl;
    PRINTF("Free memory: %.2f GB\n", (float) free_mem/1073741824.0);
    std::cout << "individuals: " << individuals << std::endl;
    std::cout << "snips: " << snps << std::endl;
    std::cout << "TileSize: " << TileSize << std::endl;
    std::cout << "Integers per row: " << IntsPerRow << " = " << BytesPerRow << "B" << std::endl;
    std::cout << "Memory allocated on GPU: " << neededGPUMemory / 1073741824.0f << "GB" << std::endl;
    std::cout << "End of Debug Output" << std::endl << std::endl;
#endif // DEBUG

    // Input buffers
    data_type* raw_x;
    data_type* raw_y;
    // Buffer for output
     int32_t* raw_val;
    // Buffer for copying back results from device to host
     int32_t* temp_D2H;

     // Alex: so verstehe ich Deinen Code; kann mir aber keinen Reim
     // daraus machen.
     if (std::is_same<data_type,uint8_t>::value && memory_constant != 1)
       ERR("'memory_constant must be 1 is case of [Alex]");
     
    const int size_of_input = BytesPerRow * TileSize * CodesPerByte / memory_constant;
    const int size_of_output = sizeof(int32_t) * TileSize * TileSize;
    // Initialization of buffers
    if (cudaMalloc((void**)&raw_x, n_streams * size_of_input) != cudaSuccess)
        ERR("GPU memory allocation failed for raw_x!");

    if (cudaMalloc((void**)&raw_y, n_streams * size_of_input) != cudaSuccess)
        ERR("GPU memory allocation failed for raw_y!");

    if (cudaMalloc((void**)&raw_val, n_streams * size_of_output) != cudaSuccess)
        ERR("GPU memory allocation failed for raw_val!");

    if (cudaMallocHost((void**)&temp_D2H, n_streams * size_of_output) != cudaSuccess)
        ERR("GPU memory allocation failed for temp_D2H!");

    //initialization of cutlass gemm operators
    using ColumnMajor = cutlass::layout::ColumnMajor;
    using RowMajor = cutlass::layout::RowMajor;
    using TensorOp = cutlass::arch::OpClassTensorOp;

/* Links:
https://github.com/NVIDIA/cutlass/blob/master/media/docs/functionality.md
https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s8t_s8n_s32n_tensor_op_s32_sm75.cu
*/
using ElementA_ = data_type;
using LayoutA_  = RowMajor;
using ElementB_ = data_type;
using LayoutB_  = ColumnMajor;
using ElementC_ = int32_t;
using LayoutC_  = RowMajor;
using ElementAccumulator_ = ElementC_;
using OperatorClass_ = TensorOp;
using ArchTag_  = cutlass::arch::Sm75;
using ThreadblockShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
ElementAccumulator_>::ThreadblockShape;
using WarpShape_ = typename cutlass::gemm::device::DefaultGemmConfiguration<
OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,
ElementAccumulator_>::WarpShape ;
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
using Operator_ = typename cutlass::gemm::device::DefaultGemmConfiguration< OperatorClass_, ArchTag_, ElementA_, ElementB_, ElementC_,    ElementAccumulator_>::Operator;
const bool IsBetaZero = false;

/* Links:
https://github.com/NVIDIA/cutlass/blob/master/media/docs/functionality.md
https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s8t_s8n_s32n_tensor_op_s32_sm75.cu
*/


using CutlassGemm = cutlass::gemm::device::Gemm<
    ElementA_,        // Data-type of A matrix
    LayoutA_,         // Layout of A matrix
    ElementB_,        // Data-type of B matrix
    LayoutB_,      // Layout of B matrix
    ElementC_,              // Data-type of C matrix
    LayoutC_,        // Layout of C matrix
    ElementAccumulator_,
    OperatorClass_,            // tag indicating Tensor Cores,          // tag indicating Tensor Cores
    ArchTag_,
    ThreadblockShape_,
    WarpShape_,
    InstructionShape_,
    EpilogueOutputOp_ ,
    ThreadblockSwizzle_,
    Stages,
    AlignmentA,
    AlignmentB,
    SplitKSerial,
    cutlass::arch::CustomOp,
    IsBetaZero  >;


    // Define a CUTLASS GEMM type
    CutlassGemm gemm_operator;

    timer_stop = std::chrono::high_resolution_clock::now();
#ifdef DEBUG
    std::cout << "Initialization time: " << ((float) std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start0).count())/1000000 << "\n";
#endif
    // time measurement variables
    float copy_time = 0;
    float compute_time = 0;

    
    // Main loop
    // Iterates over upper-triangular matrix (because our matrix is symmetric)
    #ifdef DO_PARALLEL
    #pragma omp parallel for num_threads(n_streams <  1+(individuals-1)/TileSize ? n_streams : 1+(individuals-1)/TileSize) schedule(dynamic)   
    #endif
    for (int64_t i = 0; i < individuals; i += TileSize) {
        int threadidx = omp_get_thread_num();
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate ( & stream );
        if(err != cudaSuccess)
            ERR("Stream couldn't be created");
        
        // Pointer to the first element of current rows
        unsigned int *x = (CGM + i * IntsPerRow);
        data_type *x_dev = raw_x + threadidx * TileSize * BytesPerRow ;
        // Number of rows in matrix
        size_t const rows_left = individuals - i;
        // Size x of current tile
        size_t const x_tile_size = TileSize < rows_left ? TileSize : rows_left;

        timer_start = std::chrono::high_resolution_clock::now();
        //2-bit to 8-bit transformation and copy from host to device
       /* thrust::transform(thrust::cuda::par.on(stream), reinterpret_cast<unsigned char*>(x), reinterpret_cast<unsigned char*>(x) + x_tile_size * BytesPerRow,
            thrust::cuda::pointer<dev_data_type>(reinterpret_cast<dev_data_type*>(raw_x + threadidx * TileSize * BytesPerRow * CodesPerByte / memory_constant )), unpack_data<data_type>());
        */    cudaMemcpyAsync(x_dev, x, x_tile_size * BytesPerRow,cudaMemcpyHostToDevice, stream);

        timer_stop = std::chrono::high_resolution_clock::now();
        copy_time += std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();

        // Inner loop
        for (int64_t j = i; j < individuals; j += TileSize) {

            // Same as above with y
            size_t const columns_left = individuals - j;
            size_t const y_tile_size = TileSize < columns_left ? TileSize : columns_left;

            timer_start = std::chrono::high_resolution_clock::now();
            //2-bit to 8-bit transformation and copy from host to device
            unsigned int* y = (CGM + j * IntsPerRow);
            data_type *y_dev = raw_y + threadidx * TileSize * BytesPerRow ;

       /*     thrust::transform(thrust::cuda::par.on(stream), reinterpret_cast<unsigned char*>(y), reinterpret_cast<unsigned char*>(y) + y_tile_size * BytesPerRow,
                thrust::cuda::pointer<dev_data_type>(reinterpret_cast<dev_data_type*>(raw_y + threadidx *  TileSize * BytesPerRow * CodesPerByte /memory_constant)), unpack_data<data_type>());
        */    cudaMemcpyAsync(y_dev, y, y_tile_size * BytesPerRow, cudaMemcpyHostToDevice, stream);
            
            timer_stop = std::chrono::high_resolution_clock::now();
            copy_time += std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
            cudaStreamSynchronize(stream);

            //initialize gemm arguments
            CutlassGemm::Arguments args(
                { int(x_tile_size), int(y_tile_size), int(BytesPerRow * CodesPerByte  / compressed) },  // Gemm Problem dimensions
                { x_dev, int(BytesPerRow * CodesPerByte / compressed) },    // Tensor-ref for source matrix A
                { y_dev, int(BytesPerRow * CodesPerByte  / compressed) },    // Tensor-ref for source matrix B
                { raw_val + threadidx * TileSize * TileSize, int(y_tile_size) },  // Tensor-ref for source matrix C
                { raw_val + threadidx * TileSize * TileSize, int(y_tile_size) },  // Tensor-ref for destination matrix D (may be different memory than source C matrix)
                { 1, 0 });                      // scalars used in gemm (beta is 0, because we don't wand to add C to the multiplication)

            // compute Multiplication
            timer_start = std::chrono::high_resolution_clock::now();
            cutlass::Status status = gemm_operator(args, nullptr, stream);
            cudaStreamSynchronize(stream);

        	// check if calculation were successfull
            if (status != cutlass::Status::kSuccess)
                ERR("Computation failed.");

            timer_stop = std::chrono::high_resolution_clock::now();
            compute_time += std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();

            timer_start = std::chrono::high_resolution_clock::now();
            // Copy results back to host
        	if (cudaMemcpyAsync(temp_D2H + threadidx * TileSize * TileSize, raw_val + threadidx * TileSize * TileSize, TileSize * TileSize * sizeof(int32_t), cudaMemcpyDeviceToHost, stream) != cudaSuccess)
                ERR("Copying of results from device to host failed!");

            cudaStreamSynchronize(stream);
            // Loop over tile and store values in output matrix 
            #ifdef DO_PARALLEL
            #pragma omp parallel for num_threads(CORES) schedule(static)   
            #endif
            for (int64_t di = 0; di < x_tile_size; ++di)
            {
                for (int64_t dj = 0; dj < y_tile_size; ++dj)
                {
                    // Get result
                    const auto Mij = *(temp_D2H +threadidx * TileSize * TileSize + dj + di * y_tile_size ); 
                    // Create pointers to the output matrix (because it is symmetric we use two pointers)
                    double* ans0 = ans + (i + di),
                        * ans1 = ans + (i + di) * individuals;
                    // Store result in ouput matrix
                    ans0[(j + dj) * individuals] = ans1[j + dj] = (double) Mij;
                }
               // PRINTF("\n");
            }

            timer_stop = std::chrono::high_resolution_clock::now();
            copy_time += std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();

   
        }
        cudaStreamDestroy(stream);
    }

#ifdef DEBUG
    // Time output
    std::cout << "Timing: Copy time + compute time = total time" << std::endl;
    std::cout << copy_time << "us + " << compute_time << "us = " << (copy_time + compute_time) / 1000000 << "s" << std::endl;
#endif // DEBUG


    // Free memory
    cudaFree(raw_x);
    cudaFree(raw_y);
    cudaFree(raw_val);
    cudaFreeHost(temp_D2H);

}
#endif // PACKEDGPU_CUH



using uint4b_t = cutlass::uint4b_t;
/*
struct pack_data
{

    Uint operator()(const uint512_t& input) const
    {
        Uint result = 0;
        for(int i = 0; i<CodesPerUnit; i++) 
            result |= ((Uint)(input >> CodesPerUnit * i) &0x03 ) << 2*i;
        return result;
    };
};
*/

static SEXP matrix_start_Intern(Uint snps, Uint individuals, SEXP VARIABLE_IS_NOT_USED file){
  SEXP Code;
  PROTECT(Code = CreateEmptyCodeVector(snps, individuals, MY_METHOD));
  UNPROTECT(1);
  return Code;
}


static void crossprodIntern(Uint *CGM, Uint snps, Uint individuals,
          double *ans) {
  // TILE_SIZE is adjusted automatically this is just an upper bound
  const size_t TILE_SIZE = 2048 / BitsPerCode; // Alex: bitte #define verwenden fuer 2048

#ifdef DEBUG
  std::chrono::time_point<std::chrono::high_resolution_clock> timer_start, timer_stop;
  timer_start = std::chrono::high_resolution_clock::now();
#endif

  // force multiples of 32 byte
  Uint* cuda_CGM; // Alex: naechste Zeile ueberpruefen
  const size_t BytesPerIndiv = UnitsPerIndiv256(snps) * BytesPerUnit;
    // (1 + ((1 + (snps - 1) / CodesPerByte) - 1) / 32) * 32;
  cudaMallocHost((void**)&cuda_CGM, individuals * BytesPerIndiv);
  MEMCPY(cuda_CGM, CGM, individuals * BytesPerIndiv);

#ifdef DEBUG
  timer_stop = std::chrono::high_resolution_clock::now();
  float compute_time = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();  
  PRINTF("Initial copy time: %.4f\n", (float)compute_time /1000000);
#endif

  gpuCrossprodIntern(cuda_CGM, snps, individuals, ans, TILE_SIZE);
  cudaFreeHost(cuda_CGM);

#ifdef DEBUG 
  timer_stop = std::chrono::high_resolution_clock::now();
   compute_time = std::chrono::duration_cast<std::chrono::microseconds>(timer_stop - timer_start).count();
  PRINTF("Aggregated time: %.4f\n", (float) compute_time /1000000);
#endif
}



#endif 
