
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2020 -- 2021  Martin Schlather, Alexander Freudenberg

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


// This file supplies the internal PTX-ISA algorithm for 2-bit matrix multiplication 
// An extra file is needed to make the package compilable on devices with compute capability 7.0 or lower 

#ifndef miraculix_mmagpuIntern_H
#define miraculix_mmagpuIntern_H 1
#include <unistd.h>

#include <inttypes.h>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

// CUTLASS includes
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/system/cuda/memory.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/discard_iterator.h>

#include <chrono>
#include <cstring>


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
    "bar.warp.sync 0xffffffff;\n\t"
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

//////////////////////////////////////////////////
// DO NOT MOVE OR DELETE INCLUDES OR CHANGE ORDER
// very nasty compile errors caused by redefinitions
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <cublasLt.h>
#include <cublas_v2.h>

#include <omp.h>
#include "error.h"
#include "MX.h"
#include "intrinsics.h"
#include "IntrinsicsBase.h"
#include "xport_import.h"
#include "align.h"
#include "haplogeno.h"
#include "Haplo.h"
#include <inttypes.h>

// These constants might be needed in future releases for an 8-bit fallback option on older devices. 

// // Define data type on device
// typedef std::conditional< std::is_same<data_type,uint8_t>::value,
//                             unsigned int,
//                             unsigned short>::type dev_data_type;

// bool use2bit = true;
// const int memory_constant = std::is_same<data_type,uint8_t>::value ? 1 : (use2bit ? 4 : 2); 
//     int compressed = (use2bit ? 2 : 1); 

#define n_streams 10L
#define COMPRESSION_GPU 2L
#define MEMORY_FACTOR 4L
#define PADDIM 4L


using data_type = cutlass::uint4b_t;



void err_check(const char* string){
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) printf("%s %s\n", string, cudaGetErrorString(err)); 
}
__global__ static void print_kernel(int32_t* ptr){
    printf("Value on GPU is %d \n", *ptr);
}

// ------------------------------------------------------------
// Operator that unpacks 2-bit intermediary data into 8-bits 
// NOT NEEDED for 2-bit compression algorithms which unpack at calculation time

// template<typename T> struct unpack_data{};

// template<>
// struct unpack_data<uint8_t>
// {
//     __host__ __device__ unsigned int operator()(const unsigned char& i) const
//     {
//         unsigned int result = 0;
//         // unpack 2-bit into 8-bit values
//         result |= ((i >> 6) & 0x03) << 0;
//         result |= ((i >> 4) & 0x03) << 8;
//         result |= ((i >> 2) & 0x03) << 16;
//         result |= ((i >> 0) & 0x03) << 24;

//         return result;
//     };
// };

// template<>
// struct unpack_data<cutlass::uint4b_t>
// {
//     __host__ __device__ unsigned short operator()(const unsigned char& i) const
//     {
//      unsigned short result = 0;

//        // unpack 2-bit into 4-bit values of a 16bit integer
//        result |= ((i >> 6) & 0x03) << 0;
//        result |= ((i >> 4) & 0x03) << 4;
//        result |= ((i >> 2) & 0x03) << 8;
//        result |= ((i >> 0) & 0x03) << 12;

//        return result;
//     };
// };



#endif 