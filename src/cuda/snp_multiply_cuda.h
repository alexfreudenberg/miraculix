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

This file extends the CUTLASS library by a template specialization. CUTLASS is released by NVIDIA under BSD 3-Clause License, with the following license notice:

BSD 3-Clause License

Copyright (c) 2020-2023, NVIDIA CORPORATION & AFFILIATES
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


// This file supplies the internal PTX-ISA algorithm for 2-bit unsigned integer matrix multiplication 
// An extra file is needed to make the package compilable on devices with compute capability 7.0 or lower 

#pragma once

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

#include <thrust/device_vector.h>
#include <thrust/transform.h>

// Struct for format conversion in Thrust
typedef struct replace_functor replace_functor;

// Function declarations


/**
 * Computes the crossproduct of the SNP matrix with itself.
 * 
 * @param snp_matrix A pointer to the SNP matrix in compressed 2bit format.
 * @param snps The number of SNPs in the SNP matrix.
 * @param indiv The number of individuals in the SNP matrix.
 * @param ans A pointer to a preallocated array of size indiv * indiv to store the result.
 * @param is_plink_format Specifies if the SNP matrix is in PLINK binary format or already converted to 2bit encoding.
 * 
 * @return Returns an integer indicating the success (0) or failure (1) of the computation.
 * 
 * @note The SNP matrix is assumed to be ordered in SNP-major format.
 *       The computation is parallelized across streams on the GPU using a custom CUTLASS extension.
 */
int snp_crossproduct(unsigned char *snp_matrix, int snps, int indiv,
                     double *ans, bool is_plink_format);

/**
 *  \brief Convert PLINK storage format to 2bit encoding on the GPU.
 *
 *  Uses the thrust library to perform a byte-sized look-up of PLINK values in a conversion table.
 */
void device_convert_plink_2bit(unsigned char *d_block, unsigned int extent,
                               replace_functor functor);


extern "C" {

/**
 *  \brief C - Wrapper for snp_crossproduct.
 *
 *  Refer to the documentation of snp_crossproduct for details.
 */
int snp_multiply_gpu(unsigned char *snp_matrix, int snps, int indiv,
                     double *ans, bool is_plink_format);
                     
}

// This structure adds a two-bit unsigned integer (u2b) matrix multiplication
// microkernel to the cutlass namespace

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

// Set-up Plink-to-2bit conversion table
__constant__ unsigned char d_conversion_table[256]  = {0x00,0xff,0x01,0x02,0xff,0xff,0xff,0xff,0x04,0xff,0x05,0x06,0x08,0xff,0x09,0x0a,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0x10,0xff,0x11,0x12,0xff,0xff,0xff,0xff,0x14,0xff,0x15,0x16,0x18,0xff,0x19,0x1a,0x20,0xff,0x21,0x22,0xff,0xff,0xff,0xff,0x24,0xff,0x25,0x26,0x28,0xff,0x29,0x2a,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0x40,0xff,0x41,0x42,0xff,0xff,0xff,0xff,0x44,0xff,0x45,0x46,0x48,0xff,0x49,0x4a,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0x50,0xff,0x51,0x52,0xff,0xff,0xff,0xff,0x54,0xff,0x55,0x56,0x58,0xff,0x59,0x5a,0x60,0xff,0x61,0x62,0xff,0xff,0xff,0xff,0x64,0xff,0x65,0x66,0x68,0xff,0x69,0x6a,0x80,0xff,0x81,0x82,0xff,0xff,0xff,0xff,0x84,0xff,0x85,0x86,0x88,0xff,0x89,0x8a,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0x90,0xff,0x91,0x92,0xff,0xff,0xff,0xff,0x94,0xff,0x95,0x96,0x98,0xff,0x99,0x9a,0xa0,0xff,0xa1,0xa2,0xff,0xff,0xff,0xff,0xa4,0xff,0xa5,0xa6,0xa8,0xff,0xa9,0xaa};;


struct replace_functor {
    __device__
    int operator()(unsigned char value) const {
        return d_conversion_table[value];
    }
};
