/*
 Authors 
 Alexander Freudenberg, alexander.freudenberg@uni-mannheim.de
 
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



// This file supplies the internal PTX-ISA algorithm for 2-bit matrix multiplication 
// An extra file is needed to make the package compilable on devices with compute capability 7.0 or lower 

#ifndef miraculix_mmagpuIntern_H
#define miraculix_mmagpuIntern_H 1

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
