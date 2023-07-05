/*
 Alexander Freudenberg, alexander.freudenberg@stads.de

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

-------------------------------------------------------------------------------

This software package includes code developed by NVIDIA CORPORATION AFFILIATES, released under the
BSD 3-Clause License. The following is the license notice for that code:

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


 This file contains work which was originally published by NVIDIA CORPORATION &
AFFILIATES as part of the cutlass library and has been modified by Alexander
Freudenberg.
*/

#pragma once

#if defined(__CUDACC_RTC__)
#include <cuda/std/cassert>
#else
#include <assert.h>
#endif

#include <stdarg.h>

// CUTLASS includes
#include "cutlass/arch/mma.h"
#include "cutlass/arch/wmma.h"
#include "cutlass/array.h"
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/gemm/thread/mma.h"
#include "cutlass/gemm/threadblock/default_mma.h"
#include "cutlass/gemm/warp/mma_simt.h"
#include "cutlass/layout/matrix.h"
#include "cutlass/numeric_types.h"
#include "cutlass/subbyte_reference.h"

// This structure is used to handle a GPU object for genotype matrix
// multiplication. It contains genotype matrix data and allocated device memory
// to which double-precision matrices will be transferred.

struct GPU_gemm_storage {
  uint8_t *d_plink; // Pointer to device copy of SNP matrix in plink format
  uint8_t *d_plink_transposed; // Pointer to device copy ot transposed SNP
                               // matrix in plink format
  double *d_f;      // Pointer to device copy of vector of allele frequencies
  double *d_unit;   // Pointer to device vector of 1s
  double *d_B;      // Pointer to device vector of B
  double *d_C;      // Pointer to device vector of C
  double *d_D;      // Pointer to device vector of D
  long size_buffer; // Size of matrices B and C
  int snps;         // Number of SNPs
  int indiv;        // Number of individuals
  int device;       // Device to perform the operations on
};

// Function declarations

int plink2gpu(char *plink, char *plink_transposed, int snps,
              int indiv, double *f, int n, void **GPU_obj);
int dgemm_compressed_gpu(bool transa, void *GPU_obj, int n, double *B, int ldb,
                         int centered, int normalized, double *C, int ldc);
int freegpu(void **GPU_obj);


// Below we augment the cutlass namespace with templates for
// genotype matrix multiplication

// This part adds a packed_double data structure which contains four
// double-precision floating point values
namespace cutlass {

    template <int N>
    struct packed_double {

        /// Number of bits
        static int const kBits = N * sizeof_bits<double>::value;

        /// Storage type
        using Storage = double;

        //
        // Data members
        //

        Storage storage[N];

        //
        // Methods
        //

        /// No operation
        CUTLASS_HOST_DEVICE
        packed_double()
        {
            for (int i = 0; i < N; i++)
                storage[i] = 0.0;
        }

        CUTLASS_HOST_DEVICE
        packed_double(int value)
        {
            for (int i = 0; i < N; i++)
                storage[i] = (double)value;
        }
        CUTLASS_HOST_DEVICE
        packed_double(unsigned value)
        {
            for (int i = 0; i < N; i++)
                storage[i] = (double)value;
        }
        CUTLASS_HOST_DEVICE
        packed_double(double value)
        {
            for (int i = 0; i < N; i++)
                storage[i] = value;
        }

        /// Equality
        CUTLASS_HOST_DEVICE
        bool operator==(packed_double const& rhs) const
        {
            for (int i = 0; i < N; i++)
                if (rhs.storage[i] != storage[i])
                    return false;
            return true;
        }

        /// Inequality
        CUTLASS_HOST_DEVICE
        bool operator!=(packed_double const& rhs) const
        {
            for (int i = 0; i < N; i++)
                if (rhs.storage[i] != storage[i])
                    return true;
            return false;
        }

        /// Less than or equal
        CUTLASS_HOST_DEVICE
        bool operator<=(packed_double const& rhs) const
        {
            assert(0);
            return false;
        }

        /// Less than
        CUTLASS_HOST_DEVICE
        bool operator<(packed_double const& rhs) const
        {
            assert(0);
            return false;
        }

        /// Greater than or equal
        CUTLASS_HOST_DEVICE
        bool operator>=(packed_double const& rhs) const
        {
            assert(0);
            return false;
        }

        /// Greater than
        CUTLASS_HOST_DEVICE
        bool operator>(packed_double const& rhs) const
        {
            assert(0);
            return false;
        }
    };

    using u4f64_t = packed_double<4>;

    template <>
    struct sizeof_bits<u4f64_t> {
        static int const value = 256;
    };

} // namespace cutlass

// This structure adds a mixed-input datatype matrix multiplication microkernel
// of 2-bit integers and double-precision floats to the cutlass namespace

namespace cutlass {
namespace arch {

    template <typename LayoutA, typename LayoutB,
        typename LayoutC>
    struct Mma<
        cutlass::gemm::GemmShape<1, 1, 1>,
        1,
        uint8_t,
        LayoutA,
        cutlass::u4f64_t,
        LayoutB,
        double,
        LayoutC,
        OpMultiplyAdd> {

        using Shape = cutlass::gemm::GemmShape<1, 1, 1>;
        using Operator = OpMultiplyAdd;
        using ElementC = double;

        // This implementation of the microkernel probably isn't the most
        // efficient but tests showed that the GEMM operations is bound by data
        // transfers at the moment anyway
        CUTLASS_HOST_DEVICE
        void operator()(
            Array<double, 1>& d,
            Array<uint8_t, 1> const& a,
            Array<cutlass::u4f64_t, 1> const& b,
            Array<double, 1> const& c)
        {
            int8_t tmp;
            d[0] = c[0];
            for (int i = 0; i < 4; i++) {
                tmp = ((a[0] >> (2 * i)) & 0x3);
#ifndef NOPLINK
                tmp = max(tmp - 1, 0); // Convert PLINK format to 0,1,2
#endif
                d[0] += double(tmp) * b[0].storage[i];
            }
        }
    };

} // namespace arch

} // namespace cutlass

// The following part consists of specializations for a number of templates in
// the GEMM structure of cutlass The specializations are adapted to the new data
// format but are otherwise very similar to the ones in cutlass
namespace cutlass {
namespace gemm {

    namespace threadblock {
        template <
            /// Shape of threadblock-scoped matrix multiply operator (concept:
            /// GemmShape)
            typename Shape_,
            /// Shape of warp-level matrix multiply operator (concept: GemmShape)
            typename WarpShape_,
            /// Data type of A operand
            typename ElementA_,
            /// Data type of accumulator
            typename ElementC_,
            /// Layout of accumulator
            typename LayoutC_,
            /// Operation performed by GEMM
            typename Operator_>
        struct DefaultMmaCore<Shape_, WarpShape_, GemmShape<1, 1, 1>, ElementA_,
            layout::RowMajor, cutlass::u4f64_t, layout::ColumnMajor,
            ElementC_, LayoutC_, arch::OpClassSimt, 2, Operator_> {
            using Shape = Shape_;
            using WarpShape = WarpShape_;
            using InstructionShape = GemmShape<1, 1, 1>;
            using ElementA = ElementA_;
            using LayoutA = layout::RowMajor;
            using ElementB = cutlass::u4f64_t;
            using LayoutB = layout::ColumnMajor;
            using ElementC = ElementC_;
            using LayoutC = LayoutC_;
            using OperatorClass = arch::OpClassSimt;
            static int const PartitionsK = Shape::kK / WarpShape::kK;

            /// Default Operator
            using Operator = Operator_;

            /// Number of warps present
            using WarpCount = GemmShape<
                Shape::kM / WarpShape::kM,
                Shape::kN / WarpShape::kN,
                PartitionsK>;

            // Divisility requirements
            static_assert(
                !(Shape::kM % WarpShape::kM) && !(Shape::kN % WarpShape::kN),
                "Threadblock-scoped GEMM should be divisible by warp-scoped GEMM size.");

            /// Number of threads per warp
            static int const kWarpSize = warp::WarpSize<arch::OpClassSimt>::value;

            /// Number of threads total
            static int const kThreads = WarpCount::kCount * kWarpSize;

            static int const kElementsPerAccess = 1;

            //
            // Shared memory layouts
            //

            using SmemLayoutA = layout::ColumnMajor;
            using SmemLayoutB = layout::RowMajor;

            //
            // Iterators to write to shared memory
            //

            /// ThreadMap of iterator A
            using IteratorThreadMapA = transform::PitchLinearStripminedThreadMap<
                layout::PitchLinearShape<Shape::kK, Shape::kM>,
                kThreads,
                kElementsPerAccess>;

            /// Transpose the ThreadMap of iterator A
            using SmemThreadMapA = transform::TransposePitchLinearThreadMapSimt<IteratorThreadMapA>;

            /// Shared memory iterator to A operand
            using SmemIteratorA = transform::threadblock::RegularTileIterator<
                MatrixShape<Shape::kM, Shape::kK>,
                ElementA,
                SmemLayoutA,
                1,
                SmemThreadMapA // was IteratorThreadMapA
                >;

            /// ThreadMap of iterator B
            using IteratorThreadMapB = transform::PitchLinearStripminedThreadMap<
                layout::PitchLinearShape<Shape::kK, Shape::kN>,
                kThreads,
                kElementsPerAccess>;

            /// Transpose the ThreadMap of iterator A
            using SmemThreadMapB = transform::TransposePitchLinearThreadMapSimt<IteratorThreadMapB>;

            /// Shared memory iterator to B operand
            using SmemIteratorB = transform::threadblock::RegularTileIterator<
                MatrixShape<Shape::kK, Shape::kN>,
                ElementB,
                SmemLayoutB,
                0,
                SmemThreadMapB // was IteratorThreadMapA
                >;

            //
            // Warp-level matrix multiply operator
            //

            // Define the warp-level op
            static const int WarpNumThreadsM = detail::simt_get_warp_threads_m<WarpShape>();
            static const int WarpNumThreadsN = kWarpSize / WarpNumThreadsM;
            static const int ThreadTileM = WarpShape::kM / WarpNumThreadsM;
            static const int ThreadTileN = WarpShape::kN / WarpNumThreadsN;
            static_assert(!(WarpShape::kM % WarpNumThreadsM) && !(WarpShape::kN % WarpNumThreadsN),
                "WarpShape must be divisible by ThreadTile shape.");
            static const int LaneLayout = ThreadTileM > 4 && ThreadTileN > 4 ? 2 : 1;
            static const int numElementsA = 128 / sizeof_bits<ElementA>::value;
            static const int numElementsB = 1;
            static const int LaneM = cutlass::const_min(numElementsA, ThreadTileM);
            static const int LaneN = cutlass::const_min(numElementsB, ThreadTileN);

            static int const kPaddingM = detail::simt_transpose_padding(kWarpSize, Shape::kK, sizeof_bits<ElementA>::value);
            static int const kPaddingN = detail::simt_transpose_padding(kWarpSize, Shape::kK, sizeof_bits<ElementB>::value);

            static_assert(!(kPaddingM % LaneM) && !(kPaddingN % LaneN),
                "Padding must be divisible by Lane");

            // these should have max of thread tile also
            using LaneMmaShape = cutlass::gemm::GemmShape<
                LaneM,
                LaneN,
                1>;
            using Policy = cutlass::gemm::warp::MmaSimtPolicy<
                cutlass::MatrixShape<WarpNumThreadsM, WarpNumThreadsN>, // WarpShape
                cutlass::layout::RowMajorInterleaved<LaneLayout>, // LaneLayout
                LaneMmaShape>;

            using MmaWarpSimt = cutlass::gemm::warp::MmaSimt<
                WarpShape, /// Size of the Gemm problem - concept: gemm::GemmShape<> 128, 128, 8
                ElementA, /// Data type of A elements
                SmemLayoutA, /// Layout of A matrix (concept: MatrixLayout)
                ElementB, /// Data type of B elements
                SmemLayoutB, /// Layout of B matrix (concept: MatrixLayout)
                ElementC, /// Element type of C matrix
                LayoutC, /// Layout of C matrix (concept: MatrixLayout)
                Policy /// Policy describing warp-level MmaSimtOp (concept: MmaSimtOp policy)
                >;

            /// Policy used to define MmaPipelined
            using MmaPolicy = MmaPolicy<
                MmaWarpSimt,
                MatrixShape<kPaddingM, 0>, // skew for A matrix to avoid SMEM bank conflicts
                MatrixShape<0, kPaddingN>, // skew for B matrix to avoid SMEM bank conflicts
                WarpCount::kK>;
        };

    }
    namespace warp {
        /// Structure to compute the matrix product targeting CUDA cores and SIMT math instructions.
        template <
            /// Size of the Gemm problem - concept: gemm::GemmShape<>
            typename Shape_,
            /// Layout of A matrix (concept: MatrixLayout)
            typename LayoutA_,
            /// Data type of B elements
            typename LayoutB_,
            /// Layout of C matrix (concept: MatrixLayout)
            typename LayoutC_,
            /// Shape of the warp in units of thread (concept: MmaSimtPolicy)
            typename Policy_,
            /// Number of partitions along K dimension
            int PartitionsK>
        class MmaSimt<Shape_, uint8_t, LayoutA_, double, LayoutB_, double, LayoutC_, Policy_, PartitionsK> {
        public:
            /// Shape of warp-level matrix operation (concept: GemmShape)
            using Shape = Shape_;

            /// Data type of multiplicand A
            using ElementA = uint8_t;

            /// Layout of multiplicand A
            using LayoutA = LayoutA_;

            /// Data type of multiplicand B
            using ElementB = double;

            /// Layout of multiplicand B
            using LayoutB = LayoutB_;

            /// Data type of accumulator matrix C
            using ElementC = double;

            /// Layout of accumulator matrix C
            using LayoutC = LayoutC_;

            /// Shape of the warp in units of thread (concept: MmaLanePolicySimt)
            using Policy = Policy_;

            /// Indicates class of matrix operator
            using OperatorClass = arch::OpClassSimt;

            /// Hard-coded for now
            using ArchTag = arch::Sm61;

            /// Layout of threads
            using ThreadLayoutA = cutlass::layout::RowMajor;

            using ThreadLayoutB = cutlass::layout::ColumnMajor;

            /// Thread-level matrix multiply accumulate operator
            using ThreadMma = thread::Mma<
                GemmShape<
                    Shape::kM / Policy::WarpShape::kRow,
                    Shape::kN / Policy::WarpShape::kColumn,
                    Policy::LaneMmaShape::kK>,
                ElementA,
                ThreadLayoutA,
                ElementB,
                ThreadLayoutB,
                ElementC,
                LayoutC,
                arch::OpMultiplyAdd,
                bool>;

            /// Underlying matrix multiply operator (concept: arch::Mma)
            using ArchMmaOperator = typename ThreadMma::ArchMmaOperator;

            /// Indicates math operator
            using MathOperator = typename ArchMmaOperator::Operator;

            /// Shape of the underlying instruction
            using InstructionShape = cutlass::gemm::GemmShape<1, 1, 1>;

        public:
            /// Iterates over the A operand in memory
            using IteratorA = MmaSimtTileIterator<
                MatrixShape<Shape::kM, Policy::LaneMmaShape::kK>,
                Operand::kA,
                ElementA,
                LayoutA,
                Policy,
                PartitionsK,
                Shape::kK>;

            /// Storage for A tile
            using FragmentA = typename IteratorA::Fragment;

            /// Storage for transformed A tile
            using TransformedFragmentA = FragmentA;

            /// Iterates over the B operand in memory
            using IteratorB = MmaSimtTileIterator<
                MatrixShape<Policy::LaneMmaShape::kK, Shape::kN>,
                Operand::kB,
                ElementB,
                LayoutB,
                Policy,
                PartitionsK,
                Shape::kK>;

            /// Storage for B tile
            using FragmentB = typename IteratorB::Fragment;

            /// Storage for transformed A tile
            using TransformedFragmentB = FragmentB;

            /// Iterates over the C operand in memory
            using IteratorC = MmaSimtTileIterator<
                MatrixShape<Shape::kM, Shape::kN>,
                Operand::kC,
                ElementC,
                LayoutC,
                Policy>;

            /// Storage for C tile
            using FragmentC = typename ThreadMma::FragmentC;

        public:
            //
            // Methods
            //

            /// Ctor
            CUTLASS_DEVICE
            MmaSimt() { }

            /// Performs a warp-level matrix multiply-accumulate operation
            CUTLASS_DEVICE
            void operator()(
                FragmentC& d,
                FragmentA a,
                FragmentB b,
                FragmentC const& c, int group_idx = 0) const
            {

                ThreadMma mma;

                mma(d, a, b, c);
            }
        };

    }
    namespace thread {

        template <
            /// Size of the Gemm problem - concept: gemm::GemmShape<>
            typename Shape_,
            /// Layout of C matrix (concept: MatrixLayout)
            typename LayoutC_>
        struct Mma<
            Shape_,
            uint8_t,
            cutlass::layout::RowMajor,
            double,
            cutlass::layout::ColumnMajor,
            double,
            LayoutC_,
            arch::OpMultiplyAdd,
            bool> {

            /// Size of the Gemm problem - concept: gemm::GemmShape<>
            using Shape = Shape_;

            /// Data type of operand A
            using ElementA = uint8_t;

            /// Layout of A matrix (concept: layout::MapFunc)
            using LayoutA = cutlass::layout::RowMajor;

            /// Data type of operand B
            using ElementB = double;

            /// Layout of B matrix (concept: layout::MapFunc)
            using LayoutB = cutlass::layout::ColumnMajor;

            /// Element type of operand C
            using ElementC = double;

            /// Layout of C matrix (concept: layout::MapFunc)
            using LayoutC = LayoutC_;

            /// Underlying mathematical operator
            using Operator = arch::OpMultiplyAdd;

            /// A operand storage
            using FragmentA = Array<ElementA, Shape::kMK>;

            /// B operand storage
            using FragmentB = Array<ElementB, Shape::kKN>;

            /// C operand storage
            using FragmentC = Array<ElementC, Shape::kMN>;

            /// Underlying matrix multiply operator (concept: arch::Mma)
            //  Use 1x1x4 IDP4A sequence for bulk of computation
            using ArchMmaOperator = arch::Mma<
                cutlass::gemm::GemmShape<1, 1, 1>,
                1,
                ElementA,
                LayoutA,
                ElementB,
                LayoutB,
                ElementC,
                LayoutC,
                arch::OpMultiplyAdd>;

            //
            // Methods
            //

            /// Computes a matrix product D = A * B + C
            CUTLASS_HOST_DEVICE
            void operator()(
                FragmentC& D,
                FragmentA const& A,
                FragmentB const& B,
                FragmentC const& C)
            {

                TensorRef<ElementC, LayoutC> d(
                    reinterpret_cast<ElementC*>(&D), LayoutC::packed({ Shape::kM, Shape::kN }));

                // Copy accumulators
                D = C;

                /// Use 1x1x4 IDP4A sequence for bulk of computation
                ArchMmaOperator mma;

                // Compute matrix product
                CUTLASS_PRAGMA_UNROLL
                for (int k = 0; k < Shape::kK / ArchMmaOperator::Shape::kK; ++k) {

                    CUTLASS_PRAGMA_UNROLL
                    for (int n = 0; n < Shape::kN; ++n) {

                        CUTLASS_PRAGMA_UNROLL
                        for (int m = 0; m < Shape::kM; ++m) {
                            MatrixCoord mn(m, n);

                            Array<ElementA, 4> const* ptr_A = reinterpret_cast<Array<ElementA, 4> const*>(&A);
                            Array<ElementB, 4> const* ptr_B = reinterpret_cast<Array<ElementB, 4> const*>(&B);

                            Array<ElementC, 1> tmp = reinterpret_cast<Array<ElementC, 1>&>(d.at(mn));

                            Array<ElementA, 4> tmpArray = ptr_A[0];
                            ElementA tmpA = tmpArray.at(0);
                            // printf("(%d, %d, %d) ", tmpA == ElementA(0), n, m);

                            mma(
                                tmp,
                                ptr_A[m * Shape::kK / ArchMmaOperator::Shape::kK + k],
                                ptr_B[n * Shape::kK / ArchMmaOperator::Shape::kK + k],
                                tmp);

                            d.at(mn) = reinterpret_cast<ElementC&>(tmp);
                        }
                    }
                }
            }
        };
    }
}
}
