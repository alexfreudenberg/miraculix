
/*
 Authors
 Alexander Freudenberg


 Copyright (C) 2022 Alexander Freudenberg

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

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <string>

// CUTLASS includes
#include "cutlass/arch/mma.h"
#include "cutlass/arch/wmma.h"
#include "cutlass/array.h"
#include "cutlass/cutlass.h"
#include "cutlass/gemm/device/gemm.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/layout/matrix.h"
#include "cutlass/numeric_types.h"
#include "cutlass/util/host_tensor.h"
#include "cutlass/util/reference/device/gemm.h"
#include "cutlass/util/reference/host/gemm.h"
#include "cutlass/util/reference/host/tensor_compare.h"
#include "cutlass/util/reference/host/tensor_copy.h"
#include "cutlass/util/reference/host/tensor_fill.h"
#include "cutlass/util/reference/host/tensor_reduce.h"
#include "cutlass/util/tensor_view_io.h"

#include "matvec_u8x4f64.h"
#include "reference_u2d.h"


#ifdef STANDALONE
int main()
{
    using ColumnMajor = cutlass::layout::ColumnMajor;
    using RowMajor = cutlass::layout::RowMajor;

    using ElementA_ = cutlass::u4f64_t;
    
    using LayoutA_ = RowMajor;
    
    const int length_m = 64;
    const int length_n = 64;
    const int length_k = 64;
    cutlass::gemm::GemmCoord problem_size(length_m, length_n, length_k);

    // Initialize tensors using CUTLASS helper functions
    cutlass::HostTensor<ElementA_, LayoutA_> tensor_a(
        problem_size.mk());
  
    for(int i = 0; i < 4; i++) {
        cutlass::u4f64_t* tmp = tensor_a.host_data(); 
        tmp[0].storage[i] = (double) i;
    }
     for(int i = 0; i < 4; i++) {
        cutlass::u4f64_t* tmp = tensor_a.host_data(); 
        printf("%.2lf %.2lf\n", tmp[0].storage[i], (double) i);
    }
    printf("\n");
    return 0;
}
#endif