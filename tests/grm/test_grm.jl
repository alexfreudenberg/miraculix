#  Authors 
#  Alexander Freudenberg, alexander.freudenberg@stads.de

#  Copyright (C) 2022-2023 Alexander Freudenberg

#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

# =====================================================
# This file is heavily WIP
# =====================================================
using Base, Base.Libc.Libdl;
using Distances; 
using SparseArrays;
using LinearAlgebra;
using Test
using .Threads: @threads


## Test GRM functionality
T = UInt32;
N_ROW = 256 * 100; # Number of SNPs
N_COL = 10_000; # Number of individuals

M = rand(Vector{UInt32}(0:2), (N_ROW,N_COL));
ANS = zeros(Float64, (N_COL, N_COL));

function pack_twobit(::Type{T}, M::Matrix{T}, n_row, n_col) where{T}
    BITS = sizeof(T) * 8;
    CODES_PER_UNIT = Int(BITS/2); 
    n_vec = Int(ceil(n_row*2/BITS));
    M_packed = zeros(UInt32, (n_vec, n_col)); # packed copy of M
    
    # Pack M into M_packed
    @inbounds @views @threads for i = 1:n_col
        for j = 1:n_vec
            for k = 1:min(CODES_PER_UNIT, n_row - (j-1) * CODES_PER_UNIT) # counts from 0 to minimum of CODES_PER_UNIT and remaining number of rows
                M_packed[j,i] |= (M[ (j-1) * CODES_PER_UNIT + k ,i]  << (2 * (k-1))); # consecutively shift CODES_PER_UNIT entries of M by two bits and OR them into M_packed
            end
        end
    end
    return M_packed
end
function crossprod(::Type{T}, M::Matrix, n, k) where{T}
    VEC = zeros(T,n);
    Mi = zeros(T,n);
    COUNTS = zeros(UInt32, n);
    RESULT = zeros(UInt32,(k,k));
    @inbounds @views @threads for i = 1:k
        copyto!(Mi, view(M,:, i));
        @inbounds for j in i:k
            broadcast!(&, VEC, Mi, view(M, :, j));
            broadcast!(count_ones, COUNTS, VEC);
            RESULT[i,j] = RESULT[j,i] = sum(COUNTS);
        end
    end
    return RESULT;
end
@btime M_packed = pack_twobit(T, M, N_ROW, N_COL);
M_packed = pack_twobit(T, M, N_ROW, N_COL);
@btime crossprod(T, M_packed, size(M_packed,1), N_COL);

 ccall((:crossprod_mmagpu, "./src/miraculixjl.so"),Cvoid,(Ptr{Cuint}, Cuint, Cuint, Ptr{Cdouble}), M_packed, UInt32(N_ROW), UInt32(N_COL), ANS);

D = BLAS.gemm('T','N',Matrix{Float32}(M),Matrix{Float32}(M));
@test sum(ANS-D)<1e-4