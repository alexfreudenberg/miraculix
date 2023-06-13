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
using Base;
using Random;
using Distances; 
using SparseArrays;
using LinearAlgebra;
using Test
using .Threads: @threads

# =====================
# Global definitions
# =====================

ROOT_DIR = string(@__DIR__) * "/../.."
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"

tol = 1e-1;
Random.seed!(0);

# Remove commit message verbosity
ENV["PRINT_LEVEL"] = "1";

include(MODULE_PATH)

# =====================
# Auxiliary functions
# =====================

function pack_twobit(::Type{T}, M::Matrix{T}, n_snps, n_indiv) where{T}
    @assert (n_indiv, n_snps) == size(M)
    BITS = sizeof(T) * 8;
    CODES_PER_UNIT = Int(BITS/2); 
    n_vec = Int(ceil(n_indiv*2/BITS));
    M_packed = zeros(T, (n_vec, n_snps)); # packed copy of M
    # Pack M into M_packed
    @inbounds @views @threads for i = 1:n_snps
        for j = 1:n_vec
            for k = 1:min(CODES_PER_UNIT, n_indiv - (j-1) * CODES_PER_UNIT) # counts from 0 to minimum of CODES_PER_UNIT and remaining number of rows
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

# =====================
# Main
# =====================

println("Load library and set options")
miraculix.set_library_path(LIBRARY_PATH)
miraculix.load_shared_library()


## Test GRM functionality
T = UInt8;
n_snps = 2048 * 4 + 423; # Number of SNPs
n_indiv = 2048 * 2 + 103; # Number of individuals 

M = rand(Vector{T}(0:2), (n_snps, n_indiv));
M = ones(T, (n_snps, n_indiv));

M_packed = pack_twobit(T, M, n_indiv, n_snps);
println(size(M_packed))

ANS = miraculix.grm.compute(M_packed, n_snps, n_indiv)

D = BLAS.gemm('T','N',Matrix{Float64}(M),Matrix{Float64}(M));
deviation = sum(abs.(ANS-D))
println(deviation)
println(ANS[1:10])
println(D[1:10])
@test deviation<1e-4