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
using MKL;
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
BLAS.set_num_threads(parse(Int,ENV["OMP_NUM_THREADS"]))

include(MODULE_PATH)

# =====================
# Auxiliary functions
# =====================

function pack_twobit(::Type{T}, M::Matrix{T}, n_snps::Int64, n_indiv::Int64) where{T}
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
println("Check if routine returns right results")

@testset "Correctness" begin
    for n_snps in Vector{Int64}([1e4, 5e4, 1e5])
        for n_indiv in Vector{Int64}([2e3, 15e3])
            println("n_snps: ",n_snps, " n_indiv: ", n_indiv)
            println("Simulation")
            @time M = rand(Vector{T}(0:2), (n_snps, n_indiv));
            println("Packing")
            @time M_packed = pack_twobit(T, M, n_indiv, n_snps);
            println("Casting")
            @time M_double = Matrix{Float64}(M)

            ANS = miraculix.grm.compute(M_packed, n_snps, n_indiv)
            println("MMAGPU:")
            @time miraculix.grm.compute(M_packed, n_snps, n_indiv)
            println("DGEMM:")
            @time D = BLAS.gemm('T','N',M_double,M_double);
            deviation = sum(abs.(ANS-D))
            @test deviation<1e-4
        end
    end
end 


@testset "Correctness in uneven dimensions" begin
    for n_snps in Vector{Int64}([953, 10251])
        for n_indiv in Vector{Int64}([752, 5343, 12433])
            M = rand(Vector{T}(0:2), (n_snps, n_indiv));
            M_packed = pack_twobit(T, M, n_indiv, n_snps);

            ANS = miraculix.grm.compute(M_packed, n_snps, n_indiv)
            @time miraculix.grm.compute(M_packed, n_snps, n_indiv)
            D = BLAS.gemm('T','N',Matrix{Float64}(M),Matrix{Float64}(M));
            deviation = sum(abs.(ANS-D))
            @test deviation<1e-4
        end
    end
end 