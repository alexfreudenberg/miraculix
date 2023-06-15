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
DATA_FILE = ROOT_DIR * "/data/xsmall.bed"
FREQ_FILE = ROOT_DIR * "/data/xsmall.freq"

tol = 1e-1;
Random.seed!(0);

# Remove commit message verbosity
ENV["PRINT_LEVEL"] = "1";

include(MODULE_PATH)



# =====================
# Main
# =====================

println("Load library and set options")
miraculix.set_library_path(LIBRARY_PATH)
miraculix.load_shared_library()


## Test GRM functionality
T = UInt8;

@time genotype_data, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE,coding_twobit = true, calc_freq = true)

#freq = miraculix.read_plink.read_freq(FREQ_FILE)
println(size(freq))

M = miraculix.grm.compute(genotype_data, n_snps, n_indiv; do_center = false, is_transposed = false, scaling_method = 0, is_plink_format = false, allele_freq = freq)

BLAS.syr!('U', Float64( -4*n_indiv), freq, M)
M = Symmetric(M, :U)

M_unpacked = miraculix.dgemm_compressed.decompress_2bit_format(genotype_data, n_indiv, n_snps)
M_centered = M_unpacked .- sum(M_unpacked, dims = 1)/size(M_unpacked,1)
Cov = transpose(M_centered) * M_centered