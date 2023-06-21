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
using CSV;
using DataFrames

# =====================
# Global definitions
# =====================

SIZE="xsmall"

ROOT_DIR = string(@__DIR__) * "/../.."
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"

DATA_DIR = ROOT_DIR * "/data"
DATA_FILE = DATA_DIR * "/$SIZE.bed"
FREQ_FILE = DATA_DIR * "/$SIZE.freq"
LD_FILE = DATA_DIR * "/$SIZE.ld"

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

if !isfile(LD_FILE)
    cd(DATA_DIR)
    run(`./plink --bfile $SIZE --r square`)
    run(`mv plink.ld $LD_FILE`)
    cd(ROOT_DIR)
end

## Test LD calculation
@time genotype_data, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE, coding_twobit = true, calc_freq = true)

@testset "Correctness" begin
    @time M = miraculix.crossproduct.ld(genotype_data, n_snps, n_indiv, is_plink_format = false, allele_freq = freq)
    
    ld_plink = Matrix(CSV.read(LD_FILE, delim = '\t', header = 0, DataFrame))

    @test isequal(size(ld_plink),size(M))
    @test sum(abs.(ld_plink - M)) < 0.1
end
