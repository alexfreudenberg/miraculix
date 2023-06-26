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

using Base;
using Random;
using Distances; 
using LinearAlgebra;
using MKL;
using Test;
using CSV;
using DataFrames;
using BenchmarkTools;

# =====================
# Global definitions
# =====================

ROOT_DIR = string(@__DIR__) * "/../.."

MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_DIR = ROOT_DIR * "/data"

BENCHMARK_SIZES=["few_snps", "medium_snps", "many_snps"]

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1_000_000

# Remove commit message verbosity
ENV["PRINT_LEVEL"] = "1";

OMP_NUM_THREADS = ENV["OMP_NUM_THREADS"];
println("OMP threads set to $OMP_NUM_THREADS")

include(MODULE_PATH)

# =====================
# Auxiliary function
# =====================
function run_miraculix_grm(data::String)
    data_file = data * ".bed"
    plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = true, calc_freq = true)

    plink_transposed = miraculix.compressed_operations.transpose_genotype_matrix(plink, n_snps, n_indiv)

    G1 = miraculix.crossproduct.grm(plink_transposed, n_snps, n_indiv, is_plink_format = false, allele_freq = vec(freq), do_scale = false)
    G1 ./= n_snps
    CSV.write(G1, delim= '\t')
    return Nothing
end

function run_plink_grm(data::String)
    run(`./plink --bfile $data --threads $OMP_NUM_THREADS --make-rel square cov`)
end
# =====================
# Main
# =====================

println("Load library and set options")
miraculix.set_library_path(LIBRARY_PATH)
miraculix.load_shared_library()
cd(DATA_DIR)

## Benchmark
suite = BenchmarkGroup()
suite["GRM"] =  BenchmarkGroup(["GRM", "crossproduct"])
suite["LD"] =  BenchmarkGroup(["LD", "crossproduct"])


for size in BENCHMARK_SIZES
    suite["GRM"][size, "miraculix"] = @benchmarkable run_miraculix_grm($size)
    suite["GRM"][size, "plink"] = @benchmarkable run_plink_grm($size)
end

for size in BENCHMARK_SIZES
    run_miraculix_grm(size)
    run_plink_grm(size)
end
