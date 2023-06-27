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
using LinearAlgebra;
using CSV;
using DelimitedFiles;
using DataFrames;
using BenchmarkTools;
using Tables;

# =====================
# Global definitions
# =====================

ROOT_DIR = string(@__DIR__) * "/../.."

MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_DIR = ROOT_DIR * "/data"
LOG_DIR = DATA_DIR * "/logs"

BENCHMARK_SIZES_GRM=["xsmall","few_snps", "medium_snps", "many_snps"]
BENCHMARK_SIZES_LD=["xsmall","small", "medium", "large"]

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1_000_000

# Control miraculix verbosity
ENV["PRINT_LEVEL"] = "1";

OMP_NUM_THREADS = ENV["OMP_NUM_THREADS"];
println("OMP threads set to $OMP_NUM_THREADS")

include(MODULE_PATH)

# =====================
# Auxiliary functions
# =====================
function write_result(root_file_name::String, matrix::Matrix{Float64},write_format::String)
    if write_format == "tsv"
        CSV.write(root_file_name * ".grm.tsv", Tables.table(matrix), delim='\t')
    elseif write_format == "binary"
        n_bytes = size(matrix,1)^2 * 8
        open(root_file_name * ".grm.bin", "w") do io
            unsafe_write(io, pointer(matrix), n_bytes)
        end
    else 
        error("Unsupported write format supplied.")
    end
end

function run_miraculix_grm(data::String, write_format::String = "binary")
    # Create valid bed file from data string
    data_file = data * ".bed"
    # Read-in data file, convert it to two-bit and calculate allele frequencies
    @time plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = true, calc_freq = true)

    # Transpose genotype data to sample-major format for GRM calculation
    @time plink_transposed = miraculix.compressed_operations.transpose_genotype_matrix(plink, n_snps, n_indiv)

    # Calculate GRM matrix
    G1 = miraculix.crossproduct.grm(plink_transposed, n_snps, n_indiv, is_plink_format = false, allele_freq = vec(freq), do_scale = false)
    # Scale GRM matrix analoguous to PLINK
    G1 ./= n_snps

    # Write results to file 
    @time write_result(data, G1, write_format)
    return Nothing
end
function run_miraculix_ld(data::String, write_format::String = "binary")
    # Create valid bed file from data string
    data_file = data * ".bed"
    # Read-in data file, convert it to two-bit and calculate allele frequencies
    @time plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = true, calc_freq = true)

    # Calculate LD matrix
    M = miraculix.crossproduct.ld(plink, n_snps, n_indiv, is_plink_format = false, allele_freq = freq)

    # Write results to file
    @time write_result(data, M, write_format)
    # M[diagind(M)] .= 0;
    # ld_score = sum(M.^2, dims = 1)
    # ld_max = maximum(M.^2, dims = 1)

    return Nothing
end

function run_gcta_grm(data::String)
    # Run GCTA software from command line
    run(`./gcta-1.94.1 --bfile $data --thread-num $OMP_NUM_THREADS --make-grm-bin --make-grm-alg 1 --out $data `)
end
function run_gcta_ld(data::String)
    # Run GCTA software from command line
    run(`./gcta-1.94.1 --bfile $data --thread-num $OMP_NUM_THREADS --ld-score --out $data`)
end
function run_plink_grm(data::String)
    # Run PLINK software from command line
    run(`./plink --bfile $data --threads $OMP_NUM_THREADS --make-rel square cov`)
end
function run_plink_ld(data::String)
    # Run PLINK software from command line
    run(`./plink --bfile $data --threads $OMP_NUM_THREADS --r square`)
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
suite["GRM"] = BenchmarkGroup(["GRM", "crossproduct"])
suite["LD"] = BenchmarkGroup(["LD", "crossproduct"])

for size in BENCHMARK_SIZES_GRM
    suite["GRM"][size,"miraculix"] = @benchmarkable run_miraculix_grm($size) setup = (run_miraculix_grm($size))
    suite["GRM"][size,"PLINK"] = @benchmarkable run_plink_grm($size)
    suite["GRM"][size,"GCTA"] = @benchmarkable run_gcta_grm($size)
end
for size in BENCHMARK_SIZES_LD
    suite["LD"][size,"miraculix"] = @benchmarkable run_miraculix_ld($size) setup = (run_miraculix_ld($size))
    suite["LD"][size,"PLINK"] = @benchmarkable run_plink_ld($size)
    suite["LD"][size,"GCTA"] = @benchmarkable run_gcta_ld($size)
end
