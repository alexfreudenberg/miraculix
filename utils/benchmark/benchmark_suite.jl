#  Authors 
#  Alexander Freudenberg, alexander.freudenberg@stads.de

#  Copyright (C) 2023 Alexander Freudenberg

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
using Libdl;

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

# Get thread number
@assert haskey(ENV, "OMP_NUM_THREADS") "OMP_NUM_THREADS not set"
OMP_NUM_THREADS = ENV["OMP_NUM_THREADS"];
BLAS.set_num_threads(parse(Int,OMP_NUM_THREADS))
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
    @time "Reading data" plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = true, calc_freq = true, check_for_missings = false)

    # Transpose genotype data to sample-major format for GRM calculation
    @time "Transpose matrix" plink_transposed = miraculix.compressed_operations.transpose_genotype_matrix(plink, n_snps, n_indiv)

    # Calculate GRM matrix
    @time "Calculate GRM" G1 = miraculix.crossproduct.grm(plink_transposed, n_snps, n_indiv, is_plink_format = false, allele_freq = vec(freq), do_scale = false)
    # Scale GRM matrix analoguous to PLINK
    G1 ./= n_snps

    # Write results to file 
    @time "Writing result" write_result(data, G1, write_format)
    
    GC.gc()
    return Nothing
end
function run_miraculix_ld(data::String, write_format::String = "binary")
    # Create valid bed file from data string
    data_file = data * ".bed"
    # Read-in data file, convert it to two-bit and calculate allele frequencies
    @time "Reading data" plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = true, calc_freq = true)

    # Calculate LD matrix
    @time "Calculating LD" M = miraculix.crossproduct.ld(plink, n_snps, n_indiv, is_plink_format = false, allele_freq = freq)

    # Write results to file
    @time "Writing result" write_result(data, M, write_format)

    GC.gc()
    return M
end

function run_cublas_uint8_grm(data::String, libpath::String)
    # Create valid bed file from data string
    data_file = data * ".bed"
    # Read-in data file, convert it to two-bit and calculate allele frequencies
    wtime = @elapsed plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = false, calc_freq = true, check_for_missings = false)
    @debug "Time for reading data: $wtime s."

    if ~isfile(libpath)
        error("cublas_uint8 library not available.")
    end
    lib_handle = dlopen(libpath)
    compute_sym = dlsym(lib_handle, :cublas_uint8_gemm)

    # Check if memory suffices
    device_memory = 80 # in GB
    if n_indiv^2 * 4 + n_indiv * n_snps > device_memory * 1024^3 * 0.9
        @warn "Device memory not sufficient for $data."
        return Nothing
    end
    # Decompress data into uint8
    decompressed = zeros(UInt8, n_indiv, n_snps)
    for (index,entry) in pairs(IndexLinear(),plink)
        offset_decompressed = (index-1) * 4 + 1
        @inbounds for i in 0:3
            # Convert packed SNP data to Float
            genotype_uint8 = UInt8((entry >> (2*i)) & 0x03)
            # Check if there is a missing value which is coded as 1 in PLINK format
            (genotype_uint8 == UInt8(1)) && error("Missing in genotype data")
            # Convert PLINK format to 0,1,2 format
            decompressed[offset_decompressed + i] = max(0, genotype_uint8 -1) |> UInt8
        end    
    end

    M = zeros(Float64, (n_indiv, n_indiv))
    wtime = @elapsed ccall(compute_sym,  Cint,  (Ptr{UInt8}, Cint, Cint, Ptr{Float64}), decompressed, Int32(n_snps), Int32(n_indiv), M)
    @debug "Time for calculating cuBLAS uint8 crossproduct: $wtime s."

    dlclose(lib_handle)
    @assert issymmetric(M) "Result not symmetric" 

    # Scaling of centered genotype matrix
    col_sum = sum(M, dims = 1) |> vec
    one_vector = ones(Float64, (n_indiv,)) 

    wtime = @elapsed begin 
        BLAS.ger!(-1/n_indiv, col_sum, one_vector, M)
        BLAS.ger!(-1/n_indiv, one_vector, col_sum, M)
    end
    @debug "Time for rank-one updates of GRM: $wtime s."

    wtime = @elapsed M .+= sum(col_sum) ./ n_indiv^2
    @debug "Time for affine transform of GRM: $wtime s."
    
    M ./= n_snps

    wtime = @elapsed "Writing result" write_result(data, M, "binary")
    @debug "Time for writing result: $wtime s."

    return M
end # function
function run_cublas_uint8_ld(data::String, libpath::String)
    # Create valid bed file from data string
    data_file = data * ".bed"
    # Read-in data file, convert it to two-bit and calculate allele frequencies
    wtime = @elapsed plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = false, calc_freq = true, check_for_missings = false)
    @debug "Time for reading data: $wtime s."

    if ~isfile(libpath)
        error("cublas_uint8 library not available.")
    end
    lib_handle = dlopen(libpath)
    compute_sym = dlsym(lib_handle, :cublas_uint8_gemm)

    # Check if memory suffices
    device_memory = 80 # in GB
    if n_snps^2 * 4 + n_indiv * n_snps > device_memory * 1024^3 * 0.9
        @warn "Device memory not sufficient for $data."
        return Nothing
    end
    # Decompress data into uint8
    decompressed = zeros(UInt8, n_indiv, n_snps)
    for (index,entry) in pairs(IndexLinear(),plink)
        offset_decompressed = (index-1) * 4 + 1
        @inbounds for i in 0:3
            # Convert packed SNP data to Float
            genotype_uint8 = UInt8((entry >> (2*i)) & 0x03)
            # Check if there is a missing value which is coded as 1 in PLINK format
            (genotype_uint8 == UInt8(1)) && error("Missing in genotype data")
            # Convert PLINK format to 0,1,2 format
            decompressed[offset_decompressed + i] = max(0, genotype_uint8 -1) |> UInt8
        end    
    end
    decompressed = transpose(decompressed) |> Matrix;

    M = zeros(Float64, (n_snps, n_snps))
    wtime = @elapsed ccall(compute_sym,  Cint,  (Ptr{UInt8}, Cint, Cint, Ptr{Float64}), decompressed, Int32(n_indiv), Int32(n_snps), M)
    @debug "Time for calculating cuBLAS uint8 crossproduct: $wtime s."

    dlclose(lib_handle)
    @assert issymmetric(M) "Result not symmetric" 

     # Scaling of centered genotype matrix
     wtime = @elapsed begin
        BLAS.syr!('U', Float64(-4.0 * n_indiv), freq, M)
        M = Matrix(Symmetric(M, :U))
    end
    @debug "Time for rank-one update of LD: $wtime s."

    # Calculate vector of standard deviations
    wtime = @elapsed begin
        sigma_vector = reshape(sqrt.(diag(M)), (n_snps,1))
        # Devide each row and column by the vector of standard deviations
        M ./= sigma_vector
        M ./= transpose(sigma_vector)
    end
    @debug "Time for scaling LD by std devs: $wtime s."
    
    wtime = @elapsed  "Writing result" write_result(data, M, "binary")
    @debug "Time for writing values: $wtime s."

    return M
end # function

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

for problem_size in BENCHMARK_SIZES_GRM
    suite["GRM"][problem_size,"miraculix"] = @benchmarkable run_miraculix_grm($problem_size) setup = (run_miraculix_grm($problem_size))
    suite["GRM"][problem_size,"PLINK"] = @benchmarkable run_plink_grm($problem_size)
    suite["GRM"][problem_size,"GCTA"] = @benchmarkable run_gcta_grm($problem_size)
    suite["GRM"][problem_size,"cuBLAS"] = @benchmarkable run_cublas_uint8_grm($problem_size, ROOT_DIR * "/utils/benchmark/cublas_uint8.so")
end
for problem_size in BENCHMARK_SIZES_LD
    suite["LD"][problem_size,"miraculix"] = @benchmarkable run_miraculix_ld($problem_size) setup = (run_miraculix_ld($problem_size))
    suite["LD"][problem_size,"PLINK"] = @benchmarkable run_plink_ld($problem_size)
    suite["LD"][problem_size,"GCTA"] = @benchmarkable run_gcta_ld($problem_size)
end

suite["LD"]["small","cuBLAS"] = @benchmarkable run_cublas_uint8_ld("small", ROOT_DIR * "/utils/benchmark/cublas_uint8.so")