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
@timev "Loading libraries" begin
    using Random;
    using LinearAlgebra;
    using CSV;
    using ClusterManagers, Distributed;
    using DelimitedFiles;
    using DataFrames;
    using BenchmarkTools;
    using Statistics;
    using Libdl;
end

# =====================
# Global definitions
# =====================

ROOT_DIR = string(@__DIR__) * "/../.."

MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_DIR = ROOT_DIR * "/data"
LOG_DIR = DATA_DIR * "/logs"

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



# =====================
# Main
# =====================

miraculix.set_library_path(LIBRARY_PATH)
miraculix.load_shared_library()
miraculix.dgemm_compressed.set_options(use_gpu=true, verbose=0)

init_sym = dlsym(miraculix.LIBRARY_HANDLE[], :plink2compressed)
dgemm_compressed_sym = dlsym(miraculix.LIBRARY_HANDLE[], :dgemm_compressed)
free_sym = dlsym(miraculix.LIBRARY_HANDLE[], :free_compressed)


# Set hyper parameters
max_iter   = 1_000 # Maximum number of iterations
print_iter = 1e2 # Frequency of convergence information
conv_crit  = 1e-2 # Maximum norm of residual
n_devices  = 2 # Number of devices
n_traits   = 5 # Number of traits to solve for

# We assume that genotype data has been generated by the R package MoBPS 
data_file = DATA_DIR * "/small.bed"

# Read-in data from PLINK binary format
@info "Reading in data from $data_file and transpose it"
@timev "Preprocessing" begin
    # Read PLINK data and calculate allele frequencies
    wtime = @elapsed plink, freq, n_snps, n_indiv = miraculix.read_plink.read_bed(data_file, coding_twobit = true, calc_freq = true, check_for_missings = false)
    @debug "Time for reading: $wtime s."

    if (length(ARGS) > 0) && (ARGS[1] == "test")
        n_snps = 1000
        plink = plink[:,1:n_snps]
        freq = freq[1:n_snps]
    end
    
    # Transpose matrix
    wtime = @elapsed plink_transposed = miraculix.compressed_operations.transpose_genotype_matrix(plink, n_snps, n_indiv)
    @debug "Time for transposing: $wtime s."

    GC.gc()
end


obj_ref = Ref{Ptr{Cvoid}}(C_NULL)
ENV["CUDA_DEVICE"] = 0
ccall(init_sym,  Cvoid,  (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint, Ptr{Ptr{Cvoid}}), plink, C_NULL, Int32(n_snps), Int32(n_indiv), freq, Int32(n_traits), obj_ref)

obj_ref_trans = Ref{Ptr{Cvoid}}(C_NULL)
ENV["CUDA_DEVICE"] = 1
ccall(init_sym,  Cvoid,  (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint, Ptr{Ptr{Cvoid}}), C_NULL, plink_transposed, Int32(n_snps), Int32(n_indiv), freq, Int32(n_traits), obj_ref_trans)


ccall(free_sym, Cvoid, (Ptr{Ptr{Cvoid}},), obj_ref)
ccall(free_sym, Cvoid, (Ptr{Ptr{Cvoid}},), obj_ref_trans)
