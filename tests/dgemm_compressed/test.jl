
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

# This is a test file for the miraculix.jl module that provides an 
# interface to the miraculix shared library, a tool that enables 
# highly-efficient computations on compressed SNP data stored in PLINK
# format. 

# The functions being tested perform operations such as setting computation
# options, preprocessing and transposing a PLINK .bed matrix, performing 
# matrix multiplication, and freeing memory resources.
# The tests in this file aim to ensure the consistency and correctness of
# these functions.

using Statistics
using Libdl
using Test

ROOT_DIR = string(@__DIR__) * "/../.."
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_FILE = ROOT_DIR * "/data/xsmall.bed"
FREQ_FILE = ROOT_DIR * "/data/xsmall.freq"

include(MODULE_PATH)


println("Load library and set options")
miraculix.dgemm_compressed.set_library_path(LIBRARY_PATH)
miraculix.dgemm_compressed.load_shared_library()
miraculix.dgemm_compressed.set_options(use_gpu=!false)

println("Read bed file and frequencies")
genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)

println("Transpose matrix")
genotype_data_transposed = miraculix.dgemm_compressed.transpose_genotype_matrix(genotype_data, n_snps, n_indiv)

@testset "Consistency" begin
    @test size(genotype_data) == (ceil(n_indiv/4), n_snps)
    @test size(genotype_data_transposed) == (ceil(n_snps/4), n_indiv)
    @test length(freq) == n_snps
end

println("Decompress genotype data")
genotype_data_decompressed = miraculix.dgemm_compressed.decompress_genotype_data(genotype_data, n_indiv, n_snps)
genotype_data_t_decompressed = miraculix.dgemm_compressed.decompress_genotype_data(genotype_data_transposed, n_snps, n_indiv)

@testset "Transpose operation" begin
    @test size(genotype_data_decompressed) == (n_indiv, n_snps)
    @test size(genotype_data_decompressed) == size(transpose(genotype_data_t_decompressed))
    @test all(genotype_data_decompressed == transpose(genotype_data_t_decompressed))
end

println("Initialize dgemm_compressed routine")
obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 10)

n_col = 10
B = randn(Float64, n_snps, n_col)
B_T = randn(Float64, n_indiv, n_col)

println("Calculate genotype matrix multiplication")
C = miraculix.dgemm_compressed.dgemm_compressed_main(false, obj_ref, B, n_snps, n_indiv)
C_T = miraculix.dgemm_compressed.dgemm_compressed_main(true, obj_ref, B_T, n_snps, n_indiv)

println("Test results")
freq_test = mean(genotype_data_decompressed ./ 2.0, dims=1)
C_test = (genotype_data_decompressed .- 2.0 * freq_test) * B
C_t_test = transpose(genotype_data_decompressed .- 2.0 * freq_test) * B_T

@testset "Correctness of results" begin
    @test maximum(abs.(vec(freq_test) .- freq) ./ freq_test) < 1e-2
    @test maximum(abs.(C_test .- C)) < 1e-1
    @test maximum(abs.(C_t_test .- C_T)) < 1e-1
end

miraculix.dgemm_compressed.free_compressed(obj_ref)