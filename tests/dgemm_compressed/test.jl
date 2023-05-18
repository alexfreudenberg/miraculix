

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
miraculix.dgemm_compressed.set_options(use_gpu=false)

println("Read bed file and frequencies")
genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)

println("Transpose matrix")
genotype_data_transposed = miraculix.dgemm_compressed.transpose_genotype_matrix(genotype_data, n_snps, n_indiv)

println("Decompress genotype data")
genotype_data_decompressed = miraculix.dgemm_compressed.decompress_genotype_data(genotype_data, n_indiv, n_snps)
genotype_data_t_decompressed = miraculix.dgemm_compressed.decompress_genotype_data(genotype_data_transposed, n_snps, n_indiv)

println("Initialize dgemm_compressed routine")
obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 10)

n_col = 10
B = randn(Float64, n_snps, n_col)

println("Calculate genotype matrix multiplication")
C = miraculix.dgemm_compressed.dgemm_compressed_main(false, obj_ref, B, n_snps, n_indiv)


println("Test results")
freq_test = mean(genotype_data_decompressed ./ 2.0, dims=1)
@test maximum(abs.(vec(freq_test) .- freq) ./ freq_test) < 1e-2

C_test = (genotype_data_decompressed .- 2.0 * freq_test) * B
@test maximum(abs.(C_test .- C)) < 1e-1

# miraculix.dgemm_compressed.free_compressed(obj_ref)