
# This script is designed to load the dgemm_compressed.jl module defined in src/bindings/Julia, and use 
# the genotpye matrix multiplication functionality within that module to compute the solution 
# to the equation G^-1 x, where x is a vector and G is the genomic relationship matrix.

# # Usage
# To use this script, make sure that the Julia file containing the module is in the 
# same directory as this script or its path is included in Julia's LOAD_PATH. Also, 
# ensure that the data for the genomic relationship matrix G and vector x are available 
# and correctly formatted.

# The script assumes that the module contains a function for matrix multiplication, 
# and that this function can handle the inversion of the genomic relationship matrix G.

# # Output
# The script will output the solution to the equation G^-1 x as a vector. 

# # Exceptions
# - Throws an error if the required module or its matrix multiplication function 
#   cannot be found or loaded.
# - Throws an error if the genomic relationship matrix G is not invertible.
# - Throws an error if the data for G or x cannot be read or is incorrectly formatted.

using Libdl

MODULE_PATH = "../src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = "./src/miraculix/miraculix.so"
DATA_FILE = "./data/xsmall.bed"
FREQ_FILE = "./data/xsmall.freq"

include(MODULE_PATH)

miraculix.dgemm_compressed.set_library_path(LIBRARY_PATH)
miraculix.dgemm_compressed.load_shared_library()
miraculix.dgemm_compressed.set_options(use_gpu=false)

genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)


obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 10)

n_col = 10
B = randn(Float64, n_snps, n_col)

C = miraculix.dgemm_compressed.dgemm_compressed_main(false, obj_ref, B, n_indiv, n_snps)

genotype_data_decompressed = zeros(Float64, n_indiv, n_snps);

@inbounds for index in 1:length(genotype_data)
    entry = genotype_data[index]
    offset_decompressed = (index-1) * 4 + 1
    @inbounds for i in 0:3
        genotype_data_decompressed[offset_decompressed + i] = max(0, Float64((entry >> (2*i)) & 0x03)-1)
    end    
end

freq_test = mean(genotype_data_decompressed ./ 2.0, dims=1)
max_deviation = maximum(abs.(vec(freq_test) .- freq))

C_test = (genotype_data_decompressed .- 2.0 * freq_test) * B
max_deviation = maximum(abs.(C_test .- C))
