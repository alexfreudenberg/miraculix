
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

using Statistics
using Libdl

ROOT_DIR = string(@__DIR__) * "/../.."
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_FILE = ROOT_DIR * "/data/xsmall.bed"
FREQ_FILE = ROOT_DIR * "/data/xsmall.freq"

include(MODULE_PATH)

miraculix.dgemm_compressed.set_library_path(LIBRARY_PATH)
miraculix.dgemm_compressed.load_shared_library()
miraculix.dgemm_compressed.set_options(use_gpu=false)

genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)


obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 10)

n_col = 10
B = randn(Float64, n_snps, n_col)

C = miraculix.dgemm_compressed.dgemm_compressed_main(false, obj_ref, B, n_snps, n_indiv)

