
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
DATA_FILE = "./data/small.bed"
FREQ_FILE = "./data/small.freq"

include(MODULE_PATH)

miraculix.dgemm_compressed.set_library_path(LIBRARY_PATH)
miraculix.dgemm_compressed.load_shared_library()
miraculix.dgemm_compressed.set_options()

genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)



# obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 10)
# error("success")

init_sym = dlsym(miraculix.dgemm_compressed.LIBRARY_HANDLE[], :plink2compressed)
using LinearAlgebra

# Define your inputs
plink =  rand(Char, 100, 100)  # replace with your actual data
plink_transposed = plink # transpose(plink)
snps = size(plink, 1) 
indiv = size(plink, 2) * 4
f = rand(snps)  # replace with your actual data
max_n = 10  # replace with your actual data
compressed = Ref{Ptr{Cvoid}}(C_NULL)
ccall(init_sym, Cvoid, 
      (Ptr{Char}, Ptr{Char}, Cint, Cint, Ptr{Float64}, Cint, Ptr{Ptr{Cvoid}}), 
      plink, plink_transposed, snps, indiv, f, max_n, compressed)