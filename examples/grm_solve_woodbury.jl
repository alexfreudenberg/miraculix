
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


MODULE_PATH = "../src/bindings/Julia/dgemm_compressed.jl"
LIBRARY_PATH = "./src/miraculix/miraculix.so"

include(MODULE_PATH)

dgemm_compressed.set_library_path(LIBRARY_PATH)
dgemm_compressed.load_shared_library()
dgemm_compressed.set_options()

