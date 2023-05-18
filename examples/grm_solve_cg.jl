
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
using LinearAlgebra

ROOT_DIR = string(@__DIR__) * "/../"
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_FILE = ROOT_DIR * "/data/xsmall.bed"
FREQ_FILE = ROOT_DIR * "/data/xsmall.freq"

include(MODULE_PATH)

function GRM_vec(obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, snps::Int, indiv::Int)

    if length(B) != indiv
        error("Vector must be of length $indiv")
    end
    
    Zv = miraculix.dgemm_compressed.dgemm_compressed_main(true, obj_ref, B, n_snps, n_indiv)
    Gv = miraculix.dgemm_compressed.dgemm_compressed_main(false, obj_ref, Zv, n_snps, n_indiv)

    return Gv
end

miraculix.dgemm_compressed.set_library_path(LIBRARY_PATH)
miraculix.dgemm_compressed.load_shared_library()
miraculix.dgemm_compressed.set_options(use_gpu=false, verbose=0)

genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)


obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 10)

n_indiv = 10
max_iter = 100
A = diagm(randn(Float64,n_indiv).+10)
A_vec(obj_ref, x, snps, indiv) = A * x
begin 
    b = randn(Float64, n_indiv, 1)
    x = randn(Float64, n_indiv, 1)

    Gx = A_vec(obj_ref, x, n_snps, n_indiv)
    r = b - Gx
    p = r
    for i in 1:max_iter
        global x, r, p, Gp
        norm_old = norm(r)
        if norm_old < 1e-2
            break
        end
        
        Gp = A_vec(obj_ref, p, n_snps, n_indiv)
        alpha = norm_old / (transpose(p) * Gp)[1]
        println(alpha, " ", norm_old)

        x += alpha * p
        r -= alpha * Gp
        beta = norm(r)/norm_old
        p = r + beta * p
    end
end