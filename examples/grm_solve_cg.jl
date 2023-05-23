
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
# - Throws an error if the path is wrong.
# - Converges poorly if the genomic relationship matrix G is not invertible (e.g., snps > indiv).
# - Throws an error if the PLINK bed file for G cannot be read or is incorrectly formatted.

using Statistics
using Libdl
using LinearAlgebra

ROOT_DIR = string(@__DIR__) * "/../"
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_FILE = ROOT_DIR * "/data/xsmall.bed"
FREQ_FILE = ROOT_DIR * "/data/xsmall.freq"

include(MODULE_PATH)

"""
    GRM_vec(obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, snps::Int, indiv::Int)

Multiplies the genomic relationship matrix (GRM) by a vector.

# Arguments
- `obj_ref`: A Ref{Ptr{Cvoid}} that holds a pointer to the storage object containing SNP data.
- `B`: A one-column Matrix of Float64 representing the vector to be multiplied.
- `snps`: An integer representing the number of SNPs.
- `indiv`: An integer representing the number of individuals.

This function first multiplies the transposed (centered) SNP matrix by the vector `B`, and then the result with the untransposed (centered) SNP matrix. This iterative process allows for efficient calculation of the GRM-vector multiplication.

# Returns
- The result of the GRM-vector multiplication.

# Exceptions
- Throws an error if the inputs are not in the expected format or if the multiplication operation fails.
"""
function GRM_vec(obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, snps::Int, indiv::Int)

    if size(B) != (indiv,1)
        error("Vector B must be of length $indiv")
    end
    
    Zv = miraculix.dgemm_compressed.dgemm_compressed_main(true, obj_ref, B, n_snps, n_indiv)
    Gv = miraculix.dgemm_compressed.dgemm_compressed_main(false, obj_ref, Zv, n_snps, n_indiv)

    return Gv
end

# Set library path and load library
miraculix.dgemm_compressed.set_library_path(LIBRARY_PATH)
miraculix.dgemm_compressed.load_shared_library()
# Set miraculix options
miraculix.dgemm_compressed.set_options(use_gpu=!false, verbose=0)

# Read PLINK and allele frequency files
genotype_data, n_snps, n_indiv = miraculix.read_plink.read_bed(DATA_FILE)
freq = miraculix.read_plink.read_freq(FREQ_FILE)

# Initialize storage object
obj_ref = miraculix.dgemm_compressed.init_compressed(genotype_data, n_snps, n_indiv, freq, 1)

# Set hyper parameters
max_iter = 1_000 # Maximum number of iterations
print_iter = 1e2 # Frequency of convergence information
conv_crit = 1e-2 # Maximum norm of residual

# Simulate data for toy example to test for convergence of algorithm
A = diagm(randn(Float64,n_indiv)./10 .+10) 
A_vec(obj_ref, x, snps, indiv) = A * x

# Start conjugate gradient 
begin 
    b = randn(Float64, n_indiv, 1) # RHS of equation system
    x = randn(Float64, n_indiv, 1) # Start vector x0

    Gx = GRM_vec(obj_ref, x, n_snps, n_indiv) # Calculate G * x0
    r = b - Gx # Initial residual
    p = r 
    for i in 1:max_iter
        global x, r, p, Gp, alpha, beta, norm_old
        norm_old = norm(r)
        # Exit if convergence criterion is reached
        if norm_old < conv_cirt
            break
        end
        
        Gp = GRM_vec(obj_ref, p, n_snps, n_indiv)
        alpha = norm_old^2 / (transpose(p) * Gp)[1]
        (i % print_iter == 0) && println(alpha, " ", norm_old)

        x += alpha * p
        r -= alpha * Gp
        beta = norm(r)^2/norm_old^2
        p = r + beta * p
    end
    println("Final residual norm: ", norm(r))
end