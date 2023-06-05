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

# 
# 
#

using Distances;
using SparseArrays;
using LinearAlgebra;
using Test;
using Random;
using Printf;

# =====================
# Global definitions
# =====================

ROOT_DIR = string(@__DIR__) * "/../.."
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"

tol = 1e-1;
Random.seed!(0);

# Remove commit message verbosity
ENV["PRINT_LEVEL"] = "1";

include(MODULE_PATH)

# =====================
# Auxiliary functions
# =====================

"""
    simulate_sparse_pd(n::Int, density = 0.01)

Generates a sparse, lower triangular square matrix of size `n` x `n` with specified `density`. 

# Arguments
- `n::Int`: The dimension of the square matrix to be generated.
- `density::Float64=0.01`: Roughly the density of the matrix. Defaults to 0.01 if not specified, representing 1% density.
- `bias::Float64=2.0`: Mean of the diagonal elements of the matrix.

# Returns
- A sparse, positive-definite square matrix of size `n` x `n` with approx. specified density.
"""
function simulate_sparse_triangular(n::Int, density::Float64 = 0.01, bias::Float64 = 2.0)
    diag_elements = max.(randn(n) .+ bias, 0.1);
    M_sp = spdiagm(diag_elements);
    indices = rand(1:n, (Int(2 * ceil(density * n)), 2));
    # Generate a sparse, strictly diagonally dominant matrix M_sp through sampling values such that each off-diagonal rowsum is smaller than the diagonal value 1.0
    for (i,j) in eachrow(indices)
        if i >= j
            continue
        end
        M_sp[i,j] = rand(Float64) * (0.9 - sum(abs.(M_sp[i,:])));
    end

    return M_sp
end

"""
    simulate_dense_pd(n::Int)

Simulates a positive definite matrix of dimension n.

Parameters:
- n::Int: The dimension of the square matrix.

Returns:
- M: Simulated positive definite matrix of size n x n.
"""
function simulate_dense_pd(n::Int)
    R = pairwise(Euclidean(),1:n);
    M = exp.(- R ./n);
    return M
end

# =====================
# Main
# =====================

println("Load library and set options")
miraculix.set_library_path(LIBRARY_PATH)
miraculix.load_shared_library()


println("Check if routine returns right results")
@testset "Consistency" begin
    for n in Vector{Int64}([1e2,5e2,5e3])
        for ncol in [1, 5, 20]
            for density in [0.05, 0.2, 0.9]
                @printf("n: %d, ncol: %d, density: %.2f\n", n, ncol, density)
                # Simulate LHS and RHS
                M_sp = simulate_sparse_triangular(n, density);
                M = simulate_dense_pd(n);                
                B = randn(Float64, (n, ncol)) .+ 5; # Add bias to avoid accidentally correct results
                # Convert M to COO format
                I, J, V = findnz(M_sp)

                # Initialize GPU storage object from COO 
                obj_ref = miraculix.solve.sparse_init(V, Vector{Int32}(I), Vector{Int32}(J), length(I), n, ncol, false)
                # Compute the solution to M_sp X_sp = B
                X_sp = miraculix.solve.sparse_solve(obj_ref, B, n)
                # Free GPU memory
                miraculix.solve.sparse_free(obj_ref)
                @test_throws "No valid storage object" miraculix.solve.sparse_free(obj_ref)
            
                # Compute the solution to M X = B
                X = miraculix.solve.dense_solve(M, B, calc_logdet = false, oversubscribe = false)
                
                # Calculate deviations
                D = abs.(M_sp * transpose(M_sp) * X_sp - B)
                @printf("Absolute error: %.1e, maximum error: %.1e\n", norm(D), maximum(D))

                @test norm(D)/norm(B) < tol
                @test norm(M * X - B)/norm(B) < tol
            end
        end
    end
end

println("Check if routine is resilient - uncaught memory allocations would cause this to fail")
@testset "Resilience" begin
    iter = 1e2
    n = Int(1e4)
    ncol = 12
    M_sp = simulate_sparse_pd(n, 0.05);
    M = simulate_dense_pd(n);
    B = randn(Float64, (n, ncol));
    # Convert M to COO format
    I, J, V = findnz(M_sp)
     
    # Initialize GPU storage object from COO 
    obj_ref = miraculix.solve.sparse_init(V, Vector{Int32}(I), Vector{Int32}(J), length(I), n, ncol)
    # Repeat solve call a lot of times
    for _ in 1:iter
        miraculix.solve.sparse_solve(obj_ref, B, n)
        miraculix.solve.dense_solve(M, B, calc_logdet = false, oversubscribe = false)
    end
    # Compute the solution to M_sp X_sp = B and M X = B
    X_sp = miraculix.solve.sparse_solve(obj_ref, B, n)
    X, logdet_own = miraculix.solve.dense_solve(M, B, calc_logdet = true, oversubscribe = false)

    # Free GPU memory
    miraculix.solve.sparse_free(obj_ref)
    @test norm(M_sp * X_sp - B)/norm(B) < tol
    @test norm(M * X - B)/norm(B) < tol
    @test abs(logdet(M) - logdet_own) < tol 

end

println("Check if oversubscription and logdet calculation works -- this test needs to be adjusted to actual device memory available")
@testset "Oversubscription" begin
    for n in Vector{Int64}([1e4,3e4,7e4])
        println("Testing size $n")
        ncol = 1
        # Simulate LHS and RHS
        M = simulate_dense_pd(n);
        B = randn(Float64, (n, ncol));
        
        # Compute the solution to M X = B
        X = miraculix.solve.dense_solve(M, B, calc_logdet = false, oversubscribe = true)

        println("Test correctness")
        @test norm(M * X - B)/norm(B) < tol
    end
end