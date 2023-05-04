using LinearAlgebra
using SparseArrays
using Base
using DelimitedFiles
using BenchmarkTools
using DataFrames
using ProgressBars
using CSV
using Libdl
using CUDA
using Random
CUDA.set_runtime_version!("local")

#TODO: Disentagle this file

sparse_solve = function (V::Vector{Float64}, I_vec::Vector{Int}, J_vec::Vector{Int}, sp_nnz::Int, m::Int, B::Matrix{Float64})
	#
	# Provides an interface to a sparse triangular solve functionality in CUDA.
	#
	# This function takes a sparse matrix in COO (Coordinate List) format and a right-hand side matrix, and solves the linear system A * X = B, where A is the sparse matrix and X is the unknown matrix.
	#
	# Arguments
	# I_vec::Vector{Int}: A vector containing the row indices of the non-zero elements in the sparse matrix.
	# J_vec::Vector{Int}: A vector containing the column indices of the non-zero elements in the sparse matrix.
	# V::Vector{Float64}: A vector containing the non-zero values of the sparse matrix.
	# m::Int: The number of rows and columns of the sparse matrix.
	# nnz::Int: The number of non-zero elements in the sparse matrix.
	# B::Matrix{Float64}: The right-hand side matrix in the linear system to be solved.
	# Returns
	# X::Matrix{Float64}: The solution matrix X for the linear system A * X = B.
	#



	# Open shared library and define references to its symbols
	REP_FOLDER = ENV["PWD"];
	lib = dlopen(REP_FOLDER * "/src/solve_gpu.so")
	sparse_init_sym = dlsym(lib, :sparse2gpu)
	sparse_solve_sym = dlsym(lib, :dcsrtrsv_solve)
	sparse_destroy_sym = dlsym(lib, :freegpu_sparse)

	# Auxiliary structures for C functions
	X = similar(B)
	ncol_B = size(B)[2]
	GPU_obj_ref = Ref{Ptr{Cvoid}}(C_NULL)
	status_ref = Ref{Cint}(0)


	# Call init function and check for errors
	ccall(sparse_init_sym, Cvoid, (Ptr{Cdouble}, Ptr{Int32}, Ptr{Int32}, Int64, Int64, Int32, Ptr{Ptr{Cvoid}}, Ptr{Int32}), V, Vector{Int32}(I_vec), Vector{Int32}(J_vec), Int64(sp_nnz), Int64(m), Int64(ncol_B), GPU_obj_ref, status_ref)
	if status_ref[] != 0
		println("Error in init function")
		dlclose(lib)
		error(1)
	end

	# Call solve function and check for errors
	@btime ccall(sparse_solve_sym, Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Int32, Ptr{Cdouble}, Ptr{Int32}), GPU_obj_ref[], B, ncol_B, X, status_ref)
	if status_ref[] != 0
		println("Error in solve function")
		dlclose(lib)
		error(1)
	end

	# Call destroy function and check for errors
	ccall(sparse_destroy_sym, Cvoid, (Ptr{Cvoid}, Ptr{Int32}), GPU_obj_ref[], status_ref)
	if status_ref[] != 0
		println("Error in destroy function")
		dlclose(lib)

	end

    CUDA.CUDA.device_reset!()

	return X
end


# I test the function on three matrices of varying sizes and structure - these are set-up here
println("Setting test matrices")
L_anc = sparse(df_IJV.I, df_IJV.J, df_IJV.V, maximum(df_IJV.J), maximum(df_IJV.J))
ngeno = sort(unique(df_IJV.I))
L_anc = L_anc[ngeno, ngeno]
L_small = sparse(deepcopy(L_anc[1:1_000,1:1_000]))
L_small_filled = deepcopy(L_small)
for i in 1:min(size(L_small_filled)[1], 40)
    for j in 1:(i-1)
        L_small_filled[i, j] = randn(1)[1]
    end
end

# Start the calculation and get the absolute error of the solution in the equation systems 
println("Start calculation")
for S in [L_small, L_small_filled, L_anc]

	# Convert S to COO format
    I_vec, J_vec, V = findnz(S);
	df_IJV_loc = DataFrame(I = I_vec, J = J_vec, V = V)
	df_IJV_loc = df_IJV_loc[randperm(nrow(df_IJV_loc)),:]
    m = maximum(J_vec);
    nnz = length(I_vec);
    println(m, " ", nnz, " ", length(S.nzval))

	# Simulate RHS
    ncol_B = 6
    B = randn((m, ncol_B))

	# Call solve function
    @time X = sparse_solve(df_IJV_loc.V, df_IJV_loc.I, df_IJV_loc.J, nnz, m, B)

    deviation = sum(abs.(S * (transpose(S) * X) - B))
    println("Deviation ", deviation)
end 