#  Authors 
#  Alexander Freudenberg, alexander.freudenberg@stads.de

#  Copyright (C) 2023 Alexander Freudenberg

#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

module solve
import ..LIBRARY_HANDLE, ..check_storage_object
using Libdl

"""
    sparse_init(V::Vector{Float64}, I::Vector{Int32}, J::Vector{Int32}, nnz::Int64, m::Int64, max_ncol::Int64)

Initializes storage object with required data on the GPU for solving an equation system defined by a sparse symmetric, positive-definite matrix A. 

# Arguments
- `V`: Vector of matrix values in COO format.
- `I`: Vector of row indices in COO format.
- `J`: Vector of column indices in COO format.
- `nnz`: The number of non-zero values in the matrix (length of V).
- `m`: The number of rows and columns in the matrix.
- `max_ncol`: The maximum number of columns in the right-hand side (RHS) matrix in equation systems.

# Returns 
A reference to the GPU storage object (`Ref{Ptr{Cvoid}}`).

# Exceptions
Throws an error if the initialization fails.

# Note
This function is an interface to the `sparse_solve_init` function in the `miraculix.so` library.
"""
function sparse_init(V::Vector{Float64}, I::Vector{Int64}, J::Vector{Int64}, nnz::Int64, m::Int64, max_ncol::Int64, is_lower::Bool)
    obj_ref = Ref{Ptr{Cvoid}}(C_NULL)
    if (length(V), length(I), length(J)) != (nnz, nnz, nnz)
        error("Unexpected length of vectors in COO format.")
    end
    status = zeros(Int32,1)

    init_sym = dlsym(LIBRARY_HANDLE[], :sparse2gpu)
    ccall(init_sym, Cvoid, (Ptr{Float64},Ptr{Int64},Ptr{Int64}, Int64, Int64, Int64, Int32, Ptr{Ptr{Cvoid}}, Ptr{Int32}), V, I, J, nnz, m, max_ncol, Int(is_lower), obj_ref, status)

    if status[1] != 0
        println("Status ", status)
        error("Routine not successful")
    end
    check_storage_object(obj_ref)
    
    return obj_ref
end # function

"""
    sparse_solve(obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, m::Int64)

Computes the solution to the equation system defined by the matrix in obj_ref and the RHS in B on the GPU.

# Arguments
- `obj_ref`: A reference to the GPU storage object.
- `B`: Right-hand side (RHS) matrix of size m x n.
- `m`: The number of rows in the RHS matrix.

# Returns 
The solution to the equation system.

# Exceptions
Throws an error if the computation fails.

# Note
This function is an interface to the `sparse_solve_compute` function in the `miraculix.so` library.
"""
function sparse_solve(obj_ref::Ref{Ptr{Cvoid}}, transA::Char, B::Matrix{Float64}, m::Int64)
    check_storage_object(obj_ref)

    if size(B,1) != m
        error("B must have $m rows to be compatible with the sparse matrix.")
    end
    ncol = Int(size(B, 2))
    X = zeros(Float64, (m, ncol))

    status = zeros(Int32,1)
    solve_sym = dlsym(LIBRARY_HANDLE[], :dcsrtrsv_solve_gpu)
    ccall(solve_sym, Cvoid, (Ptr{Cvoid}, Char, Ptr{Float64}, Int64, Ptr{Float64}, Ptr{Int32}), obj_ref[], transA, B, ncol, X, status)

    if status[1] != 0
        println("Status ", status)
        error("Routine not successful")
    end
    check_storage_object(obj_ref)

    return X
end # function

"""
    sparse_free(obj_ref::Ref{Ptr{Cvoid}})

Frees the GPU memory in obj_ref.

# Arguments
- `obj_ref`: A reference to the GPU storage object.

# Exceptions
Throws an error if the memory release fails.

# Note
This function is an interface to the `sparse_solve_destroy` function in the `miraculix.so` library.
"""
function sparse_free(obj_ref::Ref{Ptr{Cvoid}})
    check_storage_object(obj_ref)

    status = zeros(Int32,1)
    free_sym = dlsym(LIBRARY_HANDLE[], :free_sparse_gpu)
    ccall(free_sym, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Int32}), obj_ref, status)
    
    if status[1] != 0
        println("Status ", status)
        error("Routine not successful")
    end
end # function


"""
    dense_solve(M::Matrix{Float64}, B::Matrix{Float64}, calc_logdet::Bool = true, oversubscribe::Bool = true)

Computes the solution to the dense equation system defined by the symmetric, positive-definite M and B on the GPU.

# Arguments
- `M`: The input square matrix.
- `B`: Right-hand side (RHS) matrix with compatible number of rows with M.
- `calc_logdet`: Flag indicating if the log-determinant of M should be computed and returned. Default is `true`.
- `oversubscribe`: Flag indicating if CUDA unified memory addressing should be used, which allows the matrix M to be too large to fit in the GPU conventionally. Default is `false`.

# Returns 
The solution to the equation system and the log-determinant of M if `calc_logdet` is `true`.

# Exceptions
Throws an error if M is not a square matrix, B is not of correct dimension or the computation fails.

# Note
If `oversubscribe` is `true`, this might come with a performance penalty. This function is an interface to the `dense_solve` function in the `miraculix.so` library.
"""
function dense_solve(M::Matrix{Float64}, B::Matrix{Float64}; calc_logdet::Bool = true, oversubscribe::Bool = false)
    n = size(M,1)
    if n != size(M,2)
        error("M not a square matrix.")
    end
    if n != size(B,1)
        error("B not compatible with M.")
    end

    status = zeros(Int32,1)
    ncol = size(B,2)
    X = zeros(Float64, (n, ncol))
    if calc_logdet
        logdet = zeros(Float64,1)
    else
        logdet = Ptr{Cvoid}(C_NULL)
    end

    solve_sym = dlsym(LIBRARY_HANDLE[], :potrs_solve_gpu)
    ccall(solve_sym, Cvoid, (Ptr{Float64}, Int32, Ptr{Float64}, Int32, Ptr{Float64}, Ptr{Float64}, Int32, Ptr{Int32}), M, n, B, ncol, X, logdet, oversubscribe, status)

    if calc_logdet
        return X, logdet[1]
    else
        return X
    end
end # function


end # module