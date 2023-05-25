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

module sparse_solve
import ..LIBRARY_HANDLE, check_storage_object
using Libdl

function init(V::Vector{Float64}, I::Vector{Int32}, J::Vector{Int32}, nnz::Int64, m::Int64, max_ncol::Int64)
    obj_ref = Ref{Ptr{Cvoid}}(C_NULL)
    if (length(V), length(I), length(J)) != (nnz, nnz, nzz)
        error("Unexpected length of vectors in COO format.")
    end
    status = Int32(0)

    init_sym = dlsym(LIBRARY_HANDLE[], :sparse2gpu)
    ccall(init_sym, Cvoid, (Ptr{Float64},Ptr{Int32},Ptr{Int32}, Int64, Int64, Int64, Ptr{Ptr{Cvoid}}, Ptr{Int32}), V, I, J, nnz, m, max_ncol, obj_ref, pointer(status))

    if status != 0
        error("Routine not successful")
    end
    check_storage_object(obj_ref)
    
    return obj_ref
end

function solve(obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, m::Int64)
    check_storage_object(obj_ref)

    if size(B,1) != m
        error("B must have $m rows to be compatible with the sparse matrix.")
    end
    ncols = Int32(size(B, 2))
    X = zeros(Float64, (m, ncol))

    status = Int32(0)
    solve_sym = dlsym(LIBRARY_HANDLE[], :dcsrtrsv_solve_gpu)
    ccall(solve_sym, Cvoid, (Ptr{Cvoid},Ptr{Float64}, Int32, Ptr{Float64}, Ptr{Int32}), obj_ref[], B, ncols, X, pointer(status))

    if (obj_ref[] == C_NULL) || (status != 0)
        error("Solve routine not successful")
    end

    return X
end

function free(obj_ref::Ref{Ptr{Cvoid}})
    check_storage_object(obj_ref)

    status = Int32(0)
    free_sym = dlsym(LIBRARY_HANDLE[], :free_sparse_gpu)
    ccall(free_sym, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Int32}), obj_ref, pointer(status))
    
    if status != 0
        error("Free routine not successful")
    end
end

end # module