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
import ..LIBRARY_HANDLE
using Libdl

function init(V::Vector{Float64}, I::Vector{Int32}, J::Vector{Int32}, nnz::Int64, m::Int64, max_ncol::Int64)
    obj_ref = Ref{Ptr{Cvoid}}(C_NULL)
    if (length(V), length(I), length(J)) != (nnz, nnz, nzz)
        error("Unexpected length of vectors in COO format.")
    end
    status = 0

    init_sym = dlsym(LIBRARY_HANDLE[], :sparse2gpu)
    ccall(init_sym, Cvoid, (Ptr{Float64},Ptr{Int32},Ptr{Int32}, Int64, Int64, Int64, Ptr{Ptr{Cvoid}}, Ptr{Int32}), V, I, J, nnz, m, max_ncol, obj_ref, pointer(status))
    if obj_ref[] == C_NULL
        error("Routine not successful")
    end
    
    return obj_ref
end

end # module