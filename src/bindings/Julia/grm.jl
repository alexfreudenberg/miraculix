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


module grm

import ..LIBRARY_HANDLE, ..check_storage_object, ..check_dimensions
using Libdl
using LinearAlgebra

"""
    grm(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int)

TBW
"""
function compute(plink::Matrix{UInt8}, snps::Int, indiv::Int; do_center::Bool, is_transposed::Bool, scaling_method::Int, is_plink_format::Bool = false, allele_freq::Vector{Float64} = Vector{Float64}())
    if is_transposed
        nrow, ncol = snps, indiv
    else 
        nrow, ncol = indiv, snps
    end

    check_dimensions(plink, ncol, nrow)

    compute_sym = dlsym(LIBRARY_HANDLE[], :crossprod_mmagpu)

    M = zeros(Float64, (ncol, ncol))

    ccall(compute_sym,  Cint,  (Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint), plink, Int32(nrow), Int32(ncol), M, is_plink_format)
    
    @assert issymmetric(M) "Result not symmetric" 

    if do_center
        col_sum = sum(M, dims = 1)
        M .-= col_sum ./ indiv
        M .-= transpose(col_sum) ./ indiv
        M .+= sum(col_sum) ./ indiv^2
    end
    if scaling_method != 0
       # @assert isequal(length(allele_freq), snps) "Allele frequencies need to be supplied for scaling" 
        if scaling_method == 1
            c = 2 * sum(allele_freq .* (1 - allele_freq))
            M ./=c
        elseif scaling_method == 2
            allele_freq = reshape(allele_freq, (snps,1))
            M ./= allele_freq
            M ./= transpose(allele_freq)
        else 
            error("Wrong scaling method.")
        end
    end

    return M
end


end # module 