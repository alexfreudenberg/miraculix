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


module crossproduct

import ..LIBRARY_HANDLE, ..check_storage_object, ..check_dimensions, ..check_library_handle
using Libdl
using LinearAlgebra

"""
    crossproduct(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int)

TBW
"""
function snp_crossprod(plink::Matrix{UInt8}, snps::Int, indiv::Int; is_snpmajor::Bool, is_plink_format::Bool = false)
    if is_snpmajor
        nrow, ncol = indiv, snps
    else 
        nrow, ncol = snps, indiv
    end

    check_dimensions(plink, ncol, nrow)

    check_library_handle()
    compute_sym = dlsym(LIBRARY_HANDLE[], :crossprod_mmagpu)

    M = zeros(Float64, (ncol, ncol))

    wtime = @elapsed ccall(compute_sym,  Cint,  (Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint), plink, Int32(nrow), Int32(ncol), M, is_plink_format)
    @debug "Time for calculating crossproduct: $wtime s."

    @assert issymmetric(M) "Result not symmetric" 

 

    return M
end # function


function grm(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int; is_plink_format::Bool = false, do_scale::Bool = true, allele_freq::Vector{Float64} = Vector{Float64}())
    if length(allele_freq) != snps
        error("Allele frequencies need to be equal to length of SNPs $snps.")
    end

    M = snp_crossprod(plink_transposed, snps, indiv, is_snpmajor = false, is_plink_format = is_plink_format)

    # Scaling of centered genotype matrix
    col_sum = sum(M, dims = 1) |> vec
    one_vector = ones(Float64, (indiv,)) 

    wtime = @elapsed begin 
        BLAS.ger!(-1/indiv, col_sum, one_vector, M)
        BLAS.ger!(-1/indiv, one_vector, col_sum, M)
    end
    @debug "Time for rank-one updates of GRM: $wtime s."

    wtime = @elapsed M .+= sum(col_sum) ./ indiv^2
    @debug "Time for affine transform of GRM: $wtime s."
    
    wtime = @elapsed if do_scale
        c = 2 * sum(allele_freq .* (1.0 .- allele_freq))
        M ./=c
    end
    @debug "Time for scaling GRM: $wtime s."
    
    return M
end # function

function ld(plink::Matrix{UInt8}, snps::Int, indiv::Int; is_plink_format::Bool = false, allele_freq::Vector{Float64} = Vector{Float64}())
    if length(allele_freq) != snps
        error("Allele frequencies need to be equal to length of SNPs $snps.")
    end

    M = snp_crossprod(plink, snps, indiv, is_snpmajor = true, is_plink_format = is_plink_format)

    # Scaling of centered genotype matrix
    wtime = @elapsed begin
        BLAS.syr!('U', Float64(-4.0 * indiv), allele_freq, M)
        M = Matrix(Symmetric(M, :U))
    end
    @debug "Time for rank-one update of LD: $wtime s."

    # Calculate vector of standard deviations
    wtime = @elapsed begin
        sigma_vector = reshape(sqrt.(diag(M)), (snps,1))
        # Device each row and column by the vector of standard deviations
        M ./= sigma_vector
        M ./= transpose(sigma_vector)
    end
    @debug "Time for scaling LD by std devs: $wtime s."

    return M
end # function

end # module 