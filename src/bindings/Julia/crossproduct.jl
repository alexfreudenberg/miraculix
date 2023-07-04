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
    snp_crossprod(plink::Matrix{UInt8}, snps::Int, indiv::Int; is_snpmajor::Bool, is_plink_format::Bool = false)

Computes the crossproduct of the SNP matrix with itself using a C library interface.

# Arguments
- `plink::Matrix{UInt8}`: The SNP matrix in compressed 2bit format.
- `snps::Int`: The number of SNPs in the SNP matrix.
- `indiv::Int`: The number of individuals in the SNP matrix.
- `is_snpmajor::Bool`: Indicates if the SNP matrix is in SNP-major format.
- `is_plink_format::Bool = false`: Specifies if the SNP matrix is in PLINK binary format. Defaults to `false`.

# Returns
- The SNP matrix crossproduct as a Float64 matrix of dimensions snps times snps if is_snpmajor is true and indiv times indiv else.

# Errors
- Throws an error if the dimensions of the SNP matrix are incorrect or if the library handle for the miraculix library is not initialized.

"""
function snp_crossprod(plink::Matrix{UInt8}, snps::Int, indiv::Int; is_snpmajor::Bool, is_plink_format::Bool = false)
    if is_snpmajor
        nrow, ncol = indiv, snps
    else 
        nrow, ncol = snps, indiv
    end

    check_dimensions(plink, ncol, nrow)

    check_library_handle()
    compute_sym = dlsym(LIBRARY_HANDLE[], :snp_multiply_gpu)

    M = zeros(Float64, (ncol, ncol))

    wtime = @elapsed ccall(compute_sym,  Cint,  (Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint), plink, Int32(nrow), Int32(ncol), M, is_plink_format)
    @debug "Time for calculating crossproduct: $wtime s."

    @assert issymmetric(M) "Result not symmetric" 

    return M
end # function

"""
    grm(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int; is_plink_format::Bool = false, do_scale::Bool = true, allele_freq::Vector{Float64} = Vector{Float64}())

Calculates the genomic relationship matrix (GRM) using the SNP crossproduct function.

# Arguments
- `plink_transposed::Matrix{UInt8}`: The SNP matrix in individual-major format.
- `snps::Int`: The number of SNPs in the SNP matrix.
- `indiv::Int`: The number of individuals in the SNP matrix.
- `is_plink_format::Bool = false`: Specifies if the SNP matrix is in PLINK binary format. Defaults to `false`.
- `do_scale::Bool = true`: Indicates whether to the SNP matrix is scaled by the sum of allele variances. Defaults to `true`.
- `allele_freq::Vector{Float64} = Vector{Float64}()`: The vector of allele frequencies. Needs to be supplied if `do_scale = true`. Defaults to an empty vector.  

# Returns
- The genomic relationship matrix (GRM) as a Float64 matrix.

"""
function grm(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int; is_plink_format::Bool = false, do_scale::Bool = true, allele_freq::Vector{Float64} = Vector{Float64}())
    if do_scale && length(allele_freq) != snps
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

"""
    ld(plink::Matrix{UInt8}, snps::Int, indiv::Int; is_plink_format::Bool = false, allele_freq::Vector{Float64} = Vector{Float64}())

Calculates the Linkage Disequilibrium (LD) statistic R^2 based on allele counts using the SNP crossproduct function.

# Arguments
- `plink::Matrix{UInt8}`: The SNP matrix in compressed 2bit format.
- `snps::Int`: The number of SNPs in the SNP matrix.
- `indiv::Int`: The number of individuals in the SNP matrix.
- `is_plink_format::Bool = false`: Specifies if the SNP matrix is in PLINK binary format. Defaults to `false`.
- `allele_freq::Vector{Float64} = Vector{Float64}()`: The vector of allele frequencies. Defaults to an empty vector.

# Returns
- The R^2 matrix as a Float64 matrix.

"""
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