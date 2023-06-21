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


module compressed_operations

import ..LIBRARY_HANDLE, ..check_storage_object, ..check_dimensions
export transpose_genotype_matrix


"""
    transpose_genotype_matrix(plink::Matrix{UInt8}, snps::Int, indiv::Int)

Transposes the SNP-wise genotype matrix and returns it in a UInt8 storage format.

# Arguments
- `plink`: A Matrix of UInt8 representing the PLINK .bed matrix, which is SNP-major.
- `snps`: An integer representing the number of SNPs.
- `indiv`: An integer representing the number of individuals.

This function performs a transpose operation on the PLINK .bed matrix, converting it from a SNP-major format to a sample-major format, allowing for different types of computational operations.

# Returns
- A Matrix of UInt8, which is the transposed genotype matrix in a sample-major format.

# Exceptions
- Throws an error if the PLINK .bed matrix is not in the expected format, or if the transpose operation fails.
"""
function transpose_genotype_matrix(plink::Matrix{UInt8}, snps::Int, indiv::Int)
    check_dimensions(plink, snps, indiv)
    plink_transposed = zeros(UInt8, Int(ceil(snps/4)), indiv)

    for (index, entry) in pairs(IndexCartesian(), plink)
        id_indiv = index[1]
        id_snp = index[2]
        @inbounds for i in 0:3
            new_col = Int((id_indiv-1) * 4 + i + 1)
            new_row = Int(ceil(id_snp/4))
            offset = (id_snp-1) % 4
            plink_transposed[new_row, new_col] |= ((entry >> (2 * i)) & 0x03) << (2*offset)
        end
    end
    return plink_transposed
end



# This function decompresses genotype data in PLINK format for testing purposes -- throws an error if it finds a missing value
function decompress_plink_format(plink::Matrix{UInt8}, indiv::Int, snps::Int)
    decompressed = zeros(Float64, Int(ceil(indiv/4) * 4), snps);

    for (index,entry) in pairs(IndexLinear(),plink)
        offset_decompressed = (index-1) * 4 + 1
        @inbounds for i in 0:3
            # Convert packed SNP data to Float
            genotype_float = Float64((entry >> (2*i)) & 0x03)
            # Check if there is a missing value which is coded as 1 in PLINK format
            (genotype_float == 1.0) && error("Missing in genotype data")
            # Convert PLINK format to 0,1,2 format
            decompressed[offset_decompressed + i] = max(0, genotype_float -1)
        end    
    end
    decompressed = decompressed[1:indiv,:]
    return decompressed
end
function decompress_2bit_format(genotype_data::Matrix{UInt8}, indiv::Int, snps::Int)
    decompressed = zeros(Float64, Int(ceil(indiv/4) * 4), snps);

    for (index,entry) in pairs(IndexLinear(),genotype_data)
        offset_decompressed = (index-1) * 4 + 1
        @inbounds for i in 0:3
            # Convert packed SNP data to Float
            decompressed[offset_decompressed + i]= Float64((entry >> (2*i)) & 0x03)
        end    
    end
    decompressed = decompressed[1:indiv,:]
    return decompressed
end

end # module 