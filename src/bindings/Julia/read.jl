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


# 
#   TODOs: OneBit conversion, multithreading 
#



module read_plink

using Base;
using DelimitedFiles

"""
    read_bed(file::String, coding::String="TwoBit", snpmajor::Bool=true)::Matrix{Int32}

Reads a PLINK .bed binary file and stores the compressed genotype data in a single-precision integer matrix. 

# Arguments
- `file`: A string representing the path of the .bed file to read. It is expected that supplementary .bim and .fam files are present at the same location with the same base name as the .bed file.
- `coding`: A string specifying the storage format of the genotype data. Default value is "TwoBit". 
- `snpmajor`: A boolean value to determine if the stored matrix should be transposed. If true, the matrix is transposed. Default value is true.

The .bed file is a primary representation of genotype calls at biallelic variants, see https://www.cog-genomics.org/plink/1.9/formats#bed. It must be accompanied by .bim and .fam files. The first three bytes should be 0x6c, 0x1b, and 0x01 in that order. The rest of the file is a sequence of V blocks of N/4 (rounded up) bytes each, where V is the number of variants and N is the number of samples. Each block corresponds to a marker in the .bim file, with the low-order two bits of a block's first byte storing the first sample's genotype code, the next two bits storing the second sample's code, and so on.

# Returns
- A single-precision integer matrix (`Matrix{Int32}`) holding the genotype data from the .bed file in compressed format.

# Exceptions
- Throws an error if the .bed file or its supplementary .bim and .fam files do not exist or cannot be read.
- Throws an error if the .bed file does not follow the specified format.
"""
function read_bed(file::String, coding_twobit::Bool=false, snpmajor::Bool=true)

    if ~endswith(file,".bed")
        error("File not in .bed format")
    end
    if ~isfile(replace(file, ".bed" => ".fam")) | ~isfile(replace(file,".bed" => ".bim"))
        error("Missing supplementary file .fam or .bim")
    end
    size_in_bits = sizeof(UInt8) * 8;


    io = open(file, "r");
    start_char = zeros(UInt8, 3);
    unsafe_read(io, pointer(start_char), 3);
    if start_char != [0x6c, 0x1b, 0x01];
        error("Not a correct .bed file")
    end

    fam_contents = read(replace(file, ".bed" => ".fam"), String);
    bim_contents = read(replace(file, ".bed" => ".bim"), String);
    n_indiv = eachmatch(r"(\n)", fam_contents) |> collect |> length;
    n_snps = eachmatch(r"(\n)", bim_contents) |> collect |> length;
    n_bytes_per_col = Int(ceil(n_indiv/4));

    if snpmajor
        n_row = n_bytes_per_col;
        result = zeros(UInt8, (n_row, n_snps));

        # Read bed file - this throws an error if too small
        for i = 1:n_snps
            unsafe_read(io, pointer(result, (i-1) * n_row + 1), n_bytes_per_col);
        end
        # Assert end of file
        @assert eof(io) "Too large .bed file"
        close(io);

        if coding_twobit
            # Conversion to TwoBit format
            @inbounds for i = 1:n_row, j = 1:n_indiv
                index_str = bitstring(result[i,j]);
                new_entry = UInt8(0);
                @inbounds for substr_index = 1:2:size_in_bits
                    new_entry <<= 2;
                    substr = index_str[substr_index : (substr_index +1)];
                    # For documentation of values, see https://www.cog-genomics.org/plink/1.9/formats#bed
                    if substr == "10"
                        new_entry |= 1
                    elseif substr == "11"
                        new_entry |= 2
                    elseif substr == "01"
                        error("Missings in .bed not supported")
                    end
                end
                result[i,j] = new_entry;
            end            
        end # coding
    else
        error("Not implemented yet")        
    end # snpmajor
    
    return result, n_snps, n_indiv
end #function


"""
    read_freq(file::String)::DataFrame

Reads a file with the suffix `.freq` and extracts allele frequencies.

# Arguments
- `file`: A string specifying the path of the .freq file to read.

The .freq file is expected to contain allele frequency data with the first column indicating the SNP ID and the second column containing the frequency. 
Each row represents an allele with a specific frequency in the population.

# Returns
- A Vector where each row corresponds to an allele and its frequency in the population.

# Exceptions
- Throws an error if the .freq file does not exist or cannot be read.
- Throws an error if the .freq file is not in the expected format.
"""
function read_freq(file::String)
    freq = readdlm(file)[:,2]
    return freq
end 

end #module