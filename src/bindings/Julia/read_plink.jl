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



module read_plink

using Base;
using DelimitedFiles;
using StaticArrays;
using LoopVectorization;

# This lookup table for converting PLINK binary format to 2bit format has been created with the help of the create_conversion_table function below -- see its documentation for details
const CONVERSION_TABLE_UINT8 = Vector{UInt8}(
    [0, 255, 1, 2, 255, 255, 255, 255, 4, 255, 5, 6, 8, 255, 9, 10, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 16, 255, 17, 18, 255, 255, 255, 255, 20, 255, 21, 22, 24, 255, 25, 26, 32, 255, 33, 34, 255, 255, 255, 255, 36, 255, 37, 38, 40, 255, 41, 42, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 64, 255, 65, 66, 255, 255, 255, 255, 68, 255, 69, 70, 72, 255, 73, 74, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 80, 255, 81, 82, 255, 255, 255, 255, 84, 255, 85, 86, 88, 255, 89, 90, 96, 255, 97, 98, 255, 255, 255, 255, 100, 255, 101, 102, 104, 255, 105, 106, 128, 255, 129, 130, 255, 255, 255, 255, 132, 255, 133, 134, 136, 255, 137, 138, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 144, 255, 145, 146, 255, 255, 255, 255, 148, 255, 149, 150, 152, 255, 153, 154, 160, 255, 161, 162, 255, 255, 255, 255, 164, 255, 165, 166, 168, 255, 169, 170]
)
const CONVERSION_TABLE_UINT8_STATIC = SVector{typemax(UInt8)+1,UInt8}(CONVERSION_TABLE_UINT8)

"""
    convert_plink2twobit(entry::UInt8)::UInt8

Converts a single entry from the PLINK binary format (.bed) to a compressed 2-bit binary format using a predefined conversion table. The conversion table CONVERSION_TABLE_UINT8 should be previously defined.
"""
@inline function convert_plink2twobit(entry::UInt8)::UInt8   
    @inbounds result::UInt8 = CONVERSION_TABLE_UINT8_STATIC[entry + 1];
    return result
end #function
function convert_plink2twobit(entry)   
    @inbounds result = CONVERSION_TABLE_UINT8[entry + 1];
    return result
end #function

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
function read_bed(file::String; coding_twobit::Bool=false, calc_freq::Bool=false)

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

    n_row = n_bytes_per_col;
    result = zeros(UInt8, (n_row, n_snps));

    # Read bed file - this throws an error if too small
    for i = 1:n_snps
        unsafe_read(io, pointer(result, (i-1) * n_row + 1), n_bytes_per_col);
    end
    # Assert end of file
    @assert eof(io) "Too large .bed file"
    close(io);

    # Calculate allele frequencies
    # In PLINK binary format without missings, the allele frequency of a SNP corresponds to the number of set bits in this SNP 
    if calc_freq
        Popcounts = vmapt(count_ones, result)
        freq = sum(Popcounts, dims = 1)/(2 * n_indiv) |> vec
    end
    # Convert to 2bit format if requested
    if coding_twobit
        result = vmapt(convert_plink2twobit,result)
        # Test for missing values in original bed file
        @assert all(result .!= typemax(UInt8)) "No missings in PLINK file permitted."
    end # coding

    if calc_freq
        return result, freq, n_snps, n_indiv
    else
        return result, n_snps, n_indiv
    end
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
end #function 

"""
    print_conversion_table(::Type{T}, NaN_code = typemax(T), verbose = true)

Prints a lookup table for converting the PLINK binary format (.bed) to a compressed 2-bit binary format where:

- 0 represents homozygous for the major allele,
- 1 represents heterozygous,
- 2 represents homozygous for the minor allele.
Missing values are not supported and are coded as NaN_code.

# Arguments
- T: The type of the conversion table. Only UInt8 has been thoroughly tested.
- NaN_code: The code to use for missing values. Defaults to the maximum possible value of the type T.
- verbose: If true, then the conversion for every value from 0:typemax(T) is printed explicitly.

#Returns
- The conversion table as a vector.

# Exceptions
An error will be thrown if a non-supported type T is supplied.

# Examples
print_conversion_table(UInt8, NaN_code = 255, verbose = true)
"""
function create_conversion_table(NaN_code = typemax(UInt8), verbose = true)
    max_val = typemax(UInt8)
    conversion_table = zeros(UInt8, max_val+1)
    
    substr = Vector{String}(["00"])
    size_in_bits = sizeof(UInt8) * 8;

    # Conversion to TwoBit format
    for i in 0:max_val
        index_str = bitstring(UInt8(i));
        new_entry = UInt8(0);
         @inbounds for substr_index = 1:2:size_in_bits
            new_entry <<= 2;
            substr[1] = view(index_str, substr_index : (substr_index + 1));
            # For documentation of values, see https://www.cog-genomics.org/plink/1.9/formats#bed
            if substr[1] == "10"
                new_entry |= 1
            elseif substr[1] == "11"
                new_entry |= 2
            elseif substr[1] == "01"
                new_entry = UInt8(NaN_code)
                break
            end
        end
        verbose && println(index_str, " -> ",new_entry)
        conversion_table[i + 1] = new_entry;
    end
    str_conversion_table = join(string.(conversion_table), ", ")
    
    verbose && println("[", str_conversion_table, "]")
    return conversion_table
end #function

end #module