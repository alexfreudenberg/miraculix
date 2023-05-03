#
#   This file implements a read function for plink binary files (.bed)
#   The data is subsequently converted to a miraculix TwoBit storage format
# 


# 
#   TODOs: OneBit conversion, multithreading 
#

using Base;
const TYPE = UInt128;

file = "/Users/alex/Downloads/plink_mac_20220402/plink.bed"

module read

function read_bed(::Type{T}, file::String, coding::String="TwoBit", snpmajor::Bool=true)::T where{T}

    if ~endswith(file,".bed")
        error("File not in .bed format")
    end
    if ~isfile(replace(file, ".bed" => ".fam")) | ~isfile(replace(file,".bed" => ".bim"))
        error("Missing supplementary file .fam or .bim")
    end
    const TYPESIZE = sizeof(T) * 8;


    io = open(file, "r");
    start_char = zeros(UInt8, 3);
    unsafe_read(io, pointer(start_char), 3);
    if start_char != [0x6c, 0x1b, 0x01];
        error("Not a correct .bed file")
    end

    fam_contents = read(replace(file, ".bed" => ".fam"), String);
    bim_contents = read(replace(file, ".bed" => ".bim"), String);
    n_indiv = eachmatch(r"(\n)", fam_contents) |> collect |> length;
    n_snps = eachmatch(r"(\n)", fam_contents) |> collect |> length;
    n_bytes_per_row = Int(ceil(n_indiv/4));

    if snpmajor
        n_row = Int(ceil(2 * n_snps/ TYPESIZE));
        result = zeros(T, (n_row, n_indiv));

        if coding == "TwoBit"
            # Read bed file - this throws an error if too small
            for i = 1:n_snps
                unsafe_read(io, pointer(result, (i-1) * n_row + 1), n_bytes_per_row);
            end
            # Assert end of file
            @assert eof(io) "Too large .bed file"
            close(io);

            # Conversion to TwoBit format
            @inbounds for i = 1:n_row, j = 1:n_indiv
                index_str = bitstring(result[i,j]);
                new_entry = T(0);
                @inbounds for substr_index = 1:2:TYPESIZE
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
            
        else
            error("Not implemented yet")
        end
        
    end


end #function


end #module