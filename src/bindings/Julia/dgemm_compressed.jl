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


module dgemm_compressed

using Base
using Libdl

const LIBRARY_PATH = Ref{String}()  # Reference to store the library path
const LIBRARY_HANDLE = Ref{Ptr{Cvoid}}() # Reference to store the library handle returned by Libdl


"""
    set_library_path(path::AbstractString)

Sets the global variable `LIBRARY_PATH` to the path of the shared library "miraculix.so".

# Arguments
- `path`: An AbstractString representing the full path to the shared library "miraculix.so". 

This function modifies the global variable `LIBRARY_PATH` in the module, which is then used 
by other functions to load and interact with the shared library.

# Exceptions
- Throws an error if the specified path is invalid or does not point to a "miraculix.so" library.
"""
function set_library_path(path::AbstractString)
    if !isfile(path)
        error("The specified library path does not exist: $path")
    end
    LIBRARY_PATH[] = path
end

"""
    load_shared_library()

Loads the shared library specified by the global variable `LIBRARY_PATH` using the `Libdl` module.

This function assumes that `LIBRARY_PATH` has been set to a valid path to the shared library 
"miraculix.so". If successful, the loaded library will be ready for use in the module.

# Exceptions
- Throws an error if the `LIBRARY_PATH` is not set or points to an invalid or inaccessible library.
"""
function load_shared_library()
    if isempty(LIBRARY_PATH[])
        error("Library path not set. Please call set_library_path with the path to the shared library.")
    end
    
    libpath = LIBRARY_PATH[]
    
    if !isfile(libpath)
        error("The specified library path does not exist: $libpath")
    end
    
    # Load the shared library
    LIBRARY_HANDLE[] = dlopen(libpath)
end

"""
    close_shared_library()

Closes the connection to the shared library "miraculix.so" using the `Libdl` module.

This function assumes that the shared library is currently loaded and accessible. 
It releases any resources associated with the library.

# Exceptions
- Throws an error if the shared library is not currently loaded or cannot be properly closed.
"""
function close_shared_library()
    if !isnull(LIBRARY_HANDLE[])
        dlclose(LIBRARY_HANDLE[])
        LIBRARY_HANDLE[] = C_NULL
    end
end

function check_dimensions(plink::Matrix{UInt8}, snps::Int, indiv::Int)
    if (size(plink, 1) != ceil(indiv/4)) || (size(plink, 2) != snps)
        error("Matrix has wrong dimensions: Expected $snps SNPs and $indiv samples")
    end
end

# This function decompresses genotype data in PLINK format 
function decompress_genotype_data(plink::Matrix{UInt8}, indiv::Int, snps::Int)
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

"""
    set_options(;use_gpu::Bool = false, cores::Int = 0, not_center::Bool = false, variant::Int = 0, verbose::Int = 1)

Sets computation options in the "miraculix.so" shared library for operations on SNP data of genotypes.

# Arguments
- `use_gpu`: A boolean indicating whether to use GPU acceleration for computations. If true and GPUs are available, computations are offloaded to the GPUs. Default is false.
- `cores`: An integer specifying the number of CPU cores to use for computations. If set to 0, all available cores are used. Default is 0.
- `not_center`: A boolean indicating whether to center the genotype matrix. If true, centering is turned off. Default is false.
- `variant`: An integer selecting between different implementations of the computation. The exact options depend on the version of the "miraculix.so" library. Default is 0.
- `verbose`: An integer setting the level of print details during the computations. A higher number results in more detailed output. Default is 1.

This function configures the shared library to perform efficient computations on compressed SNP data stored in PLINK .bed format.

# Exceptions
- Throws an error if the options could not be set correctly in the "miraculix.so" library.
"""
function set_options(;use_gpu::Bool = false, cores::Int = 0, not_center::Bool = false, variant::Int = 0, verbose::Int = 1)
    floatLoop = Int32(0)
    meanSubstract = Int32(0)
    ignore_missings = Int32(1)
    normalize = Int32(0)
    use_miraculix_freq = Int32(0)

    options_sym = dlsym(LIBRARY_HANDLE[], :setOptions_compressed)
    ccall(options_sym, Cvoid, 
        (Int32, Int32, Int32, Int32, Int32, Int32, Int32, Int32, Int32, Int32 ),
        Int32(use_gpu), Int32(cores), floatLoop, meanSubstract, ignore_missings, Int32(not_center), normalize, use_miraculix_freq, Int32(variant), Int32(verbose))
end

"""
    init_compressed(plink::Matrix{UInt8}, snps::Int, indiv::Int, freq::Vector{Float64}, max_ncol::Int)

Preprocesses a PLINK .bed SNP matrix for efficient computations using the "miraculix.so" shared library.

# Arguments
- `plink`: A Matrix{UInt8} representing the SNP matrix stored in PLINK .bed format.
- `snps`: An integer specifying the number of SNPs.
- `indiv`: An integer specifying the number of individuals.
- `freq`: A Vector{Float64} representing the allele frequencies.
- `max_ncol`: An integer used by the GPU (if enabled) to preallocate memory for the real-valued matrix. This parameter indicates the maximum number of columns that the real-valued matrix is allowed to have.

If the GPU usage is enabled via the `set_options` function, the SNP matrix and its transposed are copied to the GPU memory. If the GPU usage is disabled, the SNP matrix is converted to the 5codes format for optimized CPU usage.

# Returns 
- A Ref{Ptr{Cvoid}}, which stores a pointer to a storage object that holds the supplied data.

# Exceptions
- Throws an error if the PLINK .bed matrix or other inputs are not in the expected format.
"""
function init_compressed(plink::Matrix{UInt8}, snps::Int, indiv::Int, freq::Vector{Float64}, max_ncol::Int)
    obj_ref = Ref{Ptr{Cvoid}}(C_NULL)
    check_dimensions(plink, snps, indiv)
    plink_transposed = transpose_genotype_matrix(plink, snps, indiv)

    init_sym = dlsym(LIBRARY_HANDLE[], :plink2compressed)
    ccall(init_sym,  Cvoid,  (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint, Ptr{Ptr{Cvoid}}), plink, plink_transposed, Int32(snps), Int32(indiv), freq, Int32(max_ncol), obj_ref)

    return obj_ref
end


"""
    dgemm_compressed_main(transpose::Bool, obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, snps::Int, indiv::Int)

Performs matrix multiplication on the SNP data.

# Arguments
- `transpose`: A boolean indicating if the transpose of the genotype matrix is multiplied.
- `obj_ref`: A Ref{Ptr{Cvoid}} that holds a pointer to the storage object.
- `B`: A Matrix of Float64 representing the real-valued matrix.
- `snps`: An integer representing the number of SNPs.
- `indiv`: An integer representing the number of individuals.

This function performs matrix multiplication on the SNP data, either on the original or the transposed genotype matrix depending on the `transpose` parameter.

# Returns 
- A matrix C which stores the result of the genotype matrix multiplication

# Exceptions
- Throws an error if the matrix multiplication fails or the inputs are not in the expected format.
"""
function dgemm_compressed_main(transpose::Bool, obj_ref::Ref{Ptr{Cvoid}}, B::Matrix{Float64}, snps::Int, indiv::Int)
    trans = transpose ? "T" : "N"
    n_row = size(B, 1)
    n_col = size(B, 2)

    if n_row != (transpose ? indiv : snps)
        error("Matrix B is not compatible with genotype matrix of $snps SNPs and $indiv individuals when operation = $trans ")
    end

    C = zeros(Float64, transpose ? snps : indiv, n_col)

    dgemm_compressed_sym = dlsym(LIBRARY_HANDLE[], :dgemm_compressed)
    ccall(dgemm_compressed_sym, Cvoid, 
        (Ptr{Char}, Ptr{Cvoid}, Int32, Ptr{Float64}, Int32, Ptr{Float64}, Int32),
        pointer(trans), obj_ref[], n_col, B, n_row, C, size(C,1))
    
    return C
end

"""
    free_compressed(obj_ref::Ref{Ptr{Cvoid}})

Frees the memory allocated for the storage object.

# Arguments
- `obj_ref`: A Ref{Ptr{Cvoid}} that holds a pointer to the storage object.

This function releases the memory allocated to the storage object that holds the SNP data.

# Exceptions
- Throws an error if the memory cannot be properly freed.
"""
function free_compressed(obj_ref::Ref{Ptr{Cvoid}})
    free_sym = dlsym(LIBRARY_HANDLE[], :free_compressed)
    ccall(free_sym, Cvoid, (Ptr{Ptr{Cvoid}},), obj_ref)
end

end #module 