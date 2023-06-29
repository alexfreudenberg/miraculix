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


module miraculix
export LIBRARY_HANDLE, LIBRARY_PATH, check_storage_object, check_library_handle

using Base
using Libdl

const LIBRARY_PATH = Ref{String}()  # Reference to store the library path
const LIBRARY_HANDLE = Ref{Ptr{Cvoid}}(C_NULL) # Reference to store the library handle returned by Libdl

# Error handling for uninitialized pointers
function check_handle(handle::Ref{Ptr{Cvoid}})
    if handle[] == C_NULL
        error("Encountered uninitialized pointer.")
    end
end
check_storage_object(obj_ref::Ref{Ptr{Cvoid}}) = check_handle(obj_ref)
check_library_handle() = check_handle(LIBRARY_HANDLE)

# Check if genotype matrix has correct dimensions
function check_dimensions(plink::Matrix{UInt8}, snps::Int, indiv::Int)
    if (size(plink, 1) != ceil(indiv/4)) || (size(plink, 2) != snps)
        error("Matrix has wrong dimensions: Expected $snps SNPs and $indiv samples")
    end
end

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

include("compressed_operations.jl")
include("dgemm_compressed.jl")
include("read_plink.jl")
include("solve.jl")
include("crossproduct.jl")

end