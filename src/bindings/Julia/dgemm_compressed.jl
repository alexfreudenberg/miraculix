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
const LIBRARY_HANDLE = Ref{Ptr{Cvoid}}()

function set_library_path(path::AbstractString)
    LIBRARY_PATH[] = path
end


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

function close_shared_library()
    if !isnull(LIBRARY_HANDLE[])
        dlclose(LIBRARY_HANDLE[])
        LIBRARY_HANDLE[] = C_NULL
    end
end

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

function init_compressed(plink::Matrix{UInt8}, snps::Int, indiv::Int, freq::Vector{Float64}, max_ncol::Int)
    obj_ref = Ref{Ptr{Cvoid}}(C_NULL)

    init_sym = dlsym(LIBRARY_HANDLE[], :plink2compressed)
    ccall(init_sym,  Cvoid,  (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint, Ptr{Ptr{Cvoid}}), plink, plink, Int32(snps), Int32(indiv), freq, Int32(max_ncol), obj_ref)

    return obj_ref
end


end #module 