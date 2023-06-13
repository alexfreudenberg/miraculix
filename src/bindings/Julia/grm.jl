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

"""
    grm(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int)

TBW
"""
function compute(plink_transposed::Matrix{UInt8}, snps::Int, indiv::Int, is_plink_format::Bool)
    check_dimensions(plink_transposed, indiv, snps)

    compute_sym = dlsym(LIBRARY_HANDLE[], :crossprod_mmagpu)

    M = zeros(Float64, (indiv, indiv))

    ccall(compute_sym,  Cint,  (Ptr{UInt8}, Cint, Cint, Ptr{Float64}, Cint), plink_transposed, Int32(snps), Int32(indiv), M, is_plink_format)
    # void crossprod_mmagpu(char *snp_matrix, int snps, int indiv, double *ans)
    
    return M
end


end # module 