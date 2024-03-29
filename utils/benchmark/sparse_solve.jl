
#See https://github.com/JuliaSparse/Pardiso.jl for the documentation
using Pardiso;
using DataFrames;
using CSV;
using SparseArrays;
using Random;

# =====================
# Global definitions
# =====================

ROOT_DIR = string(@__DIR__) * "/../.."
MODULE_PATH = ROOT_DIR * "/src/bindings/Julia/miraculix.jl"
LIBRARY_PATH = ROOT_DIR * "/src/miraculix/miraculix.so"
DATA_FILE = ROOT_DIR * "/data/fort.dat_chol"

tol = 1e-1;

# Remove commit message verbosity
ENV["PRINT_LEVEL"] = "1";

include(MODULE_PATH)

miraculix.set_library_path(LIBRARY_PATH)
miraculix.load_shared_library()

# Read datafile
A=CSV.read(DATA_FILE,delim=' ', header = ["I", "J", "V"], DataFrame)

n = maximum(A[:,2])
ncol = 12
I = A[:,1]
J = A[:,2]
V = A[:,3]
B = randn(Float64, (n, ncol)) .+ 5; # Add bias to avoid accidentally correct results

# miraculix initialize
obj_ref = miraculix.solve.sparse_init(V, Vector{Int64}(I), Vector{Int64}(J), length(I), n, ncol, false)
@time begin 
    # Compute the solution to M_sp^T Y_sp = B
    Y_sp = miraculix.solve.sparse_solve(obj_ref, 't', B, n)
    # Compute the solution to M_sp X_sp = Y_sp
    X_sp = miraculix.solve.sparse_solve(obj_ref, 'n', Y_sp, n)
end
# Free GPU memory
miraculix.solve.sparse_free(obj_ref)

# Time Pardiso
# ps = PardisoSolver()
# M = sparse(I, J, V)

# X = solve(ps, M, B)
