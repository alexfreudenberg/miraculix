using Base, Base.Libc.Libdl;
using CUDA;

module solve
function solve_gpu{T}(A::Matrix,b::Vector)

    x = zeros(size(b));
    ccall((:main, "./solve_gpu.so"),Cvoid,())

end #function

end