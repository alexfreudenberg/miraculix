using LinearAlgebra
using SparseArrays
using Base
using DelimitedFiles
using BenchmarkTools
using DataFrames
using ProgressBars
using CSV
using Graphs
using Traceur
using TimerOutputs
using Libdl

ped_file = "/home/hk-project-test-miraculx/ma_afreuden/ICBF_data/ICBF/hpPedigree.txt"
fam_file = "/home/hk-project-test-miraculx/ma_afreuden/ICBF_data/ICBF/hpSNPcov001.fam"
inb_file = "/home/hk-project-test-miraculx/ma_afreuden/ICBF_data/ICBF/hpPedigree.txt.inb"


const to = TimerOutput()
function calc_ped_df(P_anc::DataFrame)
    P_non = filter(:geno => !, P_anc)
    P_non = rename(P_non, :inb_sire => :inb_0, :inb_dam => :inb_1, :new_index_sire => :index_0, :new_index_dam => :index_1)
    P_non = select(P_non, Not([:x1,:x2,:x3]))
    n_new_anc = sum(P_non.index_0 .!= 0) + sum(P_non.index_1 .!= 0)
    counter = 1
    P_sub = P_non
    while n_new_anc > 0
        anc_cols = names(P_sub)
        anc_cols = anc_cols[count_digits.(anc_cols) .== counter]
        anc_cols = [col for col in anc_cols if !occursin("inb",col)]
        @timeit to "Calc" n_new_anc = sum([sum(col.!= 0) for col in eachcol(P_sub[:,anc_cols])] )
        println(n_new_anc)
        (n_new_anc == 0 )&& break

        for col in ProgressBar(anc_cols)
            all(P_sub[:,col] .== 0) && continue
            number = match(r"\d+",col).match
            @timeit to "Subsetting" new_anc = filter(:new_index => in(Set(P_sub[:,col])), P_non)
            @timeit to "Selecting" new_anc = select(new_anc, :new_index => col, :index_0 => col * "0", :index_1 => col * "1", :inb_0 => "inb_" * number * "0", :inb_1 => "inb_" * number * "1" )
            @timeit to "Joining" global P_sub = leftjoin(P_sub, new_anc, on = col)

        end
        @timeit to "Coalescing" global P_sub = coalesce.(P_sub, 0)

        println(n_new_anc)
        println(counter)
        global counter += 1
        (counter > 6 ) && break
    end
end
# P_anc = CSV.read("./data/pedigree_ancestors.csv", DataFrame)
# P_non = filter(:geno => !, P_anc)
#CSV.write("./data/pedigree_df.csv", calc_ped_df(P_non))
#exit(0)

#P_anc = CSV.read("./data/pedigree_df.csv",DataFrame)
@inline function count_digits(s::String)
    count = 0
    for c in s
        if isdigit(c)
            count += 1
        end
    end
    return count
end
function calc_chol_from_df(P_anc::DataFrame)

    I_vec, J_vec, V = Vector{Int64}(), Vector{Int64}(), Vector{Float64}()
    anc_cols = names(P_anc)
    anc_cols = anc_cols[count_digits.(anc_cols) .> 0]
    anc_cols = [col for col in anc_cols if !occursin("inb",col)] 
    for row in ProgressBar(eachrow(P_anc))
        i = row.new_index
        push!(I_vec,i)
        push!(J_vec,i)
        push!(V,sqrt(0.5 - 0.25 * (row.inb_0 + row.inb_1)))

        for col in anc_cols
            j = row[col]
            (j == 0) && continue
            (j > i) && error("Not ordered")

            number = match(r"\d+",col).match
            push!(I_vec, i)
            push!(J_vec, j)
            push!(V, 0.5^length(number) * row["inb_" * number])
        end
    end
    df_IJV = DataFrame(I = I_vec, J = J_vec, V = V)
    return df_IJV
end
#df_IJV = calc_chol_from_df(P_anc)
#CSV.write("./data/pedigree_IJV.csv", df_IJV)
#exit(0)

using CSV, SparseArrays, DataFrames, LinearAlgebra, Libdl
df_IJV = CSV.read("./data/pedigree_IJV.csv", DataFrame)

S = sparse(df_IJV.I, df_IJV.J, df_IJV.V, maximum(df_IJV.J), maximum(df_IJV.J)) 
ngeno = sort(unique(df_IJV.I))
S = S[ngeno, ngeno]
S = S[1:100,1:100]

# S = spdiagm(randn(10))
# S[2:10,1] .= 0.1;
I_vec, J_vec, V = findnz(S)
ncols = 1
n = S.n
b = rand(Float64,(n,ncols));
y = zeros(Float64,(n,ncols));

ccall((:chol_solve, "./src/solve_gpu.so"), Int32, (Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Int32, Int32, Ptr{Cdouble}, Int32, Ptr{Cdouble}), V, Vector{Int32}(I_vec),Vector{Int32}(J_vec), maximum(J_vec), length(I_vec), b, Int32(ncols), y)

println("Deviation ", sum(abs.(transpose(S) * y .- b)))


exit(0)

L_nn = spzeros((n, n)) 
n = size(A_nn)[1]

println("Number of nonzeros: ", size(A_nn.nzval)[1])

println("Calculating Cholesky")
# Sparse Cholesky
ccall((:cholSparse, "./solve_gpu.so"), Int32, (Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Int32, Int32, Ptr{Cdouble}, Int32, Ptr{Cdouble}), A_nn.nzval, Vector{Int32}(A_nn.colptr), Vector{Int32}(A_nn.rowval), Int32(n), Int32(size(A_nn.nzval)[1]), b, ncol, y);

println(sum(abs.(A_nn * y - b)))


