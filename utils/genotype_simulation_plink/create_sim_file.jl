
using Base;
using Random;
using CSV;
using DataFrames


ROOT_DIR = string(@__DIR__) * "/../.."

n = Int(1e6)

probs = rand(n) .* 0.5 .+ 0.1
df = DataFrame(
    x = fill(1, n),
    label = [string("snp", i) for i in 1:n],
    lower = probs,
    upper = probs,
    y = fill(1, n),
    z = fill(1, n)
)

CSV.write(ROOT_DIR * "/data/geno.1m", df, delim=' ', header=false)
