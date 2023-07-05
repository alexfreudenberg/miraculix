#!/bin/bash
#SBATCH -t 3:00:00
##SBATCH --exclusive
#SBATCH --mem=160gb
#SBATCH --mail-type=ALL
#SBATCH --output=BENCHMARK_julia_%x_%j.out

# Start Julia benchmark script
# Command line needs to specify "miraculix" or "PLINK"
julia --threads=$OMP_NUM_THREADS utils/benchmark/run_suite.jl $1 $2