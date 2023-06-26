#!/bin/bash
#SBATCH -t 3:00:00
##SBATCH --exclusive
#SBATCH --mem=160gb
#SBATCH --mail-type=ALL
#SBATCH --output=BENCHMARK_%x_%j.out

export OMP_NUM_THREADS=64

# Start Julia benchmark script
# Command line needs to specify "miraculix" or "PLINK"
julia utils/benchmark/run_suite.jl $1