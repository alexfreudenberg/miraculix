#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 30:00
##SBATCH --exclusive
#SBATCH -c 8
#SBATCH --mem=60gb
##SBATCH --mail-type=ALL
#SBATCH --output=BENCHMARK_%x_%j.out

export OMP_NUM_THREADS=64

julia utils/benchmark/run_suite.jl