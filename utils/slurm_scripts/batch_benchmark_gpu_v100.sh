#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 1:00:00
#SBATCH -c 10
#SBATCH --mem=230gb
#SBATCH --mail-type=ALL
#SBATCH --output=BENCHMARK_%x_%j.out
module load devel/cuda/12
nvidia-smi

population_sizes=(small medium large)
export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10
    
for size in "${population_sizes[@]}"
do
    echo ""
    echo "Number of cores: $MKL_NUM_THREADS, $OMP_NUM_THREADS"
    echo "Population size: $size"
    echo "Execution mode: GPU"
    ./benchmark.out GPU ../data/$size.bed ../data/$size.freq
done