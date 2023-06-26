#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 2:00:00
#SBATCH --exclusive
#SBATCH -c 40
#SBATCH --mem=230gb
#SBATCH --mail-type=ALL
#SBATCH --output=BENCHMARK_%x_%j.out

nvidia-smi
module load devel/cuda/12

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