#!/bin/bash
#SBATCH -p gpu_4_h100
#SBATCH --gres=gpu:1
#SBATCH -t 2:00:00
#SBATCH --exclusive
#SBATCH -c 40
#SBATCH --mem=230gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.freudenberg@uni-mannheim.de
#SBATCH --output=BENCHMARK_%x_%j.out

nvidia-smi

population_sizes=(small medium large)
module load devel/cuda/12
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