#!/bin/bash
#SBATCH -p fat
#SBATCH -t 6:00:00
###SBATCH --exclusive
#SBATCH -c 40
#SBATCH --mem=2500gb
#SBATCH --mail-type=ALL
#SBATCH --output=BENCHMARK_%x_%j.out

core_numbers=(20 10 1)
population_sizes=(small medium large)
execution_modes=(MKL)

for cores in "${core_numbers[@]}"
do
    export OMP_NUM_THREADS=$cores
    export MKL_NUM_THREADS=$cores
    
    for size in "${population_sizes[@]}"
    do
        for mode in "${execution_modes[@]}"
        do
            echo ""
            echo "Number of cores: $MKL_NUM_THREADS, $OMP_NUM_THREADS"
            echo "Population size: $size"
            echo "Execution mode: $mode"
            ./benchmark.out $mode ../data/$size.bed ../data/$size.freq
        done
    done
done