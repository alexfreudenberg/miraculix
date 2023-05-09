#!/bin/bash
#SBATCH -p single
#SBATCH -t 9:00:00

#SBATCH -c 40
#SBATCH --mem=170gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.freudenberg@uni-mannheim.de
#SBATCH --output=BENCHMARK_%x_%j.out

module load devel/cuda/12
core_numbers=(20 10 1)
population_sizes=(large medium small)
execution_modes=(CPU Org)

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