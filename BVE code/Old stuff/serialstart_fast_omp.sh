#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)

# Set up batch job settings
#SBATCH --job-name=fast_sum_omp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --account=eecs587f22_class
#SBATCH --partition=standard

export OMP_NUM_THREADS=1

{ time ./fast_sum_omp > fast_sum_omp_1.txt ; } 2> fast_sum_omp_1_time.txt
