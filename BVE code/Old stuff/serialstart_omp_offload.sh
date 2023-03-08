#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)

# Set up batch job settings
#SBATCH --job-name=direct_sum_omp_offload
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00
#SBATCH --account=eecs587f22_class
#SBATCH --partition=gpu

export OMP_NUM_THREADS=1
export OMP_TARGET_OFFLOAD=MANDATORY

{ time ./direct_sum_omp_offload > direct_sum_omp_off.txt ; } 2> direct_sum_omp_off_time.txt
