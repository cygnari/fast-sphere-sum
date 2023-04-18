#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)

# Set up batch job settings
#SBATCH --job-name=direct_sum_cuda
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00
#SBATCH --account=krasny0
#SBATCH --partition=gpu

{ time ./direct_sum_cuda > direct_sum_cuda.txt ; } 2> direct_sum_cuda_time.txt
