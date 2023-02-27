#!/bin/bash
# (See https://arc-ts.umich.edu/greatlakes/user-guide/ for command details)

# Set up batch job settings
#SBATCH --job-name=fast_sum_mpi
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=1g
#SBATCH --time=00:05:00
#SBATCH --account=krasny0
#SBATCH --partition=standard


# Run your program

# (">" redirects the print output of your program, in this case to "output.txt")

# -np command sets number of processors to run program with. If you specify more processors
# using -np command than you have requisitioned using the #SBATCH --ntasks-per-node command
# your program will not run

# ./hello_world is the name of the program executable

# --bind-to core sets up mpi environment and specifies which hardware to bind MPI processes to
# you don't need to mess with this

mpirun -np 2 --bind-to core:overload-allowed ./fast_sum_mpi  > fast_sum_mpi_2.txt
