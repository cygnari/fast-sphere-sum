#!/bin/bash
#PBS -N bve_time_all
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -j n
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l place=group=rack

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

mpirun -np 16 ./driver > run_out0.txt
mpirun -np 16 ./driver > run_out1.txt
mpirun -np 16 ./driver > run_out2.txt
mpirun -np 16 ./driver > run_out3.txt
mpirun -np 16 ./driver > run_out4.txt
mpirun -np 16 ./driver > run_out5.txt
mpirun -np 16 ./driver > run_out6.txt
mpirun -np 16 ./driver > run_out7.txt
mpirun -np 16 ./driver > run_out8.txt
mpirun -np 16 ./driver > run_out9.txt