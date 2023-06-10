#!/bin/bash
#PBS -N bve_time6
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -j n
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=4:ncpus=36:mpiprocs=36
#PBS -l place=group=rack

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

mpirun -np 144 ./driver > run_out6.txt
