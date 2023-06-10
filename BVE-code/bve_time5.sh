#!/bin/bash
#PBS -N bve_time5
#PBS -A UMIC0093
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -j n
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.edu
#PBS -l select=2:ncpus=36:mpiprocs=36
#PBS -l place=group=rack

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

mpirun -np 72 ./driver > run_out5.txt
