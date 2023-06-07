#!/bin/bash
#PBS -N bve_run
#PBS -A UMIC0093
#PBS -l walltime=00:10:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M cygnari@umich.com
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l place=group=rack

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

mpirun -np 1 ./driver > run_out.txt
