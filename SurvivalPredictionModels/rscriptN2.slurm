#!/bin/bash

#SBATCH -n 1               # Number of cores requested
#SBATCH -N 1               # Ensure that all cores are on one machine
#SBATCH -t 120              # Runtime in minutes
#SBATCH -p serial_requeue  # Partition to submit to
#SBATCH --mem-per-cpu=10000  # Memory per cpu in MB (see also --mem)
#SBATCH -o hostname.out    # Standard out goes to this file
#SBATCH -e hostname.err    # Standard err goes to this file
R CMD BATCH generaten2.R

