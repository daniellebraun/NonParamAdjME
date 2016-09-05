#!/bin/bash
#SBATCH -J tophat #A single job name for the array
#SBATCH -o tophat_%A_%a.out #Standard output
#SBATCH -e tophat_%A_%a.err #Standard error
#SBATCH -p serial_requeue#Partition
#SBATCH -t 120 #Running time of 2 hours
#SBATCH --mem-per-cpu 10000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine


R CMD BATCH Run_functionFast$SLURM_ARRAY_TASK_ID.R
