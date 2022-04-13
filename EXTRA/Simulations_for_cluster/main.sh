#!/bin/bash
#SBATCH --partition=dunsonlab, scavenger --account=dunsonlab
#SBATCH --output=DPMM.out
#SBATCH --job-name=DPMM
#SBATCH --mail-user=az119@duke.edu
#SBATCH --mail-type=END,FAIL
##SBATCH --array=1-8
##SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=20
module load R
R CMD BATCH --no-save main2.R out_$SLURM_ARRAY_TASK_ID
