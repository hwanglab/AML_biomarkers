#!/bin/bash
#SBATCH --job-name=GSVA
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=5
#SBATCH --mem=64Gb

cd $SLURM_SUBMIT_DIR

./cli/run_gsva.R -i $ID 
