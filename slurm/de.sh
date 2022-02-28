#!/bin/bash
#SBATCH --job-name=DE
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128Gb

cd $SLURM_SUBMIT_DIR

./cli/de.R -i $ID
