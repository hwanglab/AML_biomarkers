#!/bin/bash
#SBATCH --job-name=GSEA
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8Gb

cd $SLURM_SUBMIT_DIR

./cli_test/GSEA.R -i flt3_cebpa_stem -v DEBUG