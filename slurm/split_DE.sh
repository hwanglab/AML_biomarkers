#!/bin/bash
#SBATCH --job-name=DE
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64Gb

cd $SLURM_SUBMIT_DIR

./cli/split_DE_clustering.R -i flt3_cebpa_stem_bc -S -v DEBUG && sbatch slurm/gsea.sh