#!/bin/bash
#SBATCH --job-name=CIBERSORTxMixture
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=10
#SBATCH --array=0-2

cd $SLURM_SUBMIT_DIR

MIXTURE=(beat_aml.txt target_data.txt tcga_data.txt)

./cli/CIBERSORTx.py -Xi $ID -m ${MIXTURE[$SLURM_ARRAY_TASK_ID]}