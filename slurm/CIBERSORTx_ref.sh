#!/bin/bash
#SBATCH --job-name=CIBERSORTxREF
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=10

cd $SLURM_SUBMIT_DIR

./cli/CIBERSORTx_ref.py \
    -i $ID \
    -c 10

if [ $? -eq 0 ]; then
   JOBID=$(sbatch --export=ALL,ID --parsable slurm/CIBERSORTx.sh)
   sbatch --dependency=afterok:$JOBID --export=ALL,ID slurm/survival.sh
fi
