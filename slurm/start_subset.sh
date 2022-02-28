#!/bin/bash
#SBATCH --job-name=dimred
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=365gb

cd $SLURM_SUBMIT_DIR

export ID=$ID

FLT3_PATS=(PARCIN PARDZW PARSGS PARUEX PARXEC PARYCG PARZVA PASBBE \
    PASDET  PASDUD PASEGC PASHWN PASIZX PASRAT PASRRB PASRTP \
    PASSUW PASTBK PASWNZ PASXJY PASXVC PATDFZ PATHZC PATIAB \
    PATKKI PATLVC PAUIIB PAUJMC PAUNVN PAUSFM PAUVEF PAUXYG \
    PAUYAE  PAUZSY PAVAGA PAVBFN PAVBZM PAVINA PAVKFD PAVNAZ \
    PAVSFB PAWAMB PAWDTX PAWTWW PAXLWH PAXMKT)

./cli/dimension_reduction.R \
    -i $ID \
    -c timepoint stemness \
    -a Diagnosis $STEMNESS \
    -S \
    -T $PANEL \
    -r $RES \
    -v DEBUG && \

#sbatch --export=ALL,ID slurm/gsva.sh && \

./cli/make_CIBERSORT_ref.R -i $ID -T "${PANEL}_SCT" && \

sbatch --export=ALL,ID slurm/CIBERSORTx_ref.sh && \

sbatch --export=ALL,ID slurm/de.sh
