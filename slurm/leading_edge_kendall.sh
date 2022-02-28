#!/bin/bash
#SBATCH --job-name=cor_leading_edge
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32Gb

ID=flt3_cebpa_stem_bc
GENES=(ZEB2 ID2 FOS FOSB JUNB JUND KLF6 NPM1)
PATHWAYS=(TAF9B_TARGET_GENES FOXE1_TARGET_GENES GTF2E2_TARGET_GENES)

for set in $PATHWAYS
do
   ./cli/gsea_leading_edge_plot.R \
      -v DEBUG \
      -i $ID \
      -G $set \
      -D PvF
done

./cli/Kendall.r -i $ID --genes ${GENES[@]}