#!/usr/local/bin/bash

cd ~/Documents/Lab/clonal_evolution/biomarkers || exit

conda activate cpdb

INDIR=data/cellphonedb_in
OUTDIR=outs/cellphonedb_results
PLOTDIR=plots

cellphonedb method statistical_analysis \
  $INDIR/metadata.tsv \
  $INDIR/data \
  --output-path=$OUTDIR \
  --threads=12 \
  --counts-data=gene_name \
  --subsampling \
  --subsampling-log false

cellphonedb plot dot_plot \
  --output-path=$OUTDIR \
  --output-name=cellphonedb_plot.pdf \
  --means-path=$OUTDIR/means.txt \
  --pvalues-path=$OUTDIR/pvalues.txt
  
conda deactivate
