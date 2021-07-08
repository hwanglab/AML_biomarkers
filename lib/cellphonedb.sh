#!/usr/local/bin/bash

conda activate cpdb

OUTDIR=$1
PLOTDIR=$2
CORES=$3

cellphonedb method statistical_analysis \
  $OUTDIR/cellphone_db_metadata.tsv \
  $OUTDIR/cellphone_db_data \
  --output-path=$OUTDIR \
  --threads=$CORES \
  --counts-data=gene_name \
  --subsampling \
  --subsampling-log false

cellphonedb plot dot_plot \
  --output-path=$OUTDIR \
  --output-name=cellphonedb_plot.pdf \
  --means-path=$OUTDIR/means.txt \
  --pvalues-path=$OUTDIR/pvalues.txt
  
conda deactivate
