#! /bin/bash

cd ~/Documents/Lab/clonal_evolution/biomarkers || exit

TOKEN=e6ac1e482bcc1c8f22f80496aca50d86
EMAIL=rds156@case.edu

OUTDIR=$PWD"/outs/cibersort_results"
SIGMAT=nonstem_clusters_sigmatrix.txt

MIX=tcga_cibersort.txt
docker run -v $PWD/cibersort_in:/src/data \
  -v $OUTDIR:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose TRUE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $OUTDIR \
  --sigmatrix $SIGMAT \
  --mixture $MIX \
  --label tcga
  
MIX=beat_aml.txt
docker run -v $PWD/cibersort_in:/src/data \
  -v $OUTDIR:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose TRUE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $OUTDIR \
  --sigmatrix $SIGMAT \
  --mixture $MIX \
  --label beatAML
  
MIX=target_data.txt
docker run -v $PWD/cibersort_in:/src/data \
  -v $OUTDIR:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose TRUE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $OUTDIR \
  --sigmatrix $SIGMAT \
  --mixture $MIX \
  --label target

