#! /bin/bash

cd ~/Documents/Lab/clonal_evolution/biomarkers

TOKEN=e6ac1e482bcc1c8f22f80496aca50d86
EMAIL=rds156@case.edu

MIX=tcga_cibersort.txt
docker run -v $PWD/cibersort_in:/src/data \
  -v $PWD/cibersort_results:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose TRUE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $PWD/cibersort_results \
  --sigmatrix nonstem_GEP.txt \
  --mixture $MIX \
  --label tcga
  
MIX=beat_aml.txt
docker run -v $PWD/cibersort_in:/src/data \
  -v $PWD/cibersort_results:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose TRUE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $PWD/cibersort_results \
  --sigmatrix nonstem_GEP.txt \
  --mixture $MIX \
  --label beatAML
  
MIX=target_data.txt
docker run -v $PWD/cibersort_in:/src/data \
  -v $PWD/cibersort_results:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose TRUE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $PWD/cibersort_results \
  --sigmatrix nonstem_GEP.txt \
  --mixture $MIX \
  --label target
