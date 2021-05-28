#! /bin/bash

# local vars

TOKEN=e6ac1e482bcc1c8f22f80496aca50d86
EMAIL=rds156@case.edu

INDIR=$PWD/cibersort_in
OUTDIR=$PWD/outs/cibersort_results

SIGMAT=nonstem_clusters_GEP.txt

LAB="${1%.*}"

# INVOCATION
# ./lib/RunCIBERSORT MIXTURE

echo Running CIBERSORTxFractions on $LAB...

docker run -v $INDIR:/src/data \
  -v $OUTDIR:/src/outdir \
  cibersortx/fractions \
  --username $EMAIL \
  --verbose FALSE \
  --token $TOKEN \
  --single_cell TRUE \
  --outdir $OUTDIR \
  --sigmatrix $SIGMAT \
  --mixture $1 \
  --label $LAB 1> /dev/null
