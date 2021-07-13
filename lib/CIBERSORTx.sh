#!/usr/bin/env bash

cd ~/Documents/Lab/clonal_evolution/biomarkers || exit

sed 1d cibersort_in/nonstem_clusters.txt > cibersort_in/nonstem_clusters_ready.txt

./lib/RunCIBERSORT.sh tcga_cibersort.txt
./lib/RunCIBERSORT.sh beat_aml.txt
./lib/RunCIBERSORT.sh target_data.txt
