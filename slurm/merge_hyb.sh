#!/bin/bash
#SBATCH --job-name=merge_GEX_HYB
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128Gb

./cli/merge_targeted_with_whole_transcriptome.r && \
    ./run_analysis.sh -i flt3_cebpa_stem_aml_res_0.8 -T AML -r 0.8 && \ 
    ./run_analysis.sh -i flt3_cebpa_stem_aml_res_1.8 -T AML -r 1.8 && \ 
    ./run_analysis.sh -i flt3_cebpa_stem_pan_cancer_res_0.8 -T PAN_CANCER -r 0.8 && \ 
    ./run_analysis.sh -i flt3_cebpa_stem_pan_cancer_res_1.8 -T PAN_CANCER -r 1.8
