#!/bin/bash
#SBATCH --job-name=survival
#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=1

./cli/survival_analysis.R \
    -i $ID \
    -I EFS -V 'c(matches("^Year|^Age|^Overall"), where(is_character))' \
    -a 'Event Free Survival Time in Days' -E
