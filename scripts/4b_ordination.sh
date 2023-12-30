#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate rscripts
Rscript /projects/0/prjs0784/sex-differences-microbiome/scripts/4a_ordinationplots.R