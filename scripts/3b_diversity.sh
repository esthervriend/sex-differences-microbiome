#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=7G
#SBATCH --time=0:30:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate rscripts
Rscript /projects/0/prjs0784/sex-differences-microbiome/scripts/3a_diversityplots.R