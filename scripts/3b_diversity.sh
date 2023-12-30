#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
Rscript -e /projects/0/prjs0784/sex-differences-microbiome/scripts/3a_diversityplots.R