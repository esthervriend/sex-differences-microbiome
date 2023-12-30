#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=16G
#SBATCH --time=50:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name menopause \
    -path /projects/0/prjs0784/sex-differences-microbiome/menopause \
    -x class \
    -n 200 \
    -t 24 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name menopause \
    -path /projects/0/prjs0784/sex-differences-microbiome/menopause \
    -x class \
    -n 200 \
    -t 24 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json \
    -permute