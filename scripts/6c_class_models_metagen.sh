#!/bin/bash
#SBATCH -c 68
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name sex_metagen \
    -path /projects/0/prjs0784/sex-differences-microbiome/sex_metagen \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name menopause_metagen \
    -path /projects/0/prjs0784/sex-differences-microbiome/menopause_metagen \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json