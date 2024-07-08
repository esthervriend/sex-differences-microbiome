#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=16G
#SBATCH --time=120:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name male_female \
    -path /projects/0/prjs0784/sex-differences-microbiome/male_female \
    -x class \
    -n 200 \
    -t 24 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name male_female \
    -path /projects/0/prjs0784/sex-differences-microbiome/male_female \
    -x class \
    -n 200 \
    -t 24 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json \
    -permute
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name male_premen \
    -path /projects/0/prjs0784/sex-differences-microbiome/male_premen \
    -x class \
    -n 200 \
    -t 24 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json 
python /projects/0/prjs0784/sex-differences-microbiome/scripts/XGBeast.py \
    -name male_postmen \
    -path /projects/0/prjs0784/sex-differences-microbiome/male_postmen \
    -x class \
    -n 200 \
    -t 24 \
    -rand_seed 4321 \
    -param /projects/0/prjs0784/sex-differences-microbiome/scripts/param_grid.json 