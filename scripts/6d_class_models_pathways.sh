#!/bin/bash
python /media/andrei/Data/Barbara/projects/sex-differences-microbiome/scripts/XGBeast.py \
    -name sex_pathways \
    -path /media/andrei/Data/Barbara/projects/sex-differences-microbiome/sex_pathways \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /media/andrei/Data/Barbara/projects/sex-differences-microbiome/scripts/param_grid.json
python /media/andrei/Data/Barbara/projects/sex-differences-microbiome/scripts/XGBeast.py \
    -name menopause_pathways \
    -path /media/andrei/Data/Barbara/projects/sex-differences-microbiome/menopause_pathways \
    -x class \
    -n 200 \
    -t 68 \
    -rand_seed 4321 \
    -param /media/andrei/Data/Barbara/projects/sex-differences-microbiome/scripts/param_grid.json