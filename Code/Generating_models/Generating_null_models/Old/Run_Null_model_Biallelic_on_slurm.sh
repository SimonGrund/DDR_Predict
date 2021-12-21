#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH -c 40
#SBATCH --account pcawg
#SBATCH -t 24:00:00

source activate R35
Rscript Code/Generating_models/Generating_null_models/Null_Model_Parralel_biallelic.R

#50621030

