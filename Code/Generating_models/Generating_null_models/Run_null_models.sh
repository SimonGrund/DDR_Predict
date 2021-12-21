#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH -c 20
#SBATCH --account pcawg
#SBATCH -t 12:00:00
#SBATCH -o job_out3/null_mod.%A.%a.txt
#SBATCH -e job_err3/null_mod.%A.%a.err
#SBATCH --array=1-67

## to run use: 
#sbatch Run_null_models.sh --array=1-67 #Nr of cases

source activate R35
#srun Rscript /home/simong/PCAWG/simong_SV/DDR_project_May_2021/Feedback_from_reviewers/Split_by_dataset/Models/Null_models_10k/Generate_null_models.R $SLURM_ARRAY_TASK_ID

Rscript Generate_null_models.R ${SLURM_ARRAY_TASK_ID}

#remainin models, 05_11: sbatch --array=1-180 Run_null_models.sh
#56783301
#08_11: 56828690 