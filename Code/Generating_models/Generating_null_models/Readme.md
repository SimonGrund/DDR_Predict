# Generating null models

This folder contains the code to permute the mutation status of the gene for each of the 535 DDR gene deficiency models. The code is designed to generate ~10.000 data permutation models for each of the 535 cases, and takes around 24 hours with 60 CPU cores, 1gb ram each. The code is designed for parallel running, and may be started on a slurm setup using the .sh script:  Run_null_models.sh

To re-run analysis:

1) first generate the 535 cases using "Generate_cases.R"
2) Generate ~10.000 null-models per case using "Run_null_models.sh""
3) update the case list by un-hashing the code in the bottom of Generate_cases.R and updating the case list accordingly
4) re-do 2 and 3 till you are happy with the number of null models. 

The generated null-models and their performance matrix is in output/null_models/

