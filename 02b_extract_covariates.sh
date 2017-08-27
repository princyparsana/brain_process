#!/bin/sh
all_cov_fn="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.brain.txt"
all_cov_pc_fn="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt"
std_cov_pc_fn="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.std_covars.PCs.brain.txt"

python 02_extract_covariates.py -all_cov_fn $all_cov_fn -all_cov_pc_fn $all_cov_pc_fn -std_cov_pc_fn $std_cov_pc_fn
