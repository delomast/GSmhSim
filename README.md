# GSmhSim
Genomic selection simulations with microhaplotypes

This is the code used to perform a series of simulations evaluating the use of low-density microhaplotype 
panels for imputation in aquaculture selective breeding programs. 

generate_random_seeds.R is an R script that generates the random seeds used by the simulations. This was run 
once prior to running the simulations. The simulations were run on an HPC using the Slurm Workload Manager 
as array jobs (one task for each iteration). The following is an explanation of what the files are.

- mh_snp_comparison_sim.R: the R script that performed the simulations
- utils.R: utility functions used by mh_snp_comparison_sim.R
- eval_vcf.py: utilty function for evaulating a VCF file called by mh_snp_comparison_sim.R
- mh_snp_MBP.sh: the shell script that ran the Pacific oyster simulation
- mh_snp_eobc.sh: the shell script that ran the eastern oyster simulation
- mh_snp_atlSalm.sh: the shell script that ran the Atlantic salmon simulation
- microhap_GS_oyster.Rproj: file used by RStudio to organize this directory as an R project

The following are utility scripts for working on the specific HPC used for this work.

- pull_githup.sh
- runAlphaImpute2.sh
- run_blupf90.sh

archive/* contains old files not used in the final simulations


This software is a "United States Government Work" under the terms of the United States Copyright Act. It was 
written as part of the authors' official duties as United States Government employees and thus cannot be 
copyrighted. This software is freely available to the public for use. 
