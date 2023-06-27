# run AlphaImpute2
# load python, activate virutal environment
# 1 output file prefix
# 2 genotype input
# 3 pedigree input
# 4 random seed
# 5 max thread for imputation
module load miniconda
source activate /project/oyster_gs_sim/ai2_2/

AlphaImpute2 -out $1 -genotypes $2 -pedigree $3 -seed $4 -maxthreads $5 # -ped_only

# quick modification to allow more cores to be used without a lot of programming
# AlphaImpute2 -out $1 -genotypes $2 -pedigree $3 -seed $4 -maxthreads $SLURM_JOB_CPUS_PER_NODE # -ped_only
