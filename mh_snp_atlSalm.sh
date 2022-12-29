#!/bin/bash
# run simulations for
# performance of different SNP panel sizes
# written for execution on Ceres
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=60G # Maximum amount of memory this job will be given
#SBATCH --time=119:00:00 # ask that the job be allowed to run for 
#SBATCH --array=1-200%75 #specify how many jobs in the array and limit number running concurrently (e.g. 1-96%40)
#SBATCH --output=arrayAtlSalm_%a.out # tell it where to store the output console text
#SBATCH -p medium # partition to request

echo "My SLURM_JOB_ID: " $SLURM_JOB_ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load r/4.1.2

# check for random seeds
if [ ! -f randSeeds.txt ]; then
    echo "randSeeds.txt does not exist."
	exit 1
fi

# get random seed
x=$(cat randSeeds.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "My random seed is: " $x

# make temp directory
mkdir "$TMPDIR"/temp"$SLURM_ARRAY_TASK_ID"

# run simulation

# randomSeed iterationNumber TemporaryLocalStorageDirectory vcfInputPath
Rscript mh_snp_comparison_sim.R $x $SLURM_ARRAY_TASK_ID $TMPDIR ../seq_data_mh/allPhased_atlSalm.vcf

# remove temp directory
rm -r "$TMPDIR"/temp"$SLURM_ARRAY_TASK_ID"

echo "Done with simulation"
