#!/bin/bash
#SBATCH --job-name=seapodym_cohort
#SBATCH --time=00:10:00           # walltime (adjust as needed)
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=51
#SBATCH --mem=25g
#SBATCH --cpus-per-task=1
#SBATCH --extra-node-info=1:*:1 
#SBATCH --distribution=*:block:* 
#SBATCH --mem-bind=local

module purge
module load gimkl/2020a

rm log_*.txt

# Run the MPI job
srun ../../bin/seapodym_cohort -s skj_fat.xml
