#!/bin/bash
#SBATCH --job-name=seapodym_cohort
#SBATCH --time=00:10:00           # walltime (adjust as needed)
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=4g
#SBATCH --cpus-per-task=1
#SBATCH --output=seapodym_%j.out
#SBATCH --error=seapodym_%j.err
##SBATCH --profile task
##SBATCH --acctg-freq=1
#SBATCH --extra-node-info=1:*:1 
#SBATCH --distribution=*:block:* 
#SBATCH --mem-bind=local

module purge
module load gimkl/2020a

rm log_*.txt

# Run the MPI job
srun ../../bin/seapodym_cohort -s skj_fat.xml
