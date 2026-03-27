#!/bin/bash
#SBATCH --job-name=seapodym_cohort
#SBATCH --time=00:10:00           # walltime (adjust as needed)
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=20g
#SBATCH --cpus-per-task=1
#SBATCH --output=seapodym_%j.out
#SBATCH --error=seapodym_%j.err

module purge
module load gimkl/2020a

# Run the MPI job
srun ../../bin/seapodym_cohort -s skj_fat.xml
