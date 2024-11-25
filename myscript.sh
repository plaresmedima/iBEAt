#!/bin/bash
# Request 16 g
#SBATCH --mem=16G
# Request 8 cores
#SBATCH --cpus-per-task=16
#SBATCH --time=95:00:00
# Email notifications to me@somedomain.com
#SBATCH --mail-user=j.s.periquito@sheffield.ac.uk
# Email notifications if the job fails
#SBATCH --mail-type=FAIL
# Rename the job's name
#SBATCH --comment=iBEAt
#SBATCH --array=0-10
unset SLURM_CPU_BIND

export SLURM_EXPORT_ENV=ALL
module load Anaconda3/2019.07

# We assume that the conda environment 'ibeat_env_2024' has already been created
source activate ibeat_env_2024
srun python main_cluster_batch.py --num $SLURM_ARRAY_TASK_ID

