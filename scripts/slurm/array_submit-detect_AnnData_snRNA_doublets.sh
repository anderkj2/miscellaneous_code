#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=detect_AnnData_snRNA_doublets-%a
#SBATCH --error=../../logs/detect_AnnData_snRNA_doublets/detect_AnnData_snRNA_doublets-%a.err
#SBATCH --output=../../logs/detect_AnnData_snRNA_doublets/detect_AnnData_snRNA_doublets-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=08:00:00
#SBATCH --array=1-76
#SBATCH --partition=compute
#SBATCH --qos=batch


# Run analysis script for each sample
bash run-detect_AnnData_snRNA_doublets.sh ${SLURM_ARRAY_TASK_ID}