#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=test_CellDIVE_image_DAPI_filtering-%a
#SBATCH --error=../../logs/CellDIVE_processing/test_CellDIVE_image_DAPI_filtering-%a.err
#SBATCH --output=../../logs/CellDIVE_processing/test_CellDIVE_image_DAPI_filtering-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=12:00:00
#SBATCH --array=18-22
#SBATCH --partition=compute
#SBATCH --qos=batch


# Run script for each sample
bash run-test_CellDIVE_image_DAPI_filtering.sh ${SLURM_ARRAY_TASK_ID}

