#!/bin/bash

#SBATCH --get-user-env
#SBATCH --job-name=test_CellDIVE_image_DAPI_filtering-%a
#SBATCH --chdir=/projects/varn-lab/USERS/anderk/miscellaneous_code/scripts
#SBATCH --error=/projects/varn-lab/USERS/anderk/miscellaneous_code/logs/test_CellDIVE_image_DAPI_filtering-%a.err
#SBATCH --output=/projects/varn-lab/USERS/anderk/miscellaneous_code/logs/test_CellDIVE_image_DAPI_filtering-%a.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Kevin.Anderson@jax.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768G
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=12:00:00
#SBATCH --array=18-22
#SBATCH --partition=compute
#SBATCH --qos=batch


# Run merge script.
bash run_test_CellDIVE_image_DAPI_filtering.sh ${SLURM_ARRAY_TASK_ID}

