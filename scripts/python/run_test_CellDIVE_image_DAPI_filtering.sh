#!/bin/bash

### Takes pyramidal OME-TIFF and QuPath measurements as input; runs QC checks for DAPI cycles, and ###
### performs filtering based on DAPI and morphology metrics (and any previous visual inspection)   ### 


### Input arguments ###
# Define the image file parameters from input file.
# NOTE: linker file has column names, so must add 1 to array ID
ARRAYID="`expr ${SLURM_ARRAY_TASK_ID} + 1`"
LINKER="/projects/varn-lab/PvR_TME/metadata/CellDIVE_panel_sample_list.txt"

# Sample ID
sample_ID=$(sed "${ARRAYID}q;d" ${LINKER} | cut -f1)

# CellDIVE run ID
CellDIVE_ID=$(sed "${ARRAYID}q;d" ${LINKER} | cut -f3)

# CellDIVE directory containing individual image files
CellDIVE_dir=$(sed "${ARRAYID}q;d" ${LINKER} | cut -f6)

# CellDIVE single channel image (used for extracting physical pixel size)
CellDIVE_image=${CellDIVE_dir}raw/${CellDIVE_ID}_1.0.4_R000_DAPI__FINAL_F.ome.tif

# CellDIVE image physical pixel size (must be extracted from single channel image using tiffcomment)
CellDIVE_p=$(/projects/varn-lab/USERS/anderk/anderk_apps/Bio-Formats/bftools/tiffcomment  ${CellDIVE_image} | grep PhysicalSizeX | cut -d ' ' -f13 | cut -d "\"" -f2 | xargs printf "%.4f\n")

# Ouptut directory (will be created if not existing)
outdir=/projects/varn-lab/PvR_TME/results/FOV-corrected_cohort/${sample_ID}/
mkdir -p ${outdir}
###


STARTTIME=`date`
echo $STARTTIME
echo "Processing CellDIVE image for $sample_ID"
echo ""

python test_CellDIVE_image_DAPI_filtering.py $sample_ID $CellDIVE_ID $CellDIVE_dir --pixel_size $CellDIVE_p $outdir

echo "CellDIVE processing completed at: `date`"



