'''
test_CellDIVE_image_DAPI_filtering.py
Takes pyramidal OME-TIFF and QuPath measurements as input; 
runs QC checks for DAPI cycles, and performs filtering based on DAPI and 
morphology metrics (and any previous visual inspection)
'''


import os
import argparse
import re
import glob
import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
from anndata import AnnData
import h5py
import hdf5plugin
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import tifffile
import json
import geojson
from pytrendseries import detecttrend
from pytrendseries import vizplot
import phenograph
from kneed import KneeLocator
from natsort import index_natsorted, natsort_keygen, natsorted
from utils import *


# Restrict the number of warnings displayed in output
import warnings
warnings.filterwarnings(action='once')



##### Parse command line arguments #####
parser = argparse.ArgumentParser(description="test_CellDIVE_image_DAPI_filtering.py - Takes pyramidal OME-TIFF and QuPath measurements as input; runs QC checks for DAPI cycles, and performs filtering based on DAPI and morphology metrics (and any previous visual inspection)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("sample", help="Sample ID to be used as prefix for output files")
parser.add_argument("run", help="CellDIVE run ID, used as a prefix for some input files")
parser.add_argument("indir", help="Directory containing metadata and raw image files for CellDIVE run '/')")
parser.add_argument("-p", "--pixel_size", type=float, action="store", help="Image physical pixel size (in um)")
parser.add_argument("outdir", help="Output directory (should include trailing '/')")

args = vars(parser.parse_args())
########################################



##### Input arguments #####
# Set parameters from command line
sample = args["sample"]
run = args["run"]
indir = args["indir"]
pixel_size = args["pixel_size"]
outdir = args["outdir"]

# Set directory for formatted data
omedir = '/projects/varn-lab/PvR_TME/formatted_data/' + sample + '/'

# Create directory for downstream analysis if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Set results directory as working directory
os.chdir(outdir)

# Extract panel ID from metadata json file
metadata_file = indir + 'design/' + 'metadata.json'

with open(metadata_file, "r") as file:
    metadata = json.load(file)

CellDIVE_panel = metadata['panel_id']

# Load formatted OME-TIFF 
ome_img = tifffile.imread(omedir + 'CellDIVE_panel_' + str(CellDIVE_panel) + '_all_cycles_aligned-FOV-corrected-' + sample + '.ome.tif')

print('Image size: ' + str(ome_img.shape))
print('')

# Load QuPath-generated measurements for annotated tissue regions
tissue_meta = pd.read_csv(omedir + sample + '_tissue_annotations.txt', sep='\t')

# Load QuPath-generated measurements for detected cells
cell_meta = pd.read_csv(omedir + sample + '_FOV-corrected_cell_detections.txt', sep='\t')

# Load metadata table containing clinical timepoint info
clinical_meta = pd.read_csv('/projects/varn-lab/PvR_TME/metadata/CellDIVE_sample_clinical_metadata-formatted.txt', sep='\t')

# Load sorted image channel info
img_info = pd.read_csv(omedir + sample + '_sorted_image_channel_info.csv', index_col=0)
###########################



##### Check DAPI quality and morphology metrics for data filtering prior to AnnData generation #####
### Compare DAPI intensities across imaging cycles
# Subset cell detections to relevant columns for plotting
tmp = cell_meta[['Detection probability','Nucleus: DAPI_1.0.4: Mean','Nucleus: DAPI_6.0.4: Mean']]

# Sort data by cell detection probability
tmp = tmp.sort_values(by=['Detection probability','Nucleus: DAPI_1.0.4: Mean'])

# Reset data frame index
tmp = tmp.reset_index(drop = True)

# Save the detection probability column and remove from dataframe
tmp_p = tmp['Detection probability']

tmp = tmp.drop(columns = ['Detection probability'])



# Initialize figure
fig = plt.figure(figsize=(6, 4), dpi=200)

# Plot density plots of mean nucleus DAPI intensity
plot = sns.kdeplot(tmp)
plot.set_xlabel('Pixel intensity')

# Show the plot
plt.show()

# Save the plot
fig.savefig(outdir + sample + '_mean_nucleus_DAPI_cycle_comparisons_density_plot.png')
plt.close()



# Initialize figure
fig = plt.figure(figsize=(6, 4), dpi=200)

# Plot histograms of mean nucleus DAPI intensity
plot = sns.histplot(tmp)
plot.set_xlabel('Pixel intensity')

# Show the plot
plt.show()

# Save the plot
fig.savefig(outdir + sample + '_mean_nucleus_DAPI_cycle_comparisons_histogram.png')
plt.close()
###


### Determine threshold for filtering based on DAPI cycle
# Define a logistic growth function
def logistic_func(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))

# Define a sum of squared error function to use for parameter optimization
def sum_of_squared_error(parameter_set):
    val = logistic_func(x_data, *parameter_set)
    return np.sum((y_data - val) ** 2.0)

# Define a function to run parameter optimization
def generate_initial_parameters():
    parameter_bounds = []
    parameter_bounds.append([0.0, np.mean(y_data[-3:])]) # search bounds for L
    parameter_bounds.append([np.max(y_data)/2, np.max(y_data)*2]) # search bounds for k
    parameter_bounds.append([np.min(y_data)/2, np.min(y_data)/2]) # search bounds for x0

    # Set a fixed seed for parameter optimization troubleshooting
    result = sp.optimize.differential_evolution(sum_of_squared_error, parameter_bounds, seed=5)
    return result.x


# Initialize output list
dapi_threshold = []

# Initialize plot with multiple subpanels
rp = 1  # num row panels
cp = 2  # num col panels
fig, ax = plt.subplots(rp, cp, figsize=(8, 4),dpi=200)

# Create lists of x,y indices for subpanels
ax_x = np.concatenate([([i]*cp) for i in np.arange(0,rp)], axis=0)
ax_y = np.array(list(np.arange(0,cp)) * rp)


# Loop through DAPI cycles for analysis
for i, dapi in enumerate(tmp.columns):
    # Generate histogram of mean nucleus DAPI intensity
    hist = sns.histplot(tmp[dapi], ax=ax[ax_y[i]])
    
    # Extract the x-coordinates for histogram bars (corresponds to left bottom corner)
    x_coords = [patch.get_x() for patch in hist.patches]

    # Extract the y-coordinates for histogram bars (corresponds to top bar height)
    y_coords = [patch.get_height() for patch in hist.patches]

    
    # Create a dataframe of the initial histogram bin x-coordinates and the correlated date-time
    # NOTE: when input to pd.to_datetime is an array of integers, 
    # the output is parsed as the number of units (defined by "unit") since the reference date
    time_conversion = pd.DataFrame({'X':x_coords,
                                    'DT':pd.to_datetime(x_coords, unit = 'D')})

    
    # Create formatted input for detecttrend functions
    # NOTE: to prevent errors, adding pseudocount so that there are no 0 "prices"
    trend_input = pd.DataFrame({'period':time_conversion['DT'],
                                'close_price':y_coords}).set_index('period')

    trend_input['close_price'] = trend_input['close_price'] + 1

    # Isolate region(s) of curve where intensity decreases by calling detecttrend function
    trends_detected = detecttrend(trend_input, trend='downtrend', limit=2, window=1000)

    # Convert start time periods back to x-coordinates
    check_from = time_conversion.merge(trends_detected[['from']], left_on='DT', right_on='from')

    # Convert end time periods back to x-coordinates
    check_to = time_conversion.merge(trends_detected[['to']], left_on='DT', right_on='to')

    
    # Calculate peak(s) of histogram curve
    peaks, properties = sp.signal.find_peaks(y_coords, height=max(y_coords)/2)

    # Identify the x-coordinate for the tallest peak
    peak_x = x_coords[peaks[properties['peak_heights'] == max(properties['peak_heights'])][0]]

    # Identify the y-coordinate for the tallest peak
    peak_y = y_coords[peaks[properties['peak_heights'] == max(properties['peak_heights'])][0]]

    # Identify the y-coordinate for the first downward trend
    first_trend_y = y_coords[np.where(x_coords == check_to.iloc[0,0])[0][0]]
    
    # If the first downward trend is before the tallest peak, extract the region between the downward trend and the peak
    # as the space to search for filtering cutoff
    # If the first downward trend is at/after the tallest peak, extract the region before the peak as the search space
    if ((peak_x > check_from.iloc[0,0]) & (abs(peak_y-first_trend_y)/peak_y > 0.1)):
        trend_start_ind = time_conversion.loc[time_conversion['X'] == check_to.iloc[0,0]].index.values[0]
        trend_stop_ind = time_conversion.loc[time_conversion['X'] == check_from.iloc[1,0]].index.values[0]
    else:
        trend_start_ind = 0
        trend_stop_ind = time_conversion.loc[time_conversion['X'] == check_from.iloc[0,0]].index.values[0]

    # Set x and y coordinates for identified upward trend
    x_data = x_coords[trend_start_ind:trend_stop_ind]
    y_data = y_coords[trend_start_ind:trend_stop_ind]

    # If the upward trend does not have at least 3 data points, set the search space as the first 3 points
    if (len(x_data) < 3):
        x_data = x_coords[0:3]
        y_data = x_coords[0:3]

    # If the first downward trend is before the tallest peak, identify the threshold as the point of maximum curvature 
    # of a logistic growth curve fit to extracted region between downward trend and peak
    # If the first downward trend is at/after the tallest peak, identify the threshold as the knee point of the curve
    # for the region before the peak
    if ((peak_x > check_from.iloc[0,0]) & (abs(peak_y-first_trend_y)/peak_y > 0.1)):
        # Run parameter optimization for curve fitting
        initial_params = generate_initial_parameters()

        # Fit a logistic curve to identified upward trend using optimized input parameters
        params, _ = sp.optimize.curve_fit(logistic_func, 
                                          x_data, 
                                          y_data,
                                          p0=initial_params)

        # Extract fitted curve parameters
        L, k, x0 = params

        # Calculate the first and second derivatives of the fitted curve
        first_derivative = (L * k * np.exp(-k * (x_data - x0))) / (1 + np.exp(-k * (x_data - x0)))**2
        second_derivative = (L * k**2 * np.exp(-k * (x_data - x0)) * (np.exp(-k * (x_data - x0)) - 1)) / (1 + np.exp(-k * (x_data - x0)))**3

        # Calculate the curvature
        curvature = np.abs(second_derivative) / (1 + first_derivative**2)**1.5

        # Find the x value where the curvature is greatest
        max_curvature_x = x_data[np.argmax(curvature)]

        # Set DAPI threshold as the point of maximum curvature
        dapi_threshold.append(max_curvature_x)
    else:
        # Calculate the knee point of curve using kneed package
        kn = KneeLocator(x_data, y_data, curve='convex', direction='increasing', S=0.1, online=True)
    
        # Set DAPI threshold as knee point
        dapi_threshold.append(kn.knee)

    # Add vertical line representing DAPI filtering threshold to histogram plot
    ax[ax_y[i]].axvline(x=dapi_threshold[i], color='r', linestyle='--', linewidth=1)

    # Set the font size for histogram axis labels and ticks
    ax[ax_y[i]].set_xlabel(dapi, fontsize=10)
    ax[ax_y[i]].set_ylabel('Count', fontsize=10)
    ax[ax_y[i]].tick_params(axis='both', which='major', labelsize=8)

# Write final plot to file
plt.tight_layout()
plt.show()
fig.savefig(outdir + sample + '_mean_nucleus_DAPI_cycle_histogram_filtering.png')
plt.close()

# Print DAPI threhsolds to output log file
print('DAPI thresholds: ' + str(dapi_threshold))
print('')

# Write DAPI thresholds to file
np.save(outdir + sample + '_mean_nucleus_DAPI_cycle_thresholds.npy',np.array(dapi_threshold))
###


### Filter cells based on DAPI intensity and morphology metrics
# Generate scatterplot showing relationship between cell circularity, nucleus circularity, and detection probability
with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": (200)}):
    sns.scatterplot(data=cell_meta, x='Nucleus: Circularity', y='Cell: Circularity', hue='Detection probability')
    plt.show()
    plt.savefig(outdir + 'CellDIVE_panel_' + str(CellDIVE_panel) + '_all_cycles_aligned-' + sample + '_cell_vs_nucleus_circularity_scatterplot.png')
    plt.close()

# Check the number of cells in data
print('Cells before filtering: ' + str(cell_meta.shape))
print('')

# Filter data
cell_meta = cell_meta[(cell_meta['Nucleus: Circularity'] >= 0.65) & 
                      (cell_meta['Cell: Circularity'] >= 0.60) & 
                      (cell_meta['Nucleus: DAPI_1.0.4: Mean'] >= dapi_threshold[0])]

# Check cell count after filtering
print('Cells after filtering: ' + str(cell_meta.shape))
print('')
###
####################################################################################################