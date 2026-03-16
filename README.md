# Miscellaneous Code

Miscellaneous code snippets and functions for data processing and/or analysis

### Note:
Included bash scripts associated with a given code script are for running code as a SLURM job array. The initial submission script includes the prefix `array_submit-`, which runs multiple instances of the script with the prefix `run_`. The `run_` script uses an input table containing sample IDs and data parameters to populate required command-line arguments for a Python or R script. 
