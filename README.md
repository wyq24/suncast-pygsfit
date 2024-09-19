# pygsfit_cp
# pygsfit_cp

`pygsfit_cp` is a radio imaging spectroscopy package dedicated to performing spectrum fitting. There are two ways to use this package:

## Method One (Recommended)

1. **Initial Setup and Exploration**
   - a. Run `/demo/pygsfit_cp.ipynb` to identify appropriate initial gas parameters and fitting ranges for each parameter.
   - b. Define and visualize the Region of Interest (ROI).

2. **Export Batch Script**
   - After obtaining the ideal fitting results, export a batch mode fitting script to perform the fitting with the same parameters.

## Method Two

1. **Direct Batch Mode**
   - If you are already familiar with your data, directly run `/demo/batch_mode_script.py`.
   - Manually input all the parameters of the ROI.

## Common Steps for Both Methods

2. **Execute Batch Script**
   - Execute the batch script to obtain the fitting results, which include:
     - A results file (`fitting_results.h5`)
     - An IDL script (`make_params_maps.pro`)
     - Several intermediate files.

3. **Generate Parameter Maps**
   - Run the IDL script (`make_params_maps.pro`) directly in `sswidl` to generate an IDL save file (`parameters_map.sav`).
   - The resulting save file can be opened in `gsfitview` or `sswidl`.
