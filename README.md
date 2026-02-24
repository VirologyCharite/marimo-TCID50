# TCID50 calculation using marimo
This notebook calculates the 50% infectious dose (TCID50) using a generalized linear model (GLM).
You can find an online version of the notebook that does not require you to install python here: https://virologycharite.github.io/marimo-TCID50/

## Input data
Prepare your data in a spreadsheet editor in the following format:

![excel table template](https://github.com/VirologyCharite/marimo-TCID50/raw/refs/heads/main/docs/table.png)

Each sample has one row. For each pre-dilution there are two columns: One in which you put the number of wells with CPE and another in which you put the total number of replicates (usually the same across an experiment). The script will discard columns with no CPE value - so just leave empty dilutions that you did not do for a given sample

You can download a sample sheet here: https://github.com/VirologyCharite/marimo-TCID50/raw/refs/heads/main/TCID50.xlsx

**Settings**
In the notebook you can further set the volume in µl of virus dilution per well and the decimal separator for input 

** Instructions in one image: **

![gui overview](https://github.com/VirologyCharite/marimo-TCID50/raw/refs/heads/main/docs/gui.jpeg)


## Calculation
For each ID the script will attempt to fit a generalized linear model with a logit link function. It assumes the number of CPE-positive wells per dilution to follow a bimomial distribution. 
The logit of the fraction of CPE+ wells per dilution is modelled as a function of the log10 dilution. The dilution is calculated as the dilution provided in the input table multiplied with 1000/V, where V is the volume/well in µL to get to TCID50/well.
PFU/mL are calculated from TCID50/mL via the conversion factor ln(2), which is derived from the poisson distribution.
For each ID a plot is created to visualize the dose-response curve.

## Output data
1. Output table
The output table contains the following columns:
1. ID
2. log_TCID50_mL: log10 transformed TCID50/mL. NaN if all or none of the wells over all dilutions have CPE.
3. detection_limit_low / detection_limit_up: The lower and upper detection limit for each ID defined as the lowest and the highest dilution
4. message: Short info about the calculation
5. log_PFU_mL: log_TCID50_mL+log10(ln(2))
The regression can occasionally report a TCID50/mL outside of the detection range. This can occur if the highest dilution has < 50% CPE or the highest dilution > 50% CPE. 
