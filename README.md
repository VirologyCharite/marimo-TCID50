# TCID50 calculation using marimo
This notebook calculates the 50% infectious dose (TCID50) using a generalized linear model (GLM).
You can find an online version of the notebook that does not require you to install python here: https://virologycharite.github.io/marimo-TCID50/

## Input data
Prepare your data in a spreadsheet editor in the following format:

| ID       | 10 | 100 | 1000 | 10000 | 100000 | 1000000 |
|----------|----|-----|------|-------|--------|---------|
| sample_1 | 6  | 6   | 6    | 4     | 1      | 0       |
| sample_2 | 6  | 6   | 3    | 0     | 0      | 0       |
| sample_3 | 6  | 5   | 2    | 0     |        |         |
| sample_4 |    |     | 6    | 6     | 5      | 2       |
| mock     | 0  | 0   | 0    | 0     | 0      | 0       |

Each sample has one row. The first column must be named `ID`. Every other column is one dilution: the header is the dilution factor as a linear number (1, 10, 100... not the log value 0, 1, 2) and the cells contain the number of wells with CPE. Leave the cells of dilutions that you did not do for a given sample empty - they are discarded.

You can download an example table here: https://github.com/VirologyCharite/marimo-TCID50/raw/refs/heads/main/TCID50_example.tsv

**Usage**

1. Copy the table including the header from your spreadsheet editor and paste it into the text box
2. Set the volume of virus dilution per well in µL and the number of replicates per dilution (the same for all samples and dilutions)
3. Press submit
4. Check the table that appears below the form. You can correct single values directly in this table, the results update immediately. Pressing submit again resets the table to the pasted data
5. Download the results table and the regression curves

The notebook stops with an error message if the input cannot be used, e.g. if the `ID` column is missing, a column header is duplicated, there are fewer than two dilution columns or a CPE count is higher than the number of replicates.


## Calculation
For each ID the script will attempt to fit a generalized linear model with a logit link function. It assumes the number of CPE-positive wells per dilution to follow a bimomial distribution. 
The logit of the fraction of CPE+ wells per dilution is modelled as a function of the log10 dilution. The dilution is calculated as the dilution provided in the input table multiplied with 1000/V, where V is the volume/well in µL to get to TCID50/well.
PFU/mL are calculated from TCID50/mL via the conversion factor ln(2), which is derived from the poisson distribution.
For each ID a plot is created to visualize the dose-response curve.

## Output data
1. Output table
The output table contains the following columns:
1. ID
2. detection_limit_low / detection_limit_up: The lower and upper detection limit for each ID defined as the lowest and the highest dilution
3. log_TCID50_mL: log10 transformed TCID50/mL. NaN if no titre could be calculated (see message).
4. log_PFU_mL: log_TCID50_mL+log10(ln(2))
5. TCID50_mL / PFU_mL: The same values on a linear scale
6. message: Short info about the calculation
    - below / above detection limit: None or all of the wells over all dilutions have CPE, or the regression reports a TCID50/mL outside of the detection range. The latter can occur if the lowest dilution has < 50% CPE or the highest dilution > 50% CPE.
    - no fit possible: The number of wells with CPE is the same at every dilution, or only one dilution is left for this ID
    - no titre: The number of wells with CPE does not decrease with the dilution. Check the order of your dilutions.
