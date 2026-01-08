# TCID50 calculation using marimo
This notebook calculates the 50% infectious dose (TCID50) using a generalized linear model (GLM).

## Input data
Prepare your data in a spreadsheet editor in the following format:
| ID       | Dilution | CPE | Replicates |   |
|----------|----------|-----|------------|---|
| Sample 1 | 1.00E+02 | 8   | 8          |   |
| Sample 1 | 1.00E+03 | 8   | 8          |   |
| Sample 1 | 1.00E+04 | 8   | 8          |   |
| Sample 1 | 1.00E+05 | 6   | 8          |   |
| Sample 1 | 1.00E+06 | 2   | 8          |   |
| Sample 1 | 1.00E+07 | 0   | 8          |   |
| Sample 2 | 1.00E+02 | 8   | 8          |   |
| Sample 2 | 1.00E+03 | 4   | 8          |   |
| Sample 2 | 1.00E+04 | 0   | 8          |   |
| Sample 2 | 1.00E+05 | 0   | 8          |   |
| Sample 2 | 1.00E+06 | 0   | 8          |   |
| Sample 2 | 1.00E+07 | 0   | 8          |   |
| Sample 3 | 1.00E+02 | 8   | 8          |   |
| Sample 3 | 1.00E+03 | 8   | 8          |   |
| Sample 3 | 1.00E+04 | 8   | 8          |   |
| Sample 3 | 1.00E+05 | 8   | 8          |   |
| Sample 3 | 1.00E+06 | 8   | 8          |   |
| Sample 3 | 1.00E+07 | 8   | 8          |   |
| Sample 4 | 1.00E+02 | 0   | 8          |   |
| Sample 4 | 1.00E+03 | 0   | 8          |   |
| Sample 4 | 1.00E+04 | 0   | 8          |   |
| Sample 4 | 1.00E+05 | 0   | 8          |   |
| Sample 4 | 1.00E+06 | 0   | 8          |   |
| Sample 4 | 1.00E+07 | 0   | 8          |   |

**Columns**
- ID: A unique ID for each sample you titrated
- Dilution: The dilution of virus that was pipetted on the well
- CPE: Number of Wells that showed CPE
- Replicates: Total number of infected wells for this dilution

**Settings**
In the notebook you can further set the volume in µl of virus dilution per well and the decimal separator for input (in case of pasted tab separated table or .csv upload)

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
