import marimo

__generated_with = "0.17.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # TCID50 calculation
    This marimo notebooks calculates the TCID50 from a CSV file of the following structure:

    ID | Dilution | CPE | Replicates 
    ---|----------|-----|----------
    1  |10        | 0   | 8
    1  |100       | 1   | 8
    1  |1000      | 6   | 8
    1  |10000     | 8   | 8
    1  |100000 | 8 | 8 
    2  |10|0|8
    ||...|

    where "ID" is a a unique value to identify the sample, "Dilution" is the dilution factor, "CPE" is the number of wells showing cyptopathic effect and "Replicates" is the total number of replicates for the sample and dilution.

    The script will calculate the logit for the CPE probablity (linear transformation of logistic curve), determine the x-axis intercept (corresponding to p(CPE)=0.5) and generate a plot for each sample showing the transformed dose-response curve ( logit(p) ~ log10(p(CPE)) ).
    For each value without a valid regression (zero or all wells with CPE) it will set the value to 20% below or above the lowest or highest dilution. This has to be taken into account when visualizing the data.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Settings
    Change the inputfile and volume variable and re-run all cells
    """
    )
    return


@app.cell
def _():
    inputfile = "20251104_NE45G_TCID50.xlsx"
    volume = 0.05 # mL
    return inputfile, volume


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Dose-Response curves for all samples
    Y-axis (p(CPE)) is [logit](https://de.wikipedia.org/wiki/Logit)-transformed (linearisation of sigmoidalen functions)
    """
    )
    return


@app.cell(hide_code=True)
def _(alt, df, mo):
    _points = (
        alt.Chart(df)
        .mark_point()
        .encode(
            x=alt.X("Dilution:Q").title("Dilution [log10]"), y=alt.Y("logit:Q")
        )
    )
    _reg = _points.transform_regression(
        "Dilution", "logit", groupby=["ID"]
    ).mark_line(size=2)
    _hline = (
        alt.Chart()
        .mark_rule(
            color="red",
        )
        .encode(y=alt.datum(0))
    )
    _chart = (
        (_points + _reg + _hline)
        .properties(width=100, height=100)
        .facet(column="ID")
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## TCID50 for all samples
    The following plot shows the log10(TCID50) for all samples. Points that fall outside the limit of detection are colored gray.
    """
    )
    return


@app.cell(hide_code=True)
def _(alt, linreg, mo, np):
    _linreg = linreg.copy()
    _linreg["TCID50/mL"] = np.log10(_linreg["TCID50/mL"])
    _chart = (
        alt.Chart(_linreg.reset_index())
        .mark_point()
        .encode(
            x="ID",
            y=alt.Y("TCID50/mL").title("TCID50/mL (log10)").scale(zero=False),
            color=alt.when(alt.datum.outside_detection == True).then(
                alt.value("lightgray")
            ),
        )
        .properties(height=100, width=20 + 12 * len(_linreg))
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _():


    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import pandas as pd
    import altair as alt
    import numpy as np
    from scipy import stats
    from sklearn.linear_model import LogisticRegression
    import os
    from pathlib import Path
    import tomllib
    return alt, mo, np, pd, stats


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(file_input, inputfile, np, pd, stats, volume):
    df = pd.read_excel(inputfile)
    df = df.set_index("ID")
    output = "".join(file_input.value.split(".")[:-1]) + "_out.csv"

    # apply continuity correction to number of wells with CPE
    cpe_cont = (df["CPE"] + 0.5) / (df["Replicates"] + 1)

    # Calculate logit (for linear regression of logistic distribution)
    df["logit"] = np.log(cpe_cont / (1 - cpe_cont))

    # Apply volume to dilution
    df["Dilution"] = df["Dilution"] / float(volume)

    lower_limit = df.groupby("ID")["Dilution"].min()
    upper_limit = df.groupby("ID")["Dilution"].max()

    # log10 transform Dilution
    df["Dilution"] = np.log10(df["Dilution"])

    # Perform linear regression for each ID seperatly
    linreg = df.groupby("ID").apply(
        lambda x: pd.Series(stats.linregress(x["Dilution"], x["logit"])._asdict())[
            ["slope", "intercept", "intercept_stderr"]
        ],
        include_groups=False,
    )

    # Add lower and upper limit of detection
    linreg["Upper limit"] = linreg.index.map(upper_limit)
    linreg["Lower limit"] = linreg.index.map(lower_limit)
    linreg["TCID50/mL"] = 10 ** (-linreg["intercept"] / linreg["slope"])

    linreg.loc[(linreg["slope"] == 0), "outside_detection"] = True
    linreg.loc[linreg["outside_detection"].isna(), "outside_detection"] = False
    linreg.loc[(linreg["intercept"] < 0) & (linreg["slope"] == 0), "TCID50/mL"] = (
        linreg["Lower limit"] * 0.8
    )
    linreg.loc[(linreg["intercept"] > 0) & (linreg["slope"] == 0), "TCID50/mL"] = (
        linreg["Upper limit"] * 1.2
    )
    linreg.to_csv(
        output,
    )
    return df, linreg


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
