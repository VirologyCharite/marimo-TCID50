import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import locale
    import io
    import traceback
    import statsmodels.api as sm
    import numpy as np
    import altair as alt
    return alt, io, locale, mo, np, pd, sm


@app.cell
def _(locale, mo):
    current_locale = locale.setlocale(locale.LC_NUMERIC)
    locale.setlocale(locale.LC_NUMERIC, "")

    decimal_separator = locale.localeconv()["decimal_point"]

    locale.setlocale(locale.LC_NUMERIC, current_locale)

    form = (
        mo.md("""
        <style>
            .form_container {{
              display: grid;
              grid-template-columns: auto auto;
              background-color: dodgerblue;
            }}
            .form_container div {{
              background-color: #f1f1f1;
              padding: 10px;
            }}
        </style>
        <div class="form_container">
            <div>
            <b>Paste in tab separated data</b>
            {text}
            </div>
            <div>
            <b>Or upload excel/csv file</b>
            {file}
            </div>
            <div style="grid-column: span 2 / span 2; ">
            <b>Settings</b><br>
            {dec}
            <br>
            {volumen}
            </div>
        </div>

    """)
        .batch(
            text=mo.ui.text_area(full_width=True),
            file=mo.ui.file(kind="area"),
            dec=mo.ui.dropdown(
                options=[".", ","],
                value=decimal_separator,
                label="Decimal separator: ",
            ),
            volumen=mo.ui.number(value=10, start=1, label="Volume/Well [µL]:"),
        )
        .form(show_clear_button=True, bordered=False)
    )
    form
    return (form,)


@app.cell
def _(form, io, mo, np, pd):
    def read_input(form):
        mo.stop(
            not form.value,
            mo.callout(
                "Neither input file nor tab separated text provided.",
                kind="danger",
            ),
        )
        mo.stop(
            form.value["text"] != "" and len(form.value["file"]) > 0,
            mo.callout(
                "Both input file and tab separated text provided.", kind="danger"
            ),
        )
        mo.stop(
            form.value["text"] == "" and len(form.value["file"]) == 0,
            mo.callout(
                "Neither input file nor tab separated text provided.",
                kind="danger",
            ),
        )
        if form.value["text"] != "":
            return pd.read_table(
                io.StringIO(form.value["text"]), decimal=form.value["dec"]
            )
        else:
            if form.value["file"][0].name.endswith(".csv"):
                input = pd.read_csv(
                    io.BytesIO(form.value["file"][0].contents),
                    decimal=form.value["dec"],
                )
                return input
            if form.value["file"][0].name.endswith(".tsv"):
                input = pd.read_table(
                    io.BytesIO(form.value["file"][0].contents),
                    decimal=form.value["dec"],
                )
                return input
            if form.value["file"][0].name.endswith(".xlsx"):
                input = pd.read_excel(
                    io.BytesIO(form.value["file"][0].contents),
                )
                return input
            raise ValueError(
                "File does not have a supported file extension (.csv, .tsv, .xlsx)"
            )


    input_df = read_input(form)
    input_df["Dilution"] = input_df["Dilution"] * 1000 / form.value["volumen"]
    input_df["Dilution"] = np.log10(input_df["Dilution"])
    input_df["Fraction"] = input_df["CPE"] / input_df["Replicates"]
    return (input_df,)


@app.cell
def _(mo):
    mo.md(r"""
    **Results**
    """)
    return


@app.cell
def _(input_df, mo, np, pd, sm):
    def calculate_tcid50(
        df,
    ):
        if all(df["CPE"] == 0):
            return pd.Series(
                {
                    "log_TCID50_mL": None,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_high": df["Dilution"].max(),
                    "result": None,
                    "message": "below detection limit",
                }
            )
        if all(df["CPE"] == df["Replicates"]):
            return pd.Series(
                {
                    "log_TCID50_mL": None,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_high": df["Dilution"].max(),
                    "result": None,
                    "message": "above detection limit",
                }
            )
        X = sm.add_constant(df["Dilution"])
        y = df["CPE"] / df["Replicates"]
        model = sm.GLM(
            y, X, family=sm.families.Binomial(), freq_weights=df["Replicates"]
        )
        results = model.fit()
        beta_0, beta_1 = results.params
        tcid50 = -beta_0 / beta_1
        if tcid50 < df["Dilution"].min() + np.log10(0.99):
            return pd.Series(
                {
                    "log_TCID50_mL": tcid50,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_high": df["Dilution"].max(),
                    "result": results,
                    "message": "below detection limit",
                }
            )
        if tcid50 > df["Dilution"].max() + np.log10(1.01):
            return pd.Series(
                {
                    "log_TCID50_mL": tcid50,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_high": df["Dilution"].max(),
                    "result": results,
                    "message": "above detection limit",
                }
            )
        return pd.Series(
            {
                "log_TCID50_mL": tcid50,
                "detection_limit_low": df["Dilution"].min(),
                "detection_limit_high": df["Dilution"].max(),
                "result": results,
                "message": None,
            },
        )


    output_df = input_df.groupby("ID").apply(
        lambda x: calculate_tcid50(
            x,
        ),
        include_groups=False,
    )
    output_df["log_PFU_mL"] = output_df["log_TCID50_mL"] + np.log10(0.7)
    _table = mo.ui.table(
        data=output_df.drop("result", axis=1),
    )
    mo.output.replace(_table)
    return (output_df,)


@app.cell
def _(np, output_df, pd, sm):
    def predict(result):
        xmin = result.model.data.orig_exog["Dilution"].min()
        xmax = result.model.data.orig_exog["Dilution"].max()
        x = np.linspace(xmin, xmax, 200)
        X = sm.add_constant(x)
        y = result.predict(X)
        return pd.Series({"Dilution": x, "Fraction": y})


    predicted = (
        output_df.dropna(subset=["result"])["result"]
        .apply(lambda x: predict(x))
        .explode(column=["Dilution", "Fraction"])
        .reset_index()
    )
    return (predicted,)


@app.cell
def _(alt, input_df, pd, predicted):
    _concat_data = pd.concat(
        [
            predicted.assign(data_source="predicted"),
            input_df.assign(data_source="observed"),
        ]
    )
    _line = (
        alt.Chart(_concat_data)
        .mark_line()
        .encode(
            x=alt.X("Dilution:Q").scale(domainMin=_concat_data["Dilution"].min()),
            y=alt.Y("Fraction:Q"),
        )
        .transform_filter(alt.datum.data_source == "predicted")
    )
    _point = (
        alt.Chart(_concat_data)
        .mark_point()
        .encode(
            x=alt.X("Dilution").scale(domainMin=_concat_data["Dilution"].min()),
            y=alt.Y("Fraction"),
        )
        .transform_filter(alt.datum.data_source == "observed")
    )
    (_point + _line).properties(width=100, height=100).facet("ID", columns=5)
    return


if __name__ == "__main__":
    app.run()
