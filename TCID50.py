# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "altair==6.2.2",
#     "marimo>=0.24.2",
#     "numpy==2.5.3",
#     "pandas==3.0.5",
#     "pytest==9.1.1",
#     "statsmodels==0.15.0",
# ]
# ///

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import locale
    import io
    import traceback
    import statsmodels.api as sm
    import pytest
    import numpy as np
    import altair as alt
    from pathlib import Path

    return alt, io, locale, mo, np, pd, pytest, sm


@app.cell
def _(mo):
    mo.md("""
    Click here for usage instructions: https://github.com/VirologyCharite/marimo-TCID50
    """)
    return


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
            .form_container marimo-text-area {{
            font-size:xx-large !important;
            background-color:red;
            }}
        </style>
        <div class="form_container">
            <div style="grid-column: span 1 / span 2; ">
            <b>Paste in tab separated data</b>
            {text}
            </div>

            <div style="grid-column: span 1 / span 2; ">
            <b>Settings</b><br>
            <br>
            {volumen}
            </div>
        </div>

    """)
        .batch(
            text=mo.ui.text_area(full_width=True),
            file=mo.ui.file(kind="area"),
            volumen=mo.ui.number(value=10, start=1, label="Volume/Well [µL]:"),
        )
        .form(show_clear_button=True, bordered=False)
    )
    form
    return (form,)


@app.cell
def _(form, io, mo, np, pd):
    def read_input(form):
        """validate the form input and return the tab separated text or the uploaded file as a pandas DataFrame"""
        mo.stop(
            not form.value,
            mo.callout(
                "Neither input file nor tab separated text provided.",
                kind="danger",
            ),
        )
        mo.stop(
            form.value["text"] == "",
            mo.callout(
                "Neither input file nor tab separated text provided.",
                kind="danger",
            ),
        )
        if form.value["text"] != "":
            return pd.read_table(
                io.StringIO(form.value["text"]), 
            )

    # def validate_dataframe(df):
    #    if not all(["Dilution", "CPE", "Rep"] in df.columns):
    #        mo.callout(
    #            "Not all of the following columns exist in dataframe: Dilution, #CPE, Rep",
    #            kind="danger",
    #        )

    input_df = read_input(form)
    input_df["Dilution"] = input_df["Dilution"]*1000/form.value["volumen"]
    input_df["Dilution"] = np.log10(input_df["Dilution"])
    order = input_df["ID"].unique()
    # validate_dataframe(input_df)
    input_df = input_df.dropna()
    input_df["Fraction"] = input_df["CPE"] / input_df["Total"]
    input_df
    return input_df, order


@app.cell
def _(mo):
    mo.md(r"""
    **Results**
    """)
    return


@app.cell
def calculate_tcid50(np, pd, sm):
    def calculate_tcid50(
        df,
    ):
        if all(df["CPE"] == 0):
            return pd.Series(
                {
                    "log_TCID50_mL": None,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_up": df["Dilution"].max(),
                    "result": None,
                    "message": "below detection limit",
                }
            )
        if all(df["CPE"] == df["Total"]):
            return pd.Series(
                {
                    "log_TCID50_mL": None,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_up": df["Dilution"].max(),
                    "result": None,
                    "message": "above detection limit",
                }
            )
        X = sm.add_constant(df["Dilution"])
        y = df["CPE"] / df["Total"]
        model = sm.GLM(
            y, X, family=sm.families.Binomial(), freq_weights=df["Total"]
        )
        results = model.fit()
        beta_0, beta_1 = results.params
        tcid50 = -beta_0 / beta_1
        if tcid50 < df["Dilution"].min() + np.log10(0.99):
            return pd.Series(
                {
                    "log_TCID50_mL": tcid50,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_up": df["Dilution"].max(),
                    "result": results,
                    "message": "below detection limit",
                }
            )
        if tcid50 > df["Dilution"].max() + np.log10(1.01):
            return pd.Series(
                {
                    "log_TCID50_mL": tcid50,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_up": df["Dilution"].max(),
                    "result": results,
                    "message": "above detection limit",
                }
            )
        return pd.Series(
            {
                "log_TCID50_mL": tcid50,
                "detection_limit_low": df["Dilution"].min(),
                "detection_limit_up": df["Dilution"].max(),
                "result": results,
                "message": None,
            },
        )

    return (calculate_tcid50,)


@app.cell
def _(calculate_tcid50, input_df, mo, np, order, pd):
    output_df = (
        input_df.groupby("ID")
        .apply(
            lambda x: calculate_tcid50(
                x,
            ),
            include_groups=False,
        )
        .reset_index()
    )
    output_df["ID"] = pd.Categorical(output_df["ID"], order)
    output_df = output_df.sort_values("ID")
    output_df["log_PFU_mL"] = output_df["log_TCID50_mL"] + np.log10(np.log(2))
    output_df["PFU_mL"] = 10 ** output_df["log_PFU_mL"]
    output_df["TCID50_mL"] = 10 ** output_df["log_TCID50_mL"]
    output_df = output_df[
        [
            "ID",
            "detection_limit_low",
            "detection_limit_up",
            "log_TCID50_mL",
            "log_PFU_mL",
            "TCID50_mL",
            "PFU_mL",
            "message",
            "result",
        ]
    ]
    _table = mo.ui.table(
        data=output_df.drop("result", axis=1),
        format_mapping={"TCID50_mL": "{:.3e}", "PFU_mL": "{:.3e}"},
    )
    mo.output.replace(_table)
    return (output_df,)


@app.cell
def predict(np, pd, sm):
    def predict(result):
        xmin = result.model.data.orig_exog["Dilution"].min()
        xmax = result.model.data.orig_exog["Dilution"].max()
        x = np.linspace(xmin, xmax, 200)
        X = sm.add_constant(x)
        y = result.predict(X)
        return pd.Series({"Dilution": x, "Fraction": y})

    return (predict,)


@app.cell
def _(output_df, predict):
    predicted = (
        output_df.set_index("ID")
        .dropna(subset=["result"])["result"]
        .apply(lambda x: predict(x))
        .explode(column=["Dilution", "Fraction"])
        .reset_index()
    )
    return (predicted,)


@app.cell
def _():
    return


@app.cell
def _(alt, input_df, order, pd, predicted):
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
            x=alt.X("Dilution:Q").scale(
                domainMin=_concat_data["Dilution"].min()
            ),
            y=alt.Y("Fraction:Q"),
        )
        .transform_filter(alt.datum.data_source == "predicted")
    )
    _point = (
        alt.Chart(_concat_data)
        .mark_point()
        .encode(
            x=alt.X("Dilution").scale(
                domainMin=_concat_data["Dilution"].min()
            ),
            y=alt.Y("Fraction"),
        )
        .transform_filter(alt.datum.data_source == "observed")
    )
    (_point + _line).properties(width=100, height=100).facet(
        alt.Facet("ID").sort(order), columns=5
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    **Tests**
    """)
    return


@app.cell
def test_calculate_tcid50_below_detection_limit(calculate_tcid50, pd):
    def test_calculate_tcid50_below_detection_limit():
        df = pd.DataFrame(
            {
                "Dilution": [2, 3, 4, 5, 6],
                "CPE": [0, 0, 0, 0, 0],
                "Total": [8, 8, 8, 8, 8],
            }
        )
        result = calculate_tcid50(df)
        assert result["message"] == "below detection limit"
        assert result["log_TCID50_mL"] is None
        assert result["result"] is None
        assert result["detection_limit_low"] == 2
        assert result["detection_limit_up"] == 6

    return


@app.cell
def test_calculate_tcid50_above_detection_limit(calculate_tcid50, pd):
    def test_calculate_tcid50_above_detection_limit():
        df = pd.DataFrame(
            {
                "Dilution": [2, 3, 4, 5, 6],
                "CPE": [8, 8, 8, 8, 8],
                "Total": [8, 8, 8, 8, 8],
            }
        )
        result = calculate_tcid50(df)
        assert result["message"] == "above detection limit"
        assert result["log_TCID50_mL"] is None
        assert result["result"] is None

    return


@app.cell
def test_calculate_tcid50_typical_fit(calculate_tcid50, pd, pytest):
    def test_calculate_tcid50_typical_fit():
        df = pd.DataFrame(
            {
                "Dilution": [2, 3, 4, 5, 6],
                "CPE": [8, 8, 6, 2, 0],
                "Total": [8, 8, 8, 8, 8],
            }
        )
        result = calculate_tcid50(df)
        assert result["message"] is None
        assert result["result"] is not None
        assert result["log_TCID50_mL"] == pytest.approx(4.501, abs=0.01)
        assert (
            df["Dilution"].min()
            < result["log_TCID50_mL"]
            < df["Dilution"].max()
        )

    return


@app.cell
def test_predict_matches_fit_range(calculate_tcid50, pd, predict, pytest):
    def test_predict_matches_fit_range():
        df = pd.DataFrame(
            {
                "Dilution": [2, 3, 4, 5, 6],
                "CPE": [8, 8, 6, 2, 0],
                "Total": [8, 8, 8, 8, 8],
            }
        )
        fit = calculate_tcid50(df)
        curve = predict(fit["result"])
        assert list(curve.index) == ["Dilution", "Fraction"]
        assert len(curve["Dilution"]) == 200
        assert len(curve["Fraction"]) == 200
        assert curve["Dilution"].min() == pytest.approx(df["Dilution"].min())
        assert curve["Dilution"].max() == pytest.approx(df["Dilution"].max())
        assert ((curve["Fraction"] >= 0) & (curve["Fraction"] <= 1)).all()

    return


if __name__ == "__main__":
    app.run()
