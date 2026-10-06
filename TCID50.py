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

__generated_with = "0.25.1"
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

    return alt, io, mo, np, pd, pytest, sm


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Data input
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    Click here for usage instructions: https://github.com/VirologyCharite/marimo-TCID50
    """)
    return


@app.cell
def _(mo):
    sample_data = """ID	1	10	100	1000
    sample_1	6	6	3	0
    sample_2	4	1	0	0
    sample_3	6	6	6	4
    """
    form = (
        mo.md("""
        <style>

            .form_container div {{
              background-color: #f1f1f1;
              padding: 10px;
            }}

        </style>
        <div class="form_container">
            <div style="grid-column: span 2 / span 2; ">
            <b>Paste in tab separated data</b>
            {text}
            </div>

            <div style="grid-column: span 1 / span 2; ">
            <b>Settings</b><br>
            {volumen}
            <br>
            {replicates}
            </div>
        </div>

    """)
        .batch(
            text=mo.ui.text_area(full_width=True, value=sample_data),
            volumen=mo.ui.number(value=10, start=1, label="Volume/Well [µL]:"),
            replicates=mo.ui.number(
                value=6, start=1, label="Replicates/Dilution"
            ),
        )
        .form(show_clear_button=True, bordered=False)
    )
    form
    return (form,)


@app.cell
def _(form, io, mo, pd):
    mo.stop(
        (not form.value or form.value["text"].strip() == ""),
        mo.callout(
            "No input provided",
            kind="danger",
        ),
    )
    # pandas silently renames duplicated headers (10, 10 -> 10, 10.1), which would be read as a new dilution
    _header = form.value["text"].strip().splitlines()[0].split("\t")
    _duplicated = sorted({_h for _h in _header if _header.count(_h) > 1})
    mo.stop(
        len(_duplicated) > 0,
        mo.callout(
            f"Duplicated column headers: {', '.join(_duplicated)}",
            kind="danger",
        ),
    )
    _pasted = pd.read_table(io.StringIO(form.value["text"]))
    mo.stop(
        "ID" not in _pasted.columns,
        mo.callout("No ID column in provided data", kind="danger"),
    )
    editor = mo.ui.data_editor(
        _pasted, label="**Check and edit the submitted data**"
    )
    editor
    return (editor,)


@app.cell
def _(editor, form, mo, np, pd):
    def read_input(wide_df, replicates):
        """validate the wide table (one column per dilution) and convert it into a long pandas DataFrame"""
        df = pd.DataFrame(wide_df)
        df.columns = df.columns.astype(str)
        dilution_cols = []
        for col in df.columns:
            try:
                float(col)
                dilution_cols.append(col)
            except:
                continue
        print(f"Detected numeric columns: {dilution_cols}")
        mo.stop(
            len(dilution_cols) < 2,
            mo.callout(
                "At least two dilution columns are needed. The column headers must be numbers (1, 10, 100...) with a point as decimal separator",
                kind="danger",
            ),
        )
        mo.stop(
            any(float(col) <= 0 for col in dilution_cols),
            mo.callout(
                "Dilutions must be greater than 0. Provide the dilution as linear not log value (1, 10, 100... instead of 0, 1, 2)",
                kind="danger",
            ),
        )
        # rows added in the table are empty until they get an ID
        df = df[df["ID"].notna() & (df["ID"].astype(str).str.strip() != "")]
        df = df.melt(
            id_vars="ID",
            value_vars=dilution_cols,
            value_name="CPE",
            var_name="Dilution",
        )
        df["Dilution"] = df["Dilution"].astype(float)
        # cells edited in the table can come back as text or empty
        df["CPE"] = pd.to_numeric(df["CPE"], errors="coerce")
        mo.stop(
            df["CPE"].isna().all(),
            mo.callout("No CPE counts provided", kind="danger"),
        )
        invalid = df[(df["CPE"] < 0) | (df["CPE"] > replicates)]["ID"].unique()
        mo.stop(
            len(invalid) > 0,
            mo.callout(
                f"CPE counts must be between 0 and the replicates per dilution ({replicates}). Check: {', '.join(map(str, invalid))}",
                kind="danger",
            ),
        )
        df["Total"] = replicates
        return df

    input_df = read_input(editor.value, form.value["replicates"])
    _warnings = []
    # Check if the dilution is likely logarithmic
    if input_df["Dilution"].max() < 10:
        _warnings.append(
            "Your maximum dilution is < 10. Provide the dilution as linear not log value (1, 10, 100... instead of 0, 1, 2)"
        )
    input_df["Dilution"] = input_df["Dilution"] * 1000 / form.value["volumen"]
    input_df["Dilution"] = np.log10(input_df["Dilution"])
    order = input_df["ID"].unique()
    # validate_dataframe(input_df)
    input_df = input_df.dropna()
    input_df["Fraction"] = input_df["CPE"] / input_df["Total"]
    _pooled = input_df[input_df.duplicated(["ID", "Dilution"])]["ID"].unique()
    if len(_pooled) > 0:
        _warnings.append(
            f"Rows with the same ID are pooled into one fit: {', '.join(map(str, _pooled))}"
        )
    mo.vstack(
        [mo.callout(_w, kind="warn") for _w in _warnings]
    ) if _warnings else None
    return input_df, order


@app.cell
def _(input_df):
    input_df
    return


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
        if df["Dilution"].nunique() < 2 or df["CPE"].nunique() < 2:
            # a single dilution or the same CPE at every dilution cannot be fitted
            return pd.Series(
                {
                    "log_TCID50_mL": None,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_up": df["Dilution"].max(),
                    "result": None,
                    "message": "no fit possible: CPE does not change with dilution",
                }
            )
        X = sm.add_constant(df["Dilution"])
        y = df["CPE"] / df["Total"]
        model = sm.GLM(
            y, X, family=sm.families.Binomial(), freq_weights=df["Total"]
        )
        results = model.fit()
        beta_0, beta_1 = results.params
        if beta_1 >= 0:
            return pd.Series(
                {
                    "log_TCID50_mL": None,
                    "detection_limit_low": df["Dilution"].min(),
                    "detection_limit_up": df["Dilution"].max(),
                    "result": results,
                    "message": "no titre: CPE does not decrease with dilution",
                }
            )
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
def _(output_df, pd, predict):
    _fits = output_df.set_index("ID").dropna(subset=["result"])["result"]
    if len(_fits) > 0:
        predicted = (
            _fits.apply(lambda x: predict(x))
            .explode(column=["Dilution", "Fraction"])
            .reset_index()
        )
    else:
        # no sample could be fitted, e.g. all are below the detection limit
        predicted = pd.DataFrame(columns=["ID", "Dilution", "Fraction"])
    return (predicted,)


@app.cell(hide_code=True)
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


@app.cell
def test_calculate_tcid50_flat_response(calculate_tcid50, pd):
    def test_calculate_tcid50_flat_response():
        df = pd.DataFrame(
            {
                "Dilution": [2, 3, 4, 5, 6],
                "CPE": [4, 4, 4, 4, 4],
                "Total": [8, 8, 8, 8, 8],
            }
        )
        result = calculate_tcid50(df)
        assert (
            result["message"]
            == "no fit possible: CPE does not change with dilution"
        )
        assert result["log_TCID50_mL"] is None
        assert result["result"] is None

    return


@app.cell
def test_calculate_tcid50_inverted_response(calculate_tcid50, pd):
    def test_calculate_tcid50_inverted_response():
        df = pd.DataFrame(
            {
                "Dilution": [2, 3, 4, 5, 6],
                "CPE": [0, 2, 6, 8, 8],
                "Total": [8, 8, 8, 8, 8],
            }
        )
        result = calculate_tcid50(df)
        assert (
            result["message"]
            == "no titre: CPE does not decrease with dilution"
        )
        assert result["log_TCID50_mL"] is None

    return


if __name__ == "__main__":
    app.run()
