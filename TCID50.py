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
    from pathlib import Path
    return Path, alt, io, locale, mo, np, pd, sm


@app.cell
def _(Path, mo):
    readme = Path(mo.notebook_location() / "public" / "README.md")
    if not readme.is_file():
        readme = Path(mo.notebook_location() /  "README.md") 
    mo.stop(not readme.is_file(), "Did not find README.md")
    with open(str(readme), "r") as f:
        _readme = f.read()
    mo.accordion(
        {"Click here for instructions on how to use the notebook": mo.md(_readme)}
    )
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
            <div style="grid-column: span 2 / span 2; ">
            <b>Paste in tab separated data</b>
            {text}
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
                io.StringIO(form.value["text"]), decimal=form.value["dec"]
            )


    # def validate_dataframe(df):
    #    if not all(["Dilution", "CPE", "Rep"] in df.columns):
    #        mo.callout(
    #            "Not all of the following columns exist in dataframe: Dilution, #CPE, Rep",
    #            kind="danger",
    #        )


    input_df = read_input(form)
    order = input_df["ID"].unique()
    # validate_dataframe(input_df)
    input_df = input_df.melt(
        id_vars="ID",
    )

    input_df[["Dilution", "variable"]] = input_df["variable"].str.split(
        " ", expand=True
    )
    input_df = input_df.pivot_table(
        index=["ID", "Dilution"], columns="variable", values="value"
    ).reset_index()
    input_df = input_df.dropna()
    input_df["Dilution"] = input_df["Dilution"].astype(float)
    input_df["Dilution"] = input_df["Dilution"] * 1000 / form.value["volumen"]
    input_df["Dilution"] = np.log10(input_df["Dilution"])
    input_df = input_df[input_df["CPE"] != ""]
    input_df["CPE"] = input_df["CPE"].astype(int)
    input_df["Rep"] = input_df["Rep"].astype(int)
    input_df["Fraction"] = input_df["CPE"] / input_df["Rep"]
    input_df
    return input_df, order


@app.cell
def _(mo):
    mo.md(r"""
    **Results**
    """)
    return


@app.cell
def _(input_df, mo, np, order, pd, sm):
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
        if all(df["CPE"] == df["Rep"]):
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
        y = df["CPE"] / df["Rep"]
        model = sm.GLM(y, X, family=sm.families.Binomial(), freq_weights=df["Rep"])
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
def _(np, output_df, pd, sm):
    def predict(result):
        xmin = result.model.data.orig_exog["Dilution"].min()
        xmax = result.model.data.orig_exog["Dilution"].max()
        x = np.linspace(xmin, xmax, 200)
        X = sm.add_constant(x)
        y = result.predict(X)
        return pd.Series({"Dilution": x, "Fraction": y})


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
    (_point + _line).properties(width=100, height=100).facet(
        alt.Facet("ID").sort(order), columns=5
    )
    return


if __name__ == "__main__":
    app.run()
