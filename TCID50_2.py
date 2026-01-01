import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import locale
    import io
    return io, locale, mo, pd


@app.cell
def _():
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
def _(form, io, mo, pd):
    def read_input(form):
        if form.value["text"] != "" and len(form.value["file"]) > 0:
            raise AttributeError(
                "Both input file and tab separated text provided."
            )
        if form.value["text"] == "" and len(form.value["file"]) == 0:
            raise AttributeError(
                "Neither input file nor tab separated text provided."
            )
        if form.value["text"] != "":
            input = pd.read_table(
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
            raise AttributeError(
                "File does not have a supported file extension (.csv, .tsv, .xlsx)"
            )
    try:
        tcid50_df = read_input(form)
        mo.output.replace(tcid50_df)
    except AttributeError as e:
        mo.output.replace(mo.callout(e, "danger"))
    return


if __name__ == "__main__":
    app.run()
