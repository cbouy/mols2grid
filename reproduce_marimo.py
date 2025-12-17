import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    from io import StringIO

    import pandas as pd

    import mols2grid

    return StringIO, mols2grid, pd


@app.cell
def _():
    smiles = """SMILES	NAME
    CCO	Ethanol
    CC(=O)O	Acetic_acid
    CCC	Propane
    C1=CC=CC=C1	Benzene
    COC	Dimethyl_ether
    CCN	Ethylamine
    C(CO)O	Ethylene_glycol
    CCOCC	Diethyl_ether
    OC=O	Formic_acid
    CCOC(=O)C	Ethyl_acetate
    """
    return (smiles,)


@app.cell
def _(StringIO, pd, smiles):
    df = pd.read_csv(StringIO(smiles), delimiter="\t")
    df  # noqa: B018
    return (df,)


@app.cell
def _(df, mols2grid):
    mols2grid.display(df)


@app.cell
def _(mols2grid):
    mols2grid.get_selection()


if __name__ == "__main__":
    app.run()
