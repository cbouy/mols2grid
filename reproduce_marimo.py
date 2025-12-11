import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from rdkit import Chem
    import mols2grid
    import pandas as pd
    return Chem, mols2grid, pd


@app.cell
def _():
    smiles = """
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
def _(Chem, pd, smiles):
    sup = Chem.SmilesMolSupplierFromText(smiles , titleLine=False, delimiter="\t")
    mols = [mol for mol in sup if mol]

    df = pd.DataFrame(
        {
            "mol": mols,
            "name": [mol.GetProp('_Name') for mol in mols]
        }
    )
    df
    return (df,)


@app.cell
def _(df, mols2grid):
    mols2grid.display(df, mol_col="mol")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
