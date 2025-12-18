import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def import_libraries():
    import mols2grid
    from mols2grid import datafiles
    from rdkit.Chem import Descriptors
    import marimo as mo
    return Descriptors, datafiles, mo, mols2grid


@app.cell
def create_sliders(mo):
    solubility_range = mo.ui.range_slider(
        -10, 2, 0.5,
        debounce=True,show_value=True,full_width=True,
        label="Solubility",
    )

    molwt_range = mo.ui.range_slider(
        0, 600, 10,
        debounce=True, show_value=True, full_width=True,
        label="MolWt",
    )

    # logp_range = mo.ui.range_slider(
    #     -10, 10, 1,
    #     debounce=True, show_value=True, full_width=True,
    #     label="LogP",
    # )

    # num_h_donors_range = mo.ui.range_slider(
    #     0, 20, 1,
    #     debounce=True, show_value=True, full_width=True,
    #     label="NumHDonors",
    # )

    # num_h_acceptors_range = mo.ui.range_slider(
    #     0, 20, 1,
    #     debounce=True, show_value=True, full_width=True,
    #     label="NumHAcceptors",
    # )
    return molwt_range, solubility_range


@app.cell
def load_data_and_grid(Descriptors, datafiles, mols2grid):
    # NOTE:
    # This cell is intentionally kept independent from the sliders.
    # In marimo, cells are re-executed whenever any of their dependencies change.
    # Keeping grid creation here prevents MolGrid.from_sdf(...) from being
    # re-run on every slider update, which would reset the widget state.

    grid = mols2grid.MolGrid.from_sdf(datafiles.SOLUBILITY_SDF, size=(120, 100))
    df = grid.dataframe
    df["MolWt"] = df["mol"].apply(Descriptors.ExactMolWt)
    # df["LogP"] = df["mol"].apply(Descriptors.MolLogP)
    # df["NumHDonors"] = df["mol"].apply(Descriptors.NumHDonors)
    # df["NumHAcceptors"] = df["mol"].apply(Descriptors.NumHAcceptors)
    view = grid.display(n_items_per_page=12)
    return df, grid, view


@app.cell
def apply_filters_and_layout(df, grid, mo, molwt_range, solubility_range, view):
    filters = {
        "SOL": solubility_range.value,
        "MolWt": molwt_range.value,
        # "LogP": logp_range.value,
        # "NumHDonors": num_h_donors_range.value,
        # "NumHAcceptors": num_h_acceptors_range.value,
    }

    mask = True
    for col, (lo, hi) in filters.items():
        mask &= df[col].between(lo, hi)

    results = df.loc[mask]
    grid.filter_by_index(results.index)

    mo.vstack([
        solubility_range,
        molwt_range,
        # logp_range,
        # num_h_donors_range,
        # num_h_acceptors_range,
        view
    ])
    return


if __name__ == "__main__":
    app.run()
