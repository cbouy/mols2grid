import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def import_libraries():
    import marimo as mo

    import mols2grid
    from mols2grid.datafiles import SOLUBILITY_SDF

    return SOLUBILITY_SDF, mo, mols2grid


@app.cell
def prepare(mo):
    from rdkit.Chem.Draw import rdMolDraw2D

    def mol_to_svg(mol, width=130, height=90, opts=None):
        if mol is None:
            return mo.Html("")

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

        if opts is not None:
            drawer.SetDrawOptions(opts)

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return mo.Html(drawer.GetDrawingText())

    solubility_range = mo.ui.range_slider(
        -10,
        2,
        0.5,
        debounce=True,
        show_value=True,
        full_width=True,
        label="Solubility",
    )
    return mol_to_svg, solubility_range


@app.cell
def create_grid(SOLUBILITY_SDF, mols2grid):
    # NOTE:
    # This cell is intentionally kept independent from the sliders.
    # In marimo, cells are re-executed whenever any of their dependencies change.
    # Keeping grid creation here prevents MolGrid.from_sdf(...) from being
    # re-run on every slider update, which would reset the widget state.
    grid = mols2grid.MolGrid.from_sdf(SOLUBILITY_SDF, size=(120, 100))
    get_selection_ids = grid.get_marimo_selection()
    view = grid.display(n_items_per_page=12, selection=True)
    return get_selection_ids, grid, view


@app.cell
def filter_and_display(grid, mo, solubility_range, view):
    mask = grid.dataframe["SOL"].between(*solubility_range.value)
    results = grid.dataframe.loc[mask]
    # Same as:
    # grid.dataframe["SOL"] >= solubility_range.value[0]) & \
    #   (grid.dataframe["SOL"] <= solubility_range.value[1])

    grid.filter_by_index(results.index)

    mo.vstack([solubility_range, view])
    return (results,)


@app.cell(hide_code=True)
def display_selection(get_selection_ids, mo, mol_to_svg, results):
    # This cell displays the selected molecules in a Marimo table.
    # The `mol_to_svg` function is used to render the molecule images
    # directly in the table.
    # We filter the dataframe based on the selection state (`get_selection_ids`)
    # from the grid above.

    selected = results[results["mols2grid-id"].isin(get_selection_ids())]

    # # If you want to use custom drawing options:
    # opts = rdMolDraw2D.MolDrawOptions()
    # opts.explicitMethyl = True

    table = mo.ui.table(
        selected.reset_index(drop=True).drop(columns="img"),
        format_mapping={
            "mol": mol_to_svg
            # "mol": lambda mol: mol_to_svg(mol, opts=opts)
            # If you want to use custom drawing options
        },
        freeze_columns_left=["mols2grid-id", "mol"],
        freeze_columns_right=["SOL"],
        label="Try selecting molecules from the grid above!!",
    )

    table  # noqa: B018


if __name__ == "__main__":
    app.run()
