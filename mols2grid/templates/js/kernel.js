function add_selection(grid_id, _id, smiles) {
    SELECTION.multi_set(_id, smiles);
    let model = window.parent["_MOLS2GRID_" + grid_id];
    if (model) {
        model.set("selection", SELECTION.to_dict());
        model.save_changes();
    }
}
function del_selection(grid_id, _id) {
    SELECTION.multi_del(_id);
    let model = window.parent["_MOLS2GRID_" + grid_id];
    if (model) {
        model.set("selection", SELECTION.to_dict());
        model.save_changes();
    }
}
if (window.parent.IPython !== undefined) {
    // Jupyter notebook
    var kernel_env = "jupyter";
} else if (window.parent.google !== undefined) {
    // Google colab
    var kernel_env = "colab";
} else {
    var kernel_env = null;
}