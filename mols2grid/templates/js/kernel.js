if (window.parent.IPython !== undefined) {
    // Jupyter notebook
    var kernel = window.parent.IPython.notebook.kernel;
    kernel.execute('from mols2grid.select import register as _m2g_reg')
    function add_selection(grid_id, _id, smiles) {
        SELECTION.set(_id, smiles);
        kernel.execute("_m2g_reg.add_selection('"+grid_id+"', "+_id+",'"+smiles+"')")
    }
    function del_selection(grid_id, _id) {
        SELECTION.delete(_id);
        kernel.execute("_m2g_reg.del_selection('"+grid_id+"', "+_id+")");
    }
} else if (window.parent.google !== undefined) {
    // Google colab
    var kernel = window.parent.google.colab.kernel;
    function add_selection(grid_id, _id, smiles) {
        SELECTION.set(_id, smiles);
        (async function() {
        const result = await kernel.invokeFunction('_m2g_reg.add_selection',
                                                   [grid_id, _id, smiles], {});
        })();
    }
    function del_selection(grid_id, _id) {
        SELECTION.delete(_id);
        (async function() {
        const result = await kernel.invokeFunction('_m2g_reg.del_selection',
                                                   [grid_id, _id], {});
        })();
    }
} else {
    function add_selection(grid_id, _id, smiles) {
        SELECTION.set(_id, smiles);
    }
    function del_selection(grid_id, _id) {
        SELECTION.delete(_id);
    }
}