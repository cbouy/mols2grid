if (window.parent.IPython !== undefined) {
    // Jupyter notebook
    var kernel = window.parent.IPython.notebook.kernel;
    kernel.execute('from mols2grid.select import register as m2g_sel')
    function update_current_grid(grid_id) {
        kernel.execute("m2g_sel._update_current_grid('"+grid_id+"')");
    }
    function set_selection(_id, smiles) {
        SELECTION.set(_id, smiles);
        kernel.execute("m2g_sel._set_selection("+_id+",'"+smiles+"')")
    }
    function del_selection(_id) {
        SELECTION.delete(_id);
        kernel.execute("m2g_sel._del_selection("+_id+")");
    }
} else if (window.parent.google !== undefined) {
    // Google colab
    var kernel = window.parent.google.colab.kernel;
    function update_current_grid(grid_id) {
        (async function() {
        const result = await kernel.invokeFunction('m2g_sel._update_current_grid',
                                                   [grid_id], {});
        })();
    }
    function set_selection(_id, smiles) {
        SELECTION.set(_id, smiles);
        (async function() {
        const result = await kernel.invokeFunction('m2g_sel._set_selection',
                                                   [_id, smiles], {});
        })();
    }
    function del_selection(_id) {
        SELECTION.delete(_id);
        (async function() {
        const result = await kernel.invokeFunction('m2g_sel._del_selection',
                                                   [_id], {});
        })();
    }
} else {
    function update_current_grid(grid_id) {
        null;
    }
    function set_selection(_id, smiles) {
        SELECTION.set(_id, smiles);
    }
    function del_selection(_id) {
        SELECTION.delete(_id);
    }
}