var kernel_env = null;
if (window.parent.IPython !== undefined) {
    // Jupyter notebook
    kernel_env = "jupyter";
    var kernel = window.parent.IPython.notebook.kernel;
    kernel.execute('from mols2grid.select import register as _m2g_reg')
    function add_selection(grid_id, _id, smiles) {
        SELECTION.multi_set(_id, smiles);
        kernel.execute("_m2g_reg.add_selection('"+grid_id+"', "+JSON.stringify(_id)+","+JSON.stringify(smiles)+")");
    }
    function del_selection(grid_id, _id) {
        SELECTION.multi_del(_id);
        kernel.execute("_m2g_reg.del_selection('"+grid_id+"', "+JSON.stringify(_id)+")");
    }
} else if (window.parent.google !== undefined) {
    // Google colab
    kernel_env = "colab";
    var kernel = window.parent.google.colab.kernel;
    function add_selection(grid_id, _id, smiles) {
        SELECTION.multi_set(_id, smiles);
        (async function() {
        const result = await kernel.invokeFunction('_m2g_reg.add_selection',
                                                   [grid_id, _id, smiles], {});
        })();
    }
    function del_selection(grid_id, _id) {
        SELECTION.multi_del(_id);
        (async function() {
        const result = await kernel.invokeFunction('_m2g_reg.del_selection',
                                                   [grid_id, _id], {});
        })();
    }
} else {
    function add_selection(grid_id, _id, smiles) {
        SELECTION.multi_set(_id, smiles);
    }
    function del_selection(grid_id, _id) {
        SELECTION.multi_del(_id);
    }
}