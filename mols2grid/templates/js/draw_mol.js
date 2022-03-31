// generate images for the currently displayed molecules
RDKit.prefer_coordgen({{ prefer_coordGen | tojson }});
function draw_mol(smiles, index, template_mol) {
    var mol = RDKit.get_mol(smiles, '{"removeHs": {{ removeHs | tojson }} }');
    var svg = "";
    if (mol.is_valid()) {
        var highlights = smarts_matches[index];
        if (highlights) {
            var details = Object.assign({}, draw_opts, highlights);
            details = JSON.stringify(details);
            mol.generate_aligned_coords(template_mol, {{ prefer_coordGen | tojson }});
        } else {
            var details = json_draw_opts;
        }
        svg = mol.get_svg_with_highlights(details);
    }
    mol.delete();
    if (svg == "") {
        return '<svg width="{{ cell_width }}" height="{{ height }}" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 {{ cell_width }} {{ height }}"></svg>';
    }
    return svg;
}
listObj.on("updated", function (list) {
    var query = $('#mols2grid #searchbar').val();
    var template_mol;
    if (query === "") {
        smarts_matches = {};
        template_mol = null;
    } else {
        template_mol = RDKit.get_qmol(query);
        template_mol.set_new_coords({{ prefer_coordGen | tojson }});
    }
    $('#mols2grid .cell').each(function() {
        var $t = $(this);
        var smiles = $t.children(".data-{{ smiles_col }}").first().text();
        var index = parseInt(this.getAttribute("data-mols2grid-id"));
        var svg = draw_mol(smiles, index, template_mol);
        $t.children(".data-img").html(svg);
    });
    if (template_mol) {
        template_mol.delete();
    }
});