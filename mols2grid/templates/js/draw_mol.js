// generate images for the currently displayed molecules
RDKit.prefer_coordgen({{ prefer_coordGen | tojson }});
function draw_mol(smiles) {
    var mol = RDKit.get_mol(smiles);
    var svg = "";
    if (mol.is_valid()) {
        var highlights = smarts_matches[smiles];
        if (highlights) {
            var details = Object.assign({}, draw_opts, highlights);
            details = JSON.stringify(details);
        } else {
            var details = json_draw_opts;
        }
        svg = mol.get_svg_with_highlights(details);
    }    
    mol.delete();
    if (svg == "") {
        return "<img/>";
    }
    return svg;
}
listObj.on("updated", function (list) {
    var query = $('#mols2grid #searchbar').val();
    if (query === "") {
        smarts_matches = {};
    }
    $('#mols2grid .cell').each(function() {
        var $t = $(this);
        var smiles = $t.children(".data-{{ smiles_col }}").first().text()
        var svg = draw_mol(smiles);
        $t.children(".data-img").html(svg);
    });
});