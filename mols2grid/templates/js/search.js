function SmartsSearch(query, columns) {
    var smiles_col = columns[0];
    smarts_matches = {};
    var query = $('#mols2grid #searchbar').val();
    var qmol = RDKit.get_qmol(query);
    if (qmol.is_valid()) {
        for (var k = 0, kl = listObj.items.length; k < kl; k++) {
            var item = listObj.items[k];
            var smiles = item.values()[smiles_col]
            var mol = RDKit.get_mol(smiles);
            if (mol.is_valid()) {
                var results = JSON.parse(mol.get_substruct_match(qmol));
                if (results.atoms) {
                    item.found = true;
                    {% if onthefly and substruct_highlight %}
                    smarts_matches[smiles] = results;
                    {% endif %}
                } else {
                    item.found = false;
                }
            } else {
                item.found = false;
            }
            mol.delete();
        }
    }
    qmol.delete();
}
var search_type = "Text";
$('#mols2grid .search-btn').click(function() {
    search_type = $(this).text();
    $('#mols2grid button.search-btn.active').removeClass("active");
    $(this).addClass("active");
});
$('#mols2grid #searchbar').on("keyup", function(e) {
    var query = e.target.value;
    if (search_type === "Text") {
        smarts_matches = {};
        listObj.search(query, {{ search_cols }});
    } else {
        listObj.search(query, ["data-{{ smiles_col }}"], SmartsSearch);
    }
});