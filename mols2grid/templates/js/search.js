function SmartsSearch(query, columns) {
    var smiles_col = columns[0];
    smarts_matches = {};
    var query = $('#mols2grid #searchbar').val();
    var qmol = RDKit.get_qmol(query);
    if (qmol.is_valid()) {
        listObj.items.forEach(function (item) {
            var smiles = item.values()[smiles_col]
            var mol = RDKit.get_mol(smiles);
            if (mol.is_valid()) {
                var results = mol.get_substruct_matches(qmol);
                if (results === "\{\}") {
                    item.found = false;
                } else {
                    item.found = true;
                    {% if onthefly and substruct_highlight %}
                    results = JSON.parse(results);
                    {% if single_highlight %}
                    var highlights = results[0]
                    {% else %}
                    var highlights = {"atoms": [], "bonds": []};
                    results.forEach(function (match) {
                        highlights["atoms"].push(...match.atoms)
                        highlights["bonds"].push(...match.bonds)
                    });
                    {% endif %}
                    var index = item.values()["mols2grid-id"];
                    smarts_matches[index] = highlights;
                    {% endif %}
                }
            } else {
                item.found = false;
            }
            mol.delete();
        });
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