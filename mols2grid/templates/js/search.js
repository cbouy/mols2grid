function SmartsSearch(query, columns) {
    query = $('#mols2grid #searchbar').val();
    var smiles = columns[0];
    var qmol = RDKitModule.get_qmol(query);
    if (qmol.is_valid()) {
        for (var k = 0, kl = listObj.items.length; k < kl; k++) {
        var item = listObj.items[k];
        var mol = RDKitModule.get_mol(item.values()[smiles]);
        var results = JSON.parse(mol.get_substruct_match(qmol));
        if (results.atoms) {
            item.found = true;
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
        listObj.search(query, {{ search_cols }});
    } else {
        listObj.search(query, ["data-{{ smiles_col }}"], SmartsSearch);
    }
});