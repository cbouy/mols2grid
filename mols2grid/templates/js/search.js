function SmartsSearch(query, columns) {
    var smiles_col = columns[0];
    smarts_matches = {};
    var query = $('#mols2grid .m2g-searchbar').val();
    var qmol = RDKit.get_qmol(query);
    if (qmol.is_valid()) {
        listObj.items.forEach(function (item) {
            var smiles = item.values()[smiles_col]
            var mol = RDKit.get_mol(smiles, '{"removeHs": {{ removeHs | tojson }} }');
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
// Temporary fix for regex characters being escaped by list.js
// This extends String.replace to ignore the regex pattern used by list.js and returns
// the string unmodified. Other calls should not be affected, unless they use the exact
// same pattern and replacement value.
// TODO: remove once the issue is fixed in list.js and released
String.prototype.replace = (function(_super) {
    return function() {
        if (
            (arguments[0].toString() === '/[-[\\]{}()*+?.,\\\\^$|#]/g')
            && (arguments[1] === '\\$&')
        ) {
            if (this.length === 0) {
                return ''
            }
            return this
        }
        return _super.apply(this, arguments);
    };         
})(String.prototype.replace);

// Switch search type (Text or SMARTS)
$('#mols2grid .m2g-search-options .m2g-option').click(function() {
    search_type = $(this).text();
    $('#mols2grid .m2g-search-options .m2g-option.sel').removeClass("sel");
    $(this).addClass("sel");
});

// Searchbar update event handler
$('#mols2grid .m2g-searchbar').on("keyup", function(e) {
    var query = e.target.value;
    if (search_type === "Text") {
        smarts_matches = {};
        listObj.search(query, {{ search_cols }});
    } else {
        listObj.search(query, ["data-{{ smiles_col }}"], SmartsSearch);
    }
});