// list.js
var listObj = new List('mols2grid', {
    valueNames: {{ value_names }},
    item: {{ item_repr }},
    page: {{ n_items_per_page }},
    pagination: {
        name: "pagination",
        item: '<li class="page-item"><a class="page page-link" href="#" onclick="event.preventDefault()"></a></li>',
        innerWindow: 1,
        outerWindow: 1,
    },
});
listObj.remove("mols2grid-id", "0");
listObj.add({{ data }});
// filter
if (window.parent.mols2grid_lists === undefined) {
    window.parent.mols2grid_lists = {};
}
window.parent.mols2grid_lists[{{ grid_id | tojson }}] = listObj;

{% if selection %}
// selection
{% include 'js/molstorage.js' %}
var SELECTION = new MolStorage();
{% endif %}

{% if selection or callback %}
// kernel
{% include 'js/kernel.js' %}
{% endif %}

{% if selection and cached_selection %}
// restore checkbox state
listObj.on("updated", function (list) {
    $('#mols2grid .cell input[checked="false"]').prop("checked", false);
});
{% endif %}

// sort
{% include 'js/sort.js' %}

{% if whole_cell_style %}
// add style for whole cell
listObj.on("updated", function (list) {
    $('#mols2grid div.cell').each(function() {
        var $t = $(this);
        $t.attr({style: $t.attr('data-cellstyle')})
          .removeAttr('data-cellstyle');
    });
});
{% endif %}

{% if tooltip %}
// tooltips
$.fn.tooltip.Constructor.Default.whiteList.span = ['style']
listObj.on("updated", function (list) {
    $(function () {
        // hide previous popovers
        $('#mols2grid a.page-link').click(function(e) {
            $('.mols2grid-tooltip[data-toggle="popover"]').popover('hide')
        });
        // create new popover
        $('.mols2grid-tooltip[data-toggle="popover"]').popover({
            placement: {{ tooltip_placement }},
            trigger: {{ tooltip_trigger }},
            html: true,
            sanitize: false,
        });
    })
});
{% endif %}

{% if selection %}
// selection modifyers and export options
{% include 'js/selection_actions.js' %}
{% endif %}

{% if callback %}
// callback
{% include 'js/callback.js' %}
{% endif %}

{% if onthefly %}
// generate images for the currently displayed molecules
var draw_opts = {{ json_draw_opts }};
var json_draw_opts = JSON.stringify(draw_opts);
var smarts_matches = {};
{% endif %}

// Load RDKit
window
    .initRDKitModule()
    .then(function(RDKit) {
        console.log('RDKit version: ', RDKit.version());
        window.RDKit = RDKit;

        {% if onthefly %}
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
            $('#mols2grid .cell').each(function() {
                var $t = $(this);
                var smiles = $t.children(".data-{{ smiles_col }}").first().text()
                var svg = draw_mol(smiles);
                $t.children(".data-img").html(svg);
            });
        });
        {% endif %}

        // search bar
        {% include 'js/search.js' %}

        listObj.update();
        // resize if iframe
        if (window.frameElement) {
            window.top.fit_height(window.frameElement);
        }

    })
    .catch(() => {
        console.log("Error loading RDKit");
    });