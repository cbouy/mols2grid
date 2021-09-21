// RDKit
onRuntimeInitialized: initRDKitModule().then(function(instance) {
    RDKitModule = instance;
    console.log('RDKit version: ' + RDKitModule.version());
});
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
{% if selection %}
if (window.parent.mols2grid_iframe === undefined) {
    window.parent.mols2grid_iframe = {};
}
// selection
{% include 'js/molstorage.js' %}
var SELECTION = new MolStorage();
{% include 'js/kernel.js' %}
{% endif %}
// sort
{% include 'js/sort.js' %}
// search bar
{% include 'js/search.js' %}
{% if whole_cell_style %}
// add style for whole cell
listObj.on("updated", function (list) {
    $('#mols2grid div.cell').each(function() {
        var $t = $(this);
        $t.attr({
            style: $t.attr('data-cellstyle')
        })
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
{% include 'js/selection_actions.js' %}
{% endif %}
{%if tooltip or selection or whole_cell_style %}
listObj.update();
{% endif %}
window.parent.mols2grid_iframe["{{ grid_id }}"] = window;
