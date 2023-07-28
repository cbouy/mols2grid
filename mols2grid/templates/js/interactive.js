// list.js
var listObj = new List('mols2grid', {
    listClass: 'm2g-list',
    valueNames: {{ value_names }},
    item: {{ item_repr }},
    page: {{ n_items_per_page }},
    pagination: {
        paginationClass: "m2g-pagination",
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
SELECTION.multi_set({{ cached_selection[0] }}, {{ cached_selection[1] }});
listObj.on("updated", function (list) {
    $('#mols2grid .m2g-cell input[checked="false"]').prop("checked", false);
});
{% endif %}

// sort
{% include 'js/sort.js' %}

{% if whole_cell_style %}
// add style for whole cell
listObj.on("updated", function (list) {
    $('#mols2grid div.m2g-cell').each(function() {
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
        // Hide previous popovers.
        $('#mols2grid a.page-link').click(function(e) {
            $('.m2g-tooltip[data-toggle="popover"]').popover('hide')
        });
        // Create new popover.
        $('.m2g-tooltip[data-toggle="popover"]').popover({
            placement: {{ tooltip_placement }},
            trigger: 'manual',
            html: true,
            sanitize: false,
        });
    })
});
{% endif %}

// grid interactions (select, click, tooltip, key events)
{% include 'js/grid_interaction.js' %}

{% if onthefly %}
// generate images for the currently displayed molecules
var draw_opts = {{ json_draw_opts }};
var json_draw_opts = JSON.stringify(draw_opts);
{% endif %}
var smarts_matches = {};

// Load RDKit
window
.initRDKitModule()
.then(function(RDKit) {
    console.log('RDKit version: ', RDKit.version());
    window.RDKit = RDKit;
    window.RDKitModule = RDKit;

    // Searchbar
    {% include 'js/search.js' %}

    {% if onthefly %}
    {% include 'js/draw_mol.js' %}
    {% endif %}

    // Trigger update to activate tooltips, draw images, setup callbacks...
    listObj.update();
    
    // Set iframe height to fit content.
    fitIframe(window.frameElement);
});