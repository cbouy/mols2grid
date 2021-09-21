// add selection buttons
$("#sort-dropdown").before(`{% include 'html/checkbox_dropdown.html' %}`);
// check all
$('#btn-chkbox-all').click(function (e) {
    var current_page = parseInt($("li.page-item.active > a").text());
    var n_items = {{ n_items_per_page }};
    var last_page = parseInt($("li.page-item > a").last().text());
    for (let i = 0; i < last_page; i++) {
        listObj.show(i * n_items + 1, n_items);
        $("input:checkbox").each(function() {
            this.checked = true;
            $(this).trigger("change");
        });
    }
    listObj.show((current_page - 1) * n_items + 1, n_items);
});
// invert
$('#btn-chkbox-invert').click(function (e) {
    var current_page = parseInt($("li.page-item.active > a").text());
    var n_items = {{ n_items_per_page }};
    var last_page = parseInt($("li.page-item > a").last().text());
    for (let i = 0; i < last_page; i++) {
        listObj.show(i * n_items + 1, n_items);
        $("input:checkbox").each(function() {
        this.checked = !this.checked;
        $(this).trigger("change");
        });
    }
    listObj.show((current_page - 1) * n_items + 1, n_items);
});
// copy to clipboard
$("#btn-chkbox-copy").click(function(e) {
    navigator.clipboard.writeText(SELECTION.to_dict());
});
// export smiles
$("#btn-chkbox-dlsmi").click(function(e) {
    SELECTION.download_smi("selection.smi");
});
// update selection on checkbox click
listObj.on("updated", function (list) {
    $("input:checkbox").change(function() {
        update_current_grid({{ grid_id | tojson }});
        var _id = parseInt($(this).closest(".cell").attr("data-mols2grid-id"));
        var _smiles = $($(this).siblings(".data-{{ smiles_col }}")[0]).text();
        if (this.checked) {
            set_selection(_id, _smiles);
        } else {
            del_selection(_id);
        }
    }); 
});