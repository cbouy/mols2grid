// check all
$('#btn-chkbox-all').click(function (e) {
    var current_page = parseInt($("li.page-item.active > a").text());
    var n_items = {{ n_items_per_page }};
    var last_page = parseInt($("li.page-item > a").last().text());
    var _id = [];
    var _smiles = [];
    for (let i = 0; i < last_page; i++) {
        listObj.show(i * n_items + 1, n_items);
        $("input:checkbox").each(function() {
            this.checked = true;
            _id.push(parseInt($(this).closest(".cell").attr("data-mols2grid-id")));
            _smiles.push($($(this).siblings(".data-{{ smiles_col }}")[0]).text());
        });
    }
    listObj.show((current_page - 1) * n_items + 1, n_items);
    add_selection({{ grid_id | tojson }}, _id, _smiles);
});
// uncheck all
$('#btn-chkbox-none').click(function (e) {
    var current_page = parseInt($("li.page-item.active > a").text());
    var n_items = {{ n_items_per_page }};
    var last_page = parseInt($("li.page-item > a").last().text());
    var _id = [];
    for (let i = 0; i < last_page; i++) {
        listObj.show(i * n_items + 1, n_items);
        $("input:checkbox:checked").each(function() {
            this.checked = false;
            _id.push(parseInt($(this).closest(".cell").attr("data-mols2grid-id")));
        });
    }
    listObj.show((current_page - 1) * n_items + 1, n_items);
    del_selection({{ grid_id | tojson }}, _id);
});
// invert
$('#btn-chkbox-invert').click(function (e) {
    var current_page = parseInt($("li.page-item.active > a").text());
    var n_items = {{ n_items_per_page }};
    var last_page = parseInt($("li.page-item > a").last().text());
    var _id_add = [];
    var _id_del = [];
    var _smiles = [];
    for (let i = 0; i < last_page; i++) {
        listObj.show(i * n_items + 1, n_items);
        $("input:checkbox").each(function() {
            this.checked = !this.checked;
            var _id = parseInt($(this).closest(".cell").attr("data-mols2grid-id"));
            if (this.checked) {
                _id_add.push(_id);
                _smiles.push($(this).siblings(".data-{{ smiles_col }}").first().text());
            } else {
                _id_del.push(_id);
            }
        });
    }
    listObj.show((current_page - 1) * n_items + 1, n_items);
    del_selection({{ grid_id | tojson }}, _id_del);
    add_selection({{ grid_id | tojson }}, _id_add, _smiles);
});
// copy to clipboard
$("#btn-chkbox-copy").click(function(e) {
    navigator.clipboard.writeText(SELECTION.to_dict());
});
// export smiles
$("#btn-chkbox-dlsmi").click(function(e) {
    SELECTION.download_smi("selection.smi");
});
// export CSV
$("#btn-chkbox-dlcsv").click(function(e) {
    var sep = "\t"
    // same order as subset + tooltip
    var columns = Array.from(listObj.items[0].elm.querySelectorAll("div.data"))
                       .map(elm => elm.classList[1]);
    // remove 'data-'
    var header = columns.map(name => name.slice(5));
    // csv content
    header = ["index"].concat(header).join(sep);
    var content = header + "\n";
    for (let [index, smiles] of SELECTION.entries()) {
        var data = listObj.items[index].values();
        content += index;
        columns.forEach((key) => {
            content += sep + data[key];
        })
        content += "\n";
    }
    var a = document.createElement("a");
    var file = new Blob([content], {type: "text/csv"});
    a.href = URL.createObjectURL(file);
    a.download = "selection.csv";
    a.click();
    a.remove();
});
// update selection on checkbox click
listObj.on("updated", function (list) {
    $("input:checkbox").change(function() {
        var _id = parseInt($(this).closest(".cell").attr("data-mols2grid-id"));
        if (this.checked) {
            var _smiles = $($(this).siblings(".data-{{ smiles_col }}")[0]).text();
            add_selection({{ grid_id | tojson }}, [_id], [_smiles]);
        } else {
            del_selection({{ grid_id | tojson }}, [_id]);
        }
    }); 
});