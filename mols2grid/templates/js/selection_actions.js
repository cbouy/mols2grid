// check all
$('#btn-chkbox-all').click(function (e) {
    var _id = [];
    var _smiles = [];
    listObj.items.forEach(function (item) {
        if (item.elm) {
            item.elm.getElementsByTagName("input")[0].checked = true;
        } else {
            item.show()
            item.elm.getElementsByTagName("input")[0].checked = true;
            item.hide()
        }
        _id.push(item.values()["mols2grid-id"]);
        _smiles.push(item.values()["data-{{ smiles_col }}"]);
    });
    add_selection({{ grid_id | tojson }}, _id, _smiles);
});
// check matching
$('#btn-chkbox-match').click(function (e) {
    var _id = [];
    var _smiles = [];
    listObj.matchingItems.forEach(function (item) {
        if (item.elm) {
            item.elm.getElementsByTagName("input")[0].checked = true;
        } else {
            item.show()
            item.elm.getElementsByTagName("input")[0].checked = true;
            item.hide()
        }
        _id.push(item.values()["mols2grid-id"]);
        _smiles.push(item.values()["data-{{ smiles_col }}"]);
    });
    add_selection({{ grid_id | tojson }}, _id, _smiles);
});
// uncheck all
$('#btn-chkbox-none').click(function (e) {
    var _id = [];
    listObj.items.forEach(function (item) {
        if (item.elm) {
            item.elm.getElementsByTagName("input")[0].checked = false;
        } else {
            item.show()
            item.elm.getElementsByTagName("input")[0].checked = false;
            item.hide()
        }
        _id.push(item.values()["mols2grid-id"]);
    });
    del_selection({{ grid_id | tojson }}, _id);
});
// invert
$('#btn-chkbox-invert').click(function (e) {
    var _id_add = [];
    var _id_del = [];
    var _smiles = [];
    listObj.items.forEach(function (item) {
        if (item.elm) {
            var chkbox = item.elm.getElementsByTagName("input")[0]
            chkbox.checked = !chkbox.checked;
        } else {
            item.show()
            var chkbox = item.elm.getElementsByTagName("input")[0]
            chkbox.checked = !chkbox.checked;
            item.hide()
        }
        if (chkbox.checked) {
            _id_add.push(item.values()["mols2grid-id"]);
            _smiles.push(item.values()["data-{{ smiles_col }}"]);
        } else {
            _id_del.push(item.values()["mols2grid-id"]);
        }
    });
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