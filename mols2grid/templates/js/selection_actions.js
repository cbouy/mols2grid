listObj.on("updated", initInteraction);

// (Re)initialiuze all grid interaction every time the grid changes.
function initInteraction(list) {
    initCellClick()
    initToolTip()
    initCheckbox()
}

// Cell click handler.
function initCellClick() {
    $('#mols2grid .m2g-cell').click(function(e) {
        if ($(e.target).hasClass('m2g-info')) {
            var isVisible = $('div.popover[role=tooltip]').length
        } else if ($(e.target).hasClass('data')) {
            // Copy text when clicking a data string.
            var text = $(e.target).text()
            navigator.clipboard.writeText(text)

            // Blink the cell to indicate that the text was copied.
            $(e.target).addClass('m2g-copy-blink')
            setTimeout(function() {
                $(e.target).removeClass('m2g-copy-blink')
            }, 450)
        } else if (!$(e.target).is(':checkbox')) {
            // When clicking anywhere outside the checkbox, toggle the checkbox.
            var chkbox = $(this).find('input:checkbox')[0]
            chkbox.checked = !chkbox.checked
            $(chkbox).trigger('change')
        }
    })
}

// Show tooltip when hovering the info icon.
function initToolTip() {
    $('#mols2grid .m2g-info').mouseenter(function() {
        // Show on enter
        $(this).parent().find('.mols2grid-tooltip[data-toggle="popover"]').popover('show')
    }).mouseleave(function() {
        // Hide on leave, unless sticky.
        if (!$(this).parent().hasClass('m2g-keep-tooltip')) {
            $(this).parent().find('.mols2grid-tooltip[data-toggle="popover"]').popover('hide')
        }
    }).click(function() {
        // Toggle sticky on click.
        $(this).parent().toggleClass('m2g-keep-tooltip')

        // Hide tooltip when sticky was turned off.
        if ($(this).parent().hasClass('m2g-keep-tooltip')) {
            $(this).parent().find('.mols2grid-tooltip[data-toggle="popover"]').popover('show')
        } else if (!$(this).parent().hasClass('m2g-keep-tooltip')) {
            $(this).parent().find('.mols2grid-tooltip[data-toggle="popover"]').popover('hide')
        }
    })
}

// Update selection on checkbox click.
function initCheckbox() {
    $("input:checkbox").change(function() {
        var _id = parseInt($(this).closest(".m2g-cell").attr("data-mols2grid-id"));
        if (this.checked) {
            var _smiles = $($(this).siblings(".data-{{ smiles_col }}")[0]).text();
            add_selection({{ grid_id | tojson }}, [_id], [_smiles]);
        } else {
            del_selection({{ grid_id | tojson }}, [_id]);
        }
    });
}



/**
 * Actions
 */

// Listen to action dropdown.
$('#mols2grid .m2g-actions select').change(function(e) {
    var val = e.target.value
    switch(val) {
        case 'select-all':
            selectAll()
            break
        case 'select-matching':
            selectMatching()
            break
        case 'unselect-all':
            unselectAll()
            break
        case 'invert':
            invertSelection()
            break
        case 'copy':
            copy()
            break
        case 'save-smiles':
            saveSmiles()
            break
        case 'save-csv':
            saveCSV()
            break
    }
    $(this).val('') // Reset dropdown
})

// Check all.
function selectAll(e) {
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
};


// Check matching.
function selectMatching(e) {
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
};

// Uncheck all.
function unselectAll(e) {
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
};

// Invert selection.
function invertSelection(e) {
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
};

// Copy to clipboard.
function copy(e) {
    navigator.clipboard.writeText(SELECTION.to_dict());
};

// Export smiles.
function saveSmiles(e) {
    SELECTION.download_smi("selection.smi");
};

// Export CSV.
function saveCSV(e) {
    var sep = "\t"
    // Same order as subset + tooltip
    var columns = Array.from(listObj.items[0].elm.querySelectorAll("div.data"))
                       .map(elm => elm.classList[1])
                       .filter(name => name !== "data-img");
    // Remove 'data-' and img
    var header = columns.map(name => name.slice(5));
    // CSV content
    header = ["index"].concat(header).join(sep);
    var content = header + "\n";
    listObj.items.forEach(function (item) {
        let data = item.values();
        let index = data["mols2grid-id"];
        if (SELECTION.has(index)) {
            content += index;
            columns.forEach((key) => {
                content += sep + data[key];
            })
            content += "\n";
        }
    });
    var a = document.createElement("a");
    var file = new Blob([content], {type: "text/csv"});
    a.href = URL.createObjectURL(file);
    a.download = "selection.csv";
    a.click();
    a.remove();
};