var sort_field = "mols2grid-id";
var sort_order = "asc";
function coerce_type(str) {
    if (isNaN(str)) return null;
    var num = parseFloat(str)
    return isNaN(num) ? str : num;
}
function mols2gridSortFunction(itemA, itemB, options) {
    var x = coerce_type(itemA.values()[options.valueName]);
    var y = coerce_type(itemB.values()[options.valueName]);
    if (x > y) { 
        return 1;
    } else if (x < y) {
        return -1;
    } else {
        return 0;
    }
}
function checkboxSort(itemA, itemB, options) {
    if (itemA.elm !== undefined) {
        var checkedA = itemA.elm.firstChild.checked;
        if (itemB.elm !== undefined) {
            var checkedB = itemB.elm.firstChild.checked;
            if (checkedA && !checkedB) {
                return -1;
            } else if (!checkedA && checkedB) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return -1;
        }
    } else if (itemB.elm !== undefined) {
        return 1;
    } else {
        return 0;
    }
}
$('#mols2grid button.sort-btn').click(function(e) {
    var _field = $(this).attr("data-name");
    if (_field == sort_field) {
        $(this).removeClass("arrow-" + sort_order)
        sort_order = (sort_order === "desc") ? "asc" : "desc";
    } else {
        $('#mols2grid button.sort-btn.active').removeClass("active arrow-" + sort_order);
        sort_order = "asc";
        sort_field = _field;
        $(this).addClass("active");
    }
    $(this).addClass("arrow-" + sort_order)
    if (sort_field == "checkbox") {
        listObj.sort("mols2grid-id", {order: sort_order, sortFunction: checkboxSort});
    } else {
        listObj.sort(_field, {order: sort_order, sortFunction: mols2gridSortFunction});
    }
});