var sortField = '{{ sort_by }}'
var sortOrder = 'asc'

// Sort dropdown
$('#mols2grid .m2g-sort select').change(sort)

// Sort order
$('#mols2grid .m2g-order').click(flipSort)

function sort(e) {
    if (e) {
        sortField = e.target.value
        var selectedOption = e.target.options[e.target.selectedIndex]
        var sortFieldDisplay = selectedOption.text
    }

    // Sort
    if (sortField == 'checkbox') {
        listObj.sort('mols2grid-id', {order: sortOrder, sortFunction: checkboxSort})
    } else {
        listObj.sort(sortField, {order: sortOrder, sortFunction: mols2gridSortFunction})
    }

    // Update UI.
    $(this).parent().find('.m2g-display').text(sortFieldDisplay)
}

// prettier-ignore
function flipSort() {
    $(this).parent().removeClass('m2d-arrow-' + sortOrder)
    sortOrder = sortOrder === 'desc' ? 'asc' : 'desc'
    $(this).parent().addClass('m2d-arrow-' + sortOrder)
    sort()
}

function mols2gridSortFunction(itemA, itemB, options) {
    var x = itemA.values()[options.valueName]
    var y = itemB.values()[options.valueName]
    if (typeof x === 'number') {
        if (isFinite(x - y)) {
            return x - y
        } else {
            return isFinite(x) ? -1 : 1
        }
    } else {
        x = x ? x.toLowerCase() : x
        y = y ? y.toLowerCase() : y
        return x < y ? -1 : x > y ? 1 : 0
    }
}
function checkboxSort(itemA, itemB, options) {
    if (itemA.elm !== undefined) {
        var checkedA = itemA.elm.querySelector('input[type=checkbox]').checked
        if (itemB.elm !== undefined) {
            var checkedB = itemB.elm.querySelector('input[type=checkbox]').checked
            if (checkedA && !checkedB) {
                return -1
            } else if (!checkedA && checkedB) {
                return 1
            } else {
                return 0
            }
        } else {
            return -1
        }
    } else if (itemB.elm !== undefined) {
        return 1
    } else {
        return 0
    }
}
