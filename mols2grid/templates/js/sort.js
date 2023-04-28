var sort_field = '{{ sort_by }}'
var sort_order = 'asc'

// Sort dropdown
$('#dropdown-sort').change(function (e) {
    var _field = e.target.value
    if (_field == sort_field) {
        $(this).removeClass('arrow-' + sort_order)
        sort_order = sort_order === 'desc' ? 'asc' : 'desc'
    } else {
        sort_order = 'asc'
        sort_field = _field
    }
    if (sort_field == 'checkbox') {
        listObj.sort('mols2grid-id', { order: sort_order, sortFunction: checkboxSort })
    } else {
        listObj.sort(_field, { order: sort_order, sortFunction: mols2gridSortFunction })
    }
})

// // Legacy sort – kept here for reference of arrow functionality
// $('#mols2grid button.sort-btn').click(function (e) {
//     var _field = $(this).attr('data-name')
//     if (_field == sort_field) {
//         $(this).removeClass('arrow-' + sort_order)
//         sort_order = sort_order === 'desc' ? 'asc' : 'desc'
//     } else {
//         $('#mols2grid button.sort-btn.active').removeClass('active arrow-' + sort_order)
//         sort_order = 'asc'
//         sort_field = _field
//         $(this).addClass('active')
//     }
//     $(this).addClass('arrow-' + sort_order)
//     if (sort_field == 'checkbox') {
//         listObj.sort('mols2grid-id', { order: sort_order, sortFunction: checkboxSort })
//     } else {
//         listObj.sort(_field, { order: sort_order, sortFunction: mols2gridSortFunction })
//     }
// })

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
        var checkedA = itemA.elm.firstChild.checked
        if (itemB.elm !== undefined) {
            var checkedB = itemB.elm.firstChild.checked
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
