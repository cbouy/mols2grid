import type { MolGrid } from "../molgrid"
import { $ } from "../query"

export interface SortOptions {
    field: string
    order: string
    columns: string[]
}

export function checkboxSort(a: any, b: any, _?: any): number | undefined {
    if (a.elm !== undefined) {
        var checkedA = a.elm.querySelector("input[type=checkbox]").checked
        if (b.elm !== undefined) {
            var checkedB = b.elm.querySelector("input[type=checkbox]").checked
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
    } else if (b.elm !== undefined) {
        return 1
    } else {
        return 0
    }
}

export function addSortingHandler(molgrid: MolGrid, sortOptions: SortOptions) {
    const identifier = molgrid.listObj.listContainer.id
    const sortSelect = <HTMLSelectElement>document.querySelector(`#${identifier} .m2g-sort select`)
    $(sortSelect).on("change", _ => {
        molgrid.sort(sortSelect, true)
    })
    $(`#${identifier} .m2g-order`).on("click", ev => {
        let $t = $(<HTMLElement>ev.target).closest(".m2g-sort")
        $t.removeClass(`m2g-arrow-${sortOptions.order}`)
        sortOptions.order = sortOptions.order === "desc" ? "asc" : "desc"
        $t.addClass(`m2g-arrow-${sortOptions.order}`)
        molgrid.sort(sortSelect, false)
    })
}
