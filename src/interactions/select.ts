import type { MolGrid } from "../molgrid"
import type { AnyModel } from "@anywidget/types"
import type { WidgetModel } from "../widget"

// Check all.
export function selectAll(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    smilesCol: string
) {
    let identifiers: number[] = []
    let smiles_array: string[] = []
    molgrid.listObj.items.forEach((item: any) => {
        if (item.elm) {
            item.elm.getElementsByTagName("input")[0].checked = true
        } else {
            item.show()
            item.elm.getElementsByTagName("input")[0].checked = true
            item.hide()
        }
        identifiers.push(item.values()["mols2grid-id"])
        smiles_array.push(item.values()[`data-${smilesCol}`])
    })
    molgrid.addSelection(model, identifiers, smiles_array)
}

// Check matching.
export function selectMatching(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    smilesCol: string
) {
    let identifiers: number[] = []
    let smiles_array: string[] = []
    molgrid.listObj.matchingItems.forEach(function (item: any) {
        if (item.elm) {
            item.elm.getElementsByTagName("input")[0].checked = true
        } else {
            item.show()
            item.elm.getElementsByTagName("input")[0].checked = true
            item.hide()
        }
        identifiers.push(item.values()["mols2grid-id"])
        smiles_array.push(item.values()[`data-${smilesCol}`])
    })
    molgrid.addSelection(model, identifiers, smiles_array)
}

// Uncheck all.
export function unselectAll(model: AnyModel<WidgetModel>, molgrid: MolGrid) {
    let identifiers: number[] = []
    molgrid.listObj.items.forEach(function (item: any) {
        if (item.elm) {
            item.elm.getElementsByTagName("input")[0].checked = false
        } else {
            item.show()
            item.elm.getElementsByTagName("input")[0].checked = false
            item.hide()
        }
        identifiers.push(item.values()["mols2grid-id"])
    })
    molgrid.delSelection(model, identifiers)
}

// Invert selection.
export function invertSelection(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    smilesCol: string
) {
    let identifiers_add: number[] = []
    let identifiers_del: number[] = []
    let smiles_array: string[] = []
    molgrid.listObj.items.forEach(function (item: any) {
        if (item.elm) {
            var chkbox = item.elm.getElementsByTagName("input")[0]
            chkbox.checked = !chkbox.checked
        } else {
            item.show()
            var chkbox = item.elm.getElementsByTagName("input")[0]
            chkbox.checked = !chkbox.checked
            item.hide()
        }
        if (chkbox.checked) {
            identifiers_add.push(item.values()["mols2grid-id"])
            smiles_array.push(item.values()[`data-${smilesCol}`])
        } else {
            identifiers_del.push(item.values()["mols2grid-id"])
        }
    })
    molgrid.delSelection(model, identifiers_del)
    molgrid.addSelection(model, identifiers_add, smiles_array)
}

// Update selection on checkbox click.
export function initCheckbox(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    smilesCol: string
) {
    $("input:checkbox")
        .off("change")
        .on("change", function (e: JQuery.ChangeEvent) {
            var identifier = parseInt(
                // @ts-expect-error
                $(e.target).closest(".m2g-cell").attr("data-mols2grid-id")
            )
            if (e.target.checked) {
                var smiles = $(
                    $(e.target).closest(".m2g-cell").children(`.data-${smilesCol}`)[0]
                ).text()
                molgrid.addSelection(model, [identifier], [smiles])
            } else {
                molgrid.delSelection(model, [identifier])
            }
        })
}
