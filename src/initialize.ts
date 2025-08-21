import initRDKitModule from "@rdkit/rdkit"
import type { RDKitModule } from "@rdkit/rdkit"
import $ from "jquery"
import { type AnyModel } from "@anywidget/types"
import { type Placement } from "@floating-ui/dom"
import { type CSSOptions, type WidgetModel } from "./widget"
import { type MolGrid } from "./molgrid"
import { type SmartsMatches } from "./rdkit/smarts"
import { type Callback } from "./interactions/callback"
import { initCellClick } from "./interactions/click"
import { initCheckbox } from "./interactions/select"
import { initKeyboard } from "./interactions/keyboard"
import { initSearch } from "./interactions/search"
import { initToolTip } from "./interactions/tooltips"
import {
    selectAll,
    selectMatching,
    unselectAll,
    invertSelection,
} from "./interactions/select"
import { clipboardCopy, saveCSV, saveSmiles } from "./export"
import { addResizeHandler } from "./interactions/resize"

export let RDKit: RDKitModule | null = null
if (window) {
    // @ts-expect-error
    if (window.RDKitModule) {
        // @ts-expect-error
        RDKit = window.RDKitModule
    } else {
        console.log("Initializing RDKit")
        // @ts-expect-error
        window.RDKitModule = RDKit = await initRDKitModule({
            locateFile: () =>
                "https://unpkg.com/@rdkit/rdkit@2025.3.4-1.0.0/dist/RDKit_minimal.wasm",
        })
        console.log("RDKit version:", RDKit?.version())
    }
}

export function initInteractions(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    supportSelection: boolean,
    smartsMatches: SmartsMatches,
    smilesCol: string,
    searchCols: string[],
    callback: Callback,
    tooltip: boolean,
    tooltipPlacement: Placement | null
) {
    let identifier = model.get("identifier")
    initCellClick(model, supportSelection, callback)
    if (tooltip) {
        initToolTip(identifier, {tooltipPlacement: tooltipPlacement})
    }
    initKeyboard(identifier)
    if (supportSelection) {
        initCheckbox(model, molgrid, smilesCol)
    }
    initSearch(molgrid, smilesCol, searchCols, smartsMatches)

    // Hide pagination if there is only one page.
    // @ts-expect-error
    if (molgrid.listObj.matchingItems.length <= molgrid.listObj.page) {
        $(`#${identifier} .m2g-pagination`).hide()
    } else {
        $(`#${identifier} .m2g-pagination`).show()
    }

    // Add a bunch of phantom cells.
    // These are used as filler to make sure that
    // no grid cells need to be resized when there's
    // not enough results to fill the row.
    $(`#${identifier} .m2g-list`).append(
        '<div class="m2g-cell m2g-phantom"></div>'.repeat(11)
    )

    // Listen to action dropdown.
    $(`#${identifier} .m2g-actions select`).on("change", function (e: JQuery.ChangeEvent) {
        var val = e.target.value
        switch (val) {
            case "select-all":
                selectAll(model, molgrid, smilesCol)
                break
            case "select-matching":
                selectMatching(model, molgrid, smilesCol)
                break
            case "unselect-all":
                unselectAll(model, molgrid)
                break
            case "invert":
                invertSelection(model, molgrid, smilesCol)
                break
            case "copy":
                clipboardCopy(molgrid)
                break
            case "save-smiles":
                saveSmiles(molgrid)
                break
            case "save-csv":
                saveCSV(molgrid)
                break
        }
        $(e.target).val("") // Reset dropdown
    })
}

export function initStyling(el: HTMLElement, css: CSSOptions) {
    el.classList.add("mols2grid-anywidget")
    // update CSS styles
    el.style.setProperty("--m2g-font-family", css.fontFamily)
    el.style.setProperty("--m2g-font-size", css.fontsize)
    el.style.setProperty("--m2g-gap", css.gap === 0 ? "0" : "-" + css.gap)
    el.style.setProperty(
        "--m2g-cell-gap",
        css.gap === 0 ? "-1px -1px 0 0" : `${css.gap}px`
    )
    el.style.setProperty("--m2g-textalign", css.textalign)
    if (css.custom) {
        let style = document.createElement("style")
        style.textContent = css.custom
        el.appendChild(style)
    }
    addResizeHandler(el, css)
}

export function initHeader(el: HTMLElement, header: string) {
    el.innerHTML = header + el.innerHTML
}