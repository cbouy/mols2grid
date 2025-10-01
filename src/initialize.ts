import initRDKitModule from "@rdkit/rdkit"
import type { RDKitModule } from "@rdkit/rdkit"
import { $ } from "./query"
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
import { initSorting, type SortOptions } from "./interactions/sort"
import { addResizeHandler } from "./interactions/resize"
import { initSelectActions } from "./interactions/select"
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

export function initOnce(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    smartsMatches: SmartsMatches,
    smilesCol: string,
    searchCols: string[],
    sortOptions: SortOptions,
) { 
    initSearch(molgrid, smilesCol, searchCols, smartsMatches)
    initSorting(molgrid, sortOptions)
    initSelectActions(model, molgrid, smilesCol)
}

export function initOnUpdate(
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    supportSelection: boolean,
    smilesCol: string,
    callback: Callback,
    tooltip: boolean,
    tooltipPlacement: Placement | null
) {
    const identifier = model.get("identifier")
    initCellClick(model, supportSelection, callback)
    initKeyboard(identifier)
    if (tooltip) {
        initToolTip(identifier, {tooltipPlacement: tooltipPlacement})
    }
    if (supportSelection) {
        initCheckbox(model, molgrid, smilesCol)
    }

    // Add a bunch of phantom cells.
    // These are used as filler to make sure that
    // no grid cells need to be resized when there's
    // not enough results to fill the row.
    $(`#${identifier} .m2g-list`).append(
        '<div class="m2g-cell m2g-phantom"></div>'.repeat(11)
    )

    // Hide pagination if there is only one page.
    // @ts-expect-error
    if (molgrid.listObj.matchingItems.length <= molgrid.listObj.page) {
        $(`#${identifier} .m2g-pagination`).hide()
    } else {
        $(`#${identifier} .m2g-pagination`).show()
    }
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
    // https://stackoverflow.com/a/75217048
    el.style.setProperty("--m2g-truncate", css.truncate ? "initial" : " ")
    el.style.setProperty("--m2g-no-truncate", css.truncate ? " " : "initial")
    if (css.custom) {
        const style = document.createElement("style")
        style.textContent = css.custom
        el.appendChild(style)
    }
    addResizeHandler(el, css)
}