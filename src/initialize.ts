import { type AnyModel } from "@anywidget/types"
import { type Placement } from "@floating-ui/dom"
import { $ } from "./query"
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
import { RDKit } from "./rdkit/loader"

export function initOnce(
    el: HTMLElement,
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    smartsMatches: SmartsMatches,
    smilesCol: string,
    searchCols: string[],
    sortOptions: SortOptions,
    preferCoordGen: boolean
) {
    RDKit.prefer_coordgen(preferCoordGen)
    initSearch(el, molgrid, smilesCol, searchCols, smartsMatches)
    initSorting(el, molgrid, sortOptions)
    initSelectActions(el, model, molgrid, smilesCol)
}

export function initOnUpdate(
    el: HTMLElement,
    model: AnyModel<WidgetModel>,
    molgrid: MolGrid,
    supportSelection: boolean,
    smilesCol: string,
    callback: Callback,
    tooltip: boolean,
    tooltipPlacement: Placement | null
) {
    initCellClick(el, model, supportSelection, callback)
    initKeyboard(el)
    if (tooltip) {
        initToolTip(el, { tooltipPlacement: tooltipPlacement })
    }
    if (supportSelection) {
        initCheckbox(el, model, molgrid, smilesCol)
    }

    // Add a bunch of phantom cells.
    // These are used as filler to make sure that
    // no grid cells need to be resized when there's
    // not enough results to fill the row.
    $(".m2g-list", el).append(
        '<div class="m2g-cell m2g-phantom"></div>'.repeat(11)
    )

    // Hide pagination if there is only one page.
    // @ts-expect-error
    if (molgrid.listObj.matchingItems.length <= molgrid.listObj.page) {
        $(".m2g-pagination", el).hide()
    } else {
        $(".m2g-pagination", el).show()
    }
}

export function initStyling(el: HTMLElement, css: CSSOptions) {
    el.classList.add("mols2grid-anywidget")
    // update CSS styles
    el.style.setProperty("--m2g-font-family", css.fontFamily)
    el.style.setProperty("--m2g-font-size", css.fontsize)
    el.style.setProperty("--m2g-border", css.border)
    el.style.setProperty("--m2g-cell-width", `${css.cellWidth}px`)
    el.style.setProperty("--m2g-pad", `${css.pad}px`)
    el.style.setProperty("--m2g-textalign", css.textalign)
    el.style.setProperty("--m2g-background-color", css.backgroundColor)
    el.style.setProperty("--m2g-hover-color", css.hoverColor)
    el.style.setProperty("--m2g-gap", css.gap === 0 ? "0px" : `-${css.gap}px`)
    el.style.setProperty(
        "--m2g-cell-gap",
        css.gap === 0 ? "-1px -1px 0 0" : `${css.gap}px`
    )
    // https://stackoverflow.com/a/75217048
    el.style.setProperty("--m2g-truncate", css.truncate ? "initial" : " ")
    el.style.setProperty("--m2g-no-truncate", css.truncate ? " " : "initial")
    addResizeHandler(el, css)
}
