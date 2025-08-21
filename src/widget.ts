import type { AnyModel, RenderProps } from "@anywidget/types"
import "@rdkit/rdkit"
import { type Placement } from "@floating-ui/dom"
import $ from "jquery"
import type List from "list.js-fixed"
import "./widget.css"
import { type ListConfig, MolGrid } from "./molgrid"
import { type SortOptions } from "./interactions/sort"
import { initInteractions } from "./initialize"
import { RDKit } from "./initialize"
import { type SmartsMatches, type SmartsOptions } from "./rdkit/smarts"
import { type MolOptions, type DrawOptions, initMolDrawing } from "./rdkit/draw"
import { makeHTML } from "./html"
import { initStyling, initHeader } from "./initialize"
import { addSortingHandler } from "./interactions/sort"
import { type Callback } from "./interactions/callback"

export interface WidgetModel {
    options: string
    selection: string
    callback_kwargs: string
    filter_mask: boolean[]
    identifier: string
}

export interface CSSOptions {
    fontFamily: string
    fontsize: string
    cellWidth: number
    gap: number
    border: string
    pad: number
    textalign: string
    backgroundColor: string
    hoverColor: string
    custom: string
}

export interface GridConfig {
    listConfig: ListConfig
    smilesCol: string
    cachedSelection: [number[], string[]]
    wholeCellStyle: boolean
    tooltip: boolean
    tooltipPlacement: Placement | null
    callback: Callback
    onTheFlyRendering: boolean
    drawOptions: DrawOptions
    smartsOptions: SmartsOptions
    searchCols: string[]
    customHeader: string
    css: CSSOptions
}

export interface WidgetOptions {
    supportSelection: boolean
    sortOptions: SortOptions
    molOptions: MolOptions
    gridConfig: GridConfig
}

function render({ model, el }: RenderProps<WidgetModel>) {
    // Render the widget's view into the el HTMLElement.
    let params: WidgetOptions = JSON.parse(model.get("options"))
    let { supportSelection, sortOptions, molOptions, gridConfig } = params
    RDKit?.prefer_coordgen(molOptions.preferCoordGen)
    let uuid = crypto.randomUUID()
    el.id = `m2g-widget-${uuid}`
    let identifier = `m2g-${uuid}`
    model.set("identifier", identifier)
    model.save_changes()
    let html = makeHTML(identifier, sortOptions.field, sortOptions.columns, supportSelection)
    el.innerHTML = html
    let molgrid = createGrid(
        el,
        model,
        supportSelection,
        sortOptions,
        molOptions,
        gridConfig
    )
    molgrid.listObj.update()
}

export function createGrid(
    el: HTMLElement,
    model: AnyModel<WidgetModel>,
    supportSelection: boolean,
    sortOptions: SortOptions,
    molOptions: MolOptions,
    gridConfig: GridConfig
) {
    let identifier = model.get("identifier")
    initStyling(el, gridConfig.css)
    if (gridConfig.customHeader) {
        initHeader(el, gridConfig.customHeader)
    }
    let gridTarget = <HTMLElement>el.querySelector(`#${identifier}`)
    let smartsMatches: SmartsMatches = new Map()
    let molgrid = new MolGrid(
        gridTarget,
        sortOptions,
        smartsMatches,
        gridConfig.smartsOptions,
        gridConfig.listConfig
    )

    // Trigger filtering function on model value change
    model.on("change:filter_mask", function () {
        molgrid.filter(model)
    })
    // Restore checkbox state
    if (gridConfig.cachedSelection) {
        molgrid.store.zipSet(...gridConfig.cachedSelection)
        molgrid.listObj.on("updated", (_: List) => {
            $(`#${identifier} .m2g-cell input[checked="false"]`).prop("checked", false)
        })
    }
    // Sorting
    addSortingHandler(molgrid, sortOptions)

    // Add style for whole cell
    if (gridConfig.wholeCellStyle) {
        molgrid.listObj.on("updated", (_: List) => {
            $(`#${identifier} div.m2g-cell`).each(function (_: number, el: HTMLElement) {
                var $t = $(el)
                $t.attr({ style: $t.attr("data-cellstyle") }).removeAttr(
                    "data-cellstyle"
                )
            })
        })
    }

    // (Re)initialize all grid interaction every time the grid changes.
    molgrid.listObj.on("updated", function (_: List) {
        initInteractions(
            model,
            molgrid,
            supportSelection,
            smartsMatches,
            gridConfig.smilesCol,
            gridConfig.searchCols,
            gridConfig.callback,
            gridConfig.tooltip,
            gridConfig.tooltipPlacement
        )
        if (gridConfig.onTheFlyRendering) {
            initMolDrawing(
                identifier,
                gridConfig.smilesCol,
                gridConfig.drawOptions,
                molOptions,
                smartsMatches
            )
        }
    })
    return molgrid
}

export default { render }
