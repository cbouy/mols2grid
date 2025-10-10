import type { AnyModel, RenderProps } from "@anywidget/types"
import "@rdkit/rdkit"
import { type Placement } from "@floating-ui/dom"
import type List from "list.js-fixed"
import "./widget.css"
import { type ListConfig, MolGrid } from "./molgrid"
import { type SortOptions } from "./interactions/sort"
import { initOnce, initOnUpdate } from "./initialize"
import { RDKit } from "./initialize"
import { type SmartsMatches, type SmartsOptions } from "./rdkit/smarts"
import { type MolOptions, type DrawOptions, initMolDrawing } from "./rdkit/draw"
import { setupHTML } from "./html"
import { type Callback } from "./interactions/callback"
import { $ } from "./query"
import {waitForElement} from "./utils"

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
    truncate: boolean
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
}

export interface WidgetOptions {
    supportSelection: boolean
    sortOptions: SortOptions
    molOptions: MolOptions
    gridConfig: GridConfig
    css: CSSOptions
    customHeader: string
}

function render({ model, el }: RenderProps<WidgetModel>) {
    // Render the widget's view into the el HTMLElement.
    const params: WidgetOptions = JSON.parse(model.get("options"))
    let {
        supportSelection,
        sortOptions,
        molOptions,
        gridConfig,
        css,
        customHeader,
    } = params
    RDKit?.prefer_coordgen(molOptions.preferCoordGen)
    const identifier = model.get("identifier")
    el.id = `widget-${identifier}`
    const container = setupHTML(
        el,
        identifier,
        sortOptions.field,
        sortOptions.columns,
        supportSelection,
        css,
        customHeader
    )
    const molgrid = createGrid(
        container,
        model,
        supportSelection,
        sortOptions,
        molOptions,
        gridConfig
    )
    waitForElement(`#${identifier} .m2g-list`).then(molgrid.listObj.update)
}

export function createGrid(
    el: HTMLElement,
    model: AnyModel<WidgetModel>,
    supportSelection: boolean,
    sortOptions: SortOptions,
    molOptions: MolOptions,
    gridConfig: GridConfig
) {
    const identifier = model.get("identifier")
    const gridTarget = <HTMLElement>el.querySelector(`#${identifier}`)
    const smartsMatches: SmartsMatches = new Map()
    const molgrid = new MolGrid(
        gridTarget,
        sortOptions,
        smartsMatches,
        gridConfig.smartsOptions,
        gridConfig.listConfig
    )

    // Restore checkbox state
    if (gridConfig.cachedSelection) {
        molgrid.store.zipSet(...gridConfig.cachedSelection)
        molgrid.listObj.on("updated", (_: List) => {
            $<HTMLInputElement>(`#${identifier} .m2g-cell input[checked="false"]`).each(
                el => el.checked = false
            )
        })
    }

    // Add style for whole cell
    if (gridConfig.wholeCellStyle) {
        molgrid.listObj.on("updated", (_: List) => {
            $(`#${identifier} div.m2g-cell`).each(el => {
                let cellstyle = el.getAttribute("data-cellstyle")
                if (cellstyle) {
                    el.setAttribute("style", cellstyle)
                }
                el.removeAttribute("data-cellstyle")
            })
        })
    }

    // Trigger filtering function on model value change
    model.on("change:filter_mask", function () {
        molgrid.filter(model)
    })

    waitForElement(`#${identifier} .m2g-cell`).then(_ => {
        // Initialize constant interactions
        initOnce(model, molgrid, smartsMatches, gridConfig.smilesCol, gridConfig.searchCols, sortOptions)
    
        // Initialize interactions that depend on the underlying data at every update
        molgrid.listObj.on("updated", function (_: List) {
            initOnUpdate(
                model,
                molgrid,
                supportSelection,
                gridConfig.smilesCol,
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
    })


    return molgrid
}

export default { render }
