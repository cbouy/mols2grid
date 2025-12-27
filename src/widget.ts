import type { AnyModel, RenderProps } from "@anywidget/types"
import { type Placement } from "@floating-ui/dom"
import type List from "list.js-fixed"
import "./widget.css"
import { type ListConfig, MolGrid } from "./molgrid"
import { type SortOptions } from "./interactions/sort"
import { initOnce, initOnUpdate } from "./initialize"
import { type SmartsMatches, type SmartsOptions } from "./rdkit/smarts"
import {
    type MolOptions,
    type DrawOptions,
    initMolDrawing,
    getEmptySvg,
} from "./rdkit/draw"
import { setupHTML } from "./html"
import { type Callback } from "./interactions/callback"
import { $ } from "./query"
import { waitForElement } from "./utils"

export interface WidgetModel {
    options: string
    selection: string
    callback_kwargs: string
    filter_mask: boolean[]
    identifier: string
    name: string
}

export interface CSSOptions {
    fontFamily: string
    fontsize: string
    border: string
    cellWidth: number
    pad: number
    textalign: string
    backgroundColor: string
    hoverColor: string
    gap: number
    truncate: boolean
    custom: string
}

export interface GridConfig {
    listConfig: ListConfig
    smilesCol: string
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
    debug: boolean
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
        debug,
    } = params
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
    createGrid(
        container,
        model,
        supportSelection,
        sortOptions,
        molOptions,
        gridConfig,
        debug
    )
}

function createGrid(
    el: HTMLElement,
    model: AnyModel<WidgetModel>,
    supportSelection: boolean,
    sortOptions: SortOptions,
    molOptions: MolOptions,
    gridConfig: GridConfig,
    debug: boolean
): MolGrid {
    const name = model.get("name")
    const smartsMatches: SmartsMatches = new Map()
    const molgrid = new MolGrid(
        el,
        sortOptions,
        smartsMatches,
        gridConfig.smartsOptions,
        gridConfig.listConfig,
        name,
        debug
    )

    // Restore checkbox state
    const selection: object = JSON.parse(model.get("selection"))
    const cachedSelection: [number[], string[]] = [[], []]
    if (Object.keys(selection).length) {
        Object.entries(selection).forEach(x => {
            cachedSelection[0].push(Number(x[0]))
            cachedSelection[1].push(x[1])
        })
    }
    if (cachedSelection.length) {
        molgrid.store.zipSet(...cachedSelection)
        molgrid.listObj.on("updated", (_: List) => {
            $<HTMLInputElement>('.m2g-cell input[checked="false"]', el).each(
                el => (el.checked = false)
            )
        })
    }

    // Add style for whole cell
    if (gridConfig.wholeCellStyle) {
        molgrid.listObj.on("updated", (_: List) => {
            $("div.m2g-cell", el).each(el => {
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

    waitForElement(el, ".m2g-cell").then(() => {
        // Initialize constant interactions
        initOnce(
            el,
            model,
            molgrid,
            smartsMatches,
            gridConfig.smilesCol,
            gridConfig.searchCols,
            sortOptions,
            molOptions.preferCoordGen
        )

        // Initialize interactions that depend on the underlying data at every update
        const placeholder = gridConfig.onTheFlyRendering
            ? getEmptySvg(gridConfig.drawOptions)
            : ""

        molgrid.listObj.on("updated", function (_: List) {
            initOnUpdate(
                el,
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
                    el,
                    gridConfig.smilesCol,
                    gridConfig.drawOptions,
                    molOptions,
                    smartsMatches,
                    placeholder
                )
            }
        })
        molgrid.listObj.update()
    })
    return molgrid
}

export default { render }
