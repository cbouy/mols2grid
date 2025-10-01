import type { AnyModel } from "@anywidget/types"
import List from "list.js-fixed"
import type { WidgetModel } from "./widget.ts"
import { MolStore } from "./molstore.ts"
import { checkboxSort, type SortOptions } from "./interactions/sort.ts"
import {
    smartsSearchFactory,
    type SmartsOptions,
    type SmartsMatches,
} from "./rdkit/smarts.ts"
import {$} from "./query.ts"

export interface ListConfig {
    valueNames: string[]
    item: string
    page: number
    data: any[]
}

export class MolGrid {
    public listObj: List
    public store: MolStore
    private sortOptions: SortOptions
    private smartsSearchFunc: (query: string, columns: any[]) => void

    public constructor(
        el: HTMLElement,
        sortOptions: SortOptions,
        smartsMatches: SmartsMatches,
        smartsOptions: SmartsOptions,
        config: ListConfig,
    ) {
        this.listObj = new List(
            el,
            {
                listClass: "m2g-list",
                valueNames: config.valueNames,
                item: config.item,
                page: config.page,
                pagination: {
                    paginationClass: "m2g-pagination",
                    item: '<li class="page-item"><a class="page page-link"></a></li>',
                    innerWindow: 1,
                    outerWindow: 1,
                },
            },
            config.data
        )
        this.store = new MolStore()
        this.sortOptions = sortOptions
        this.smartsSearchFunc = smartsSearchFactory(this, smartsOptions, smartsMatches)
    }

    public sort(el: HTMLSelectElement, updateOptions: boolean) {
        if (updateOptions) {
            this.sortOptions.field = el.value
            var sortFieldDisplay = <string>el.options[el.selectedIndex].textContent
            // Update UI.
            $(el).parent.find(".m2g-display").text = sortFieldDisplay
        }
        // Sort
        if (this.sortOptions.field == "checkbox") {
            this.listObj.sort("mols2grid-id", {
                order: this.sortOptions.order,
                sortFunction: checkboxSort,
            })
        } else {
            this.listObj.sort(this.sortOptions.field, { order: this.sortOptions.order })
        }
    }

    public filter(model: AnyModel<WidgetModel>) {
        let filter_mask = model.get("filter_mask")
        if (filter_mask) {
            this.listObj.filter((item: any) => {
                return filter_mask[item.values()["mols2grid-id"]]
            })
        }
    }

    public addSelection(
        model: AnyModel<WidgetModel>,
        identifiers: number[],
        smiles_array: string[]
    ) {
        this.store.zipSet(identifiers, smiles_array)
        model.set("selection", this.store.toDict())
        model.save_changes()
    }

    public delSelection(model: AnyModel<WidgetModel>, identifiers: number[]) {
        this.store.zipDelete(identifiers)
        model.set("selection", this.store.toDict())
        model.save_changes()
    }

    public textSearch(query: string, searchCols: string[]) {
        this.listObj.search(query, searchCols)
    }

    public smartsSearch(query: string, searchCols: string[]) {
        this.listObj.search(
            query,
            searchCols,
            // @ts-expect-error
            this.smartsSearchFunc
        )
    }
}
