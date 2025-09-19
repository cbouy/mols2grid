import { $ } from "../query"
import { type MolGrid } from "../molgrid"
import { type SmartsMatches } from "../rdkit/smarts"
import { debounce } from "../utils"

export function initSearch(
    molgrid: MolGrid,
    smilesCol: string,
    searchCols: string[],
    smartsMatches: SmartsMatches
) {
    var searchType = "Text"
    const identifier = molgrid.listObj.listContainer.id
    // Switch search type (Text or SMARTS)
    $(`#${identifier} .m2g-search-options .m2g-option`)
        .on("click", ev => {
            let $t = $(<HTMLElement>ev.target).closest(".m2g-option")
            searchType = $t.text
            $(`#${identifier} .m2g-search-options .m2g-option.sel`).removeClass("sel")
            $t.addClass("sel")
        }
    )

    // Searchbar update event handler
    $<HTMLInputElement>(`#${identifier} .m2g-searchbar`)
        .on("keyup", debounce(
            (ev: Event) => {
                let query = (<HTMLInputElement>ev.target).value
                if (searchType === "Text") {
                    smartsMatches.clear()
                    molgrid.textSearch(query, searchCols)
                } else {
                    molgrid.smartsSearch(query, [`data-${smilesCol}`])
                }
            },
            300,
        )
    )
}
