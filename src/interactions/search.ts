import { $ } from "../query"
import { type MolGrid } from "../molgrid"
import { type SmartsMatches } from "../rdkit/smarts"

export function initSearch(
    el: HTMLElement,
    molgrid: MolGrid,
    smilesCol: string,
    searchCols: string[],
    smartsMatches: SmartsMatches
) {
    var searchType = "Text"

    // Switch search type (Text or SMARTS)
    $(".m2g-search-options .m2g-option", el).on("click", ev => {
        console.log(ev.target)
        let $t = $(<HTMLElement>ev.target).closest(".m2g-option")
        searchType = $t.text
        $(".m2g-search-options .m2g-option.sel", el).removeClass("sel")
        $t.addClass("sel")
    })

    // Searchbar update event handler
    $<HTMLInputElement>(".m2g-searchbar", el).on("keyup", ev => {
        let query = (<HTMLInputElement>ev.target).value
        smartsMatches.clear()
        if (searchType === "Text") {
            molgrid.textSearch(query, searchCols)
        } else {
            molgrid.smartsSearch(query, [`data-${smilesCol}`])
        }
    })
}
