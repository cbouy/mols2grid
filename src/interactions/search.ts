import $ from "jquery"
import { type MolGrid } from "../molgrid"
import { type SmartsMatches } from "../rdkit/smarts"

export function initSearch(
    molgrid: MolGrid,
    smilesCol: string,
    searchCols: string[],
    smartsMatches: SmartsMatches
) {
    var searchType = "Text"
    let identifier = molgrid.listObj.listContainer.id
    // Switch search type (Text or SMARTS)
    $(`#${identifier} .m2g-search-options .m2g-option`).on(
        "click",
        function (e: JQuery.ClickEvent) {
            let $t = $(e.target).closest(".m2g-option")
            searchType = $t.text()
            $(`#${identifier} .m2g-search-options .m2g-option.sel`).removeClass("sel")
            $t.addClass("sel")
        }
    )

    // Searchbar update event handler
    $(`#${identifier} .m2g-searchbar`).on("keyup", function (e: JQuery.KeyUpEvent) {
        var query = e.target.value
        if (searchType === "Text") {
            smartsMatches.clear()
            molgrid.textSearch(query, searchCols)
        } else {
            molgrid.smartsSearch(query, [`data-${smilesCol}`])
        }
    })
}
