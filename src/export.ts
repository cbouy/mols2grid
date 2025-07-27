import { type MolGrid } from "./molgrid"

// Copy to clipboard.
export function clipboardCopy(molgrid: MolGrid) {
    let content = _renderCSV(molgrid, "\t")
    navigator.clipboard.writeText(content)
}

// Export smiles.
export function saveSmiles(molgrid: MolGrid) {
    var fileName = "selection.smi"
    if (molgrid.store.size) {
        // Download selected smiles
        molgrid.store.downloadSmi(fileName)
    } else {
        // Download all smiles
        molgrid.store.downloadSmi(fileName, molgrid.listObj.items)
    }
}

// Export CSV.
export function saveCSV(molgrid: MolGrid) {
    let content = _renderCSV(molgrid, ";")
    var a = document.createElement("a")
    var file = new Blob([content], { type: "text/csv" })
    a.href = URL.createObjectURL(file)
    a.download = "selection.csv"
    a.click()
    a.remove()
}

// Render CSV for export of clipboard.
function _renderCSV(molgrid: MolGrid, sep: string) {
    // Same order as subset + tooltip
    // @ts-expect-error
    let columns = Array.from(molgrid.listObj.items[0].elm.querySelectorAll("div.data"))
        .map((elm: any) => elm.classList[1])
        .filter(name => name !== "data-img")
    // Remove 'data-' and img
    let headerColumns = columns.map(name => name.slice(5))
    // CSV content
    let header = ["index"].concat(headerColumns).join(sep)
    var content = header + "\n"
    molgrid.listObj.items.forEach((item: any) => {
        let data = item.values()
        let index = data["mols2grid-id"]
        if (molgrid.store.has(index) || molgrid.store.size === 0) {
            content += index
            columns.forEach(key => {
                content += sep + data[key]
            })
            content += "\n"
        }
    })
    return content
}
