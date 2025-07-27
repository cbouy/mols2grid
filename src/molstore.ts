export class MolStore extends Map {
    zipSet(identifiers: number[], smiles_array: string[]) {
        for (let i = 0; i < identifiers.length; i++) {
            this.set(identifiers[i], smiles_array[i])
        }
    }
    zipDelete(identifiers: number[]) {
        for (let i = 0; i < identifiers.length; i++) {
            this.delete(identifiers[i])
        }
    }
    toDict() {
        var content = "{"
        for (let [key, value] of this) {
            content += key + ":" + JSON.stringify(value) + ","
        }
        content = content.length > 1 ? content.slice(0, -1) : content
        content += "}"
        return content
    }
    downloadSmi(fileName: string, allItems?: Array<any>) {
        var content = ""

        if (allItems) {
            // Gather all smiles
            for (var item of allItems) {
                var smiles = item.values()["data-SMILES"]
                var id = item.values()["mols2grid-id"]
                content += smiles + " " + id + "\n"
            }
        } else {
            // Gather selected smiles
            for (let [key, value] of this) {
                content += value + " " + key + "\n"
            }
        }

        var a = document.createElement("a")
        var file = new Blob([content], { type: "text/plain" })
        a.href = URL.createObjectURL(file)
        a.download = fileName
        a.click()
        a.remove()
    }
}
