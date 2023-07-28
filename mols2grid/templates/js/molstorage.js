class MolStorage extends Map {
    multi_set(_id, _smiles) {
        for (let i = 0; i < _id.length; i++) {
            this.set(_id[i], _smiles[i])
        }
    }
    multi_del(_id) {
        for (let i = 0; i < _id.length; i++) {
            this.delete(_id[i])
        }
    }
    to_dict() {
        var content = '{'
        for (let [key, value] of this) {
            content += key + ':' + JSON.stringify(value) + ','
        }
        content = content.length > 1 ? content.slice(0, -1) : content
        content += '}'
        return content
    }
    to_keys() {
        var content = []
        for (let [key] of this) {
            content.push(key)
        }
        return content
    }
    download_smi(fileName, allItems) {
        var content = ''

        if (allItems) {
            // Gather all smiles
            for (var item of allItems) {
                var smiles = item.values()['data-SMILES']
                var id = item.values()['mols2grid-id']
                content += smiles + ' ' + id + '\n'
            }
        } else {
            // Gather selected smiles
            for (let [key, value] of this) {
                content += value + ' ' + key + '\n'
            }
        }

        var a = document.createElement('a')
        var file = new Blob([content], { type: 'text/plain' })
        a.href = URL.createObjectURL(file)
        a.download = fileName
        a.click()
        a.remove()
    }
}
