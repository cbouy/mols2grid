class MolStorage extends Map {
    to_dict() {
        var content = "{";
        for (let [key, value] of this) {
            content += key + ":" + JSON.stringify(value) + ",";
        }
        content = content.length > 1 ? content.slice(0, -1) : content;
        content += "}";
        return content
    }
    download_smi(fileName) {
        var content = "SMILES index\n";
        for (let [key, value] of this) {
            content += value + " " + key + "\n";
        }
        var a = document.createElement("a");
        var file = new Blob([content], {type: "text/plain"});
        a.href = URL.createObjectURL(file);
        a.download = fileName;
        a.click();
    }
}