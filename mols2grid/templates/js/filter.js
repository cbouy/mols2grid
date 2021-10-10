let name = {{ grid_id | tojson }};
if (typeof mols2grid_lists !== "undefined") {
    var listObj = mols2grid_lists[name];
} else if (typeof window.parent.mols2grid_lists !== "undefined") {
    var listObj = window.parent.mols2grid_lists[name];
}
if (typeof listObj !== "undefined") {
    var ids = {{ ids }};
    listObj.filter(function (item) {
        return ids.includes(item.values()["mols2grid-id"]);
    });
}
