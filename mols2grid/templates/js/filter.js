let name = {{ grid_id | tojson }};
if ((typeof mols2grid_lists !== "undefined") && (name in mols2grid_lists)) {
    var listObj = mols2grid_lists[name];
    var kl = listObj.items.length;
    var ids = {{ ids }};
    listObj.filter(function (item) {
        return ids.includes(item.values()["mols2grid-id"]);
    });
}
