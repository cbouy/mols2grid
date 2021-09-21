var iframe = mols2grid_iframe["{{ grid_id }}"];
var listObj = iframe.listObj;
var kl = listObj.items.length;
var ids = {{ ids }};
listObj.filter(function (item) {
    return ids.includes(item.values()["mols2grid-id"]);
});
