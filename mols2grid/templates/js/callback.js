listObj.on("updated", function (list) {
    $('#mols2grid .data-img').click(function(e) {
        var data = {}
        data["mols2grid-id"] = parseInt($(this).closest(".cell")
                                               .attr("data-mols2grid-id"));
        data["img"] = this.innerHTML;
        $(this).siblings(".data").each(function() {
            let name = this.className.split(" ")
                                     .filter(cls => cls.startsWith("data-"))[0]
                                     .substring(5);
            data[name] = this.innerHTML;
        });
        {% if callback_type == "python" %}
        // trigger custom python callback
        let model = window.parent["_MOLS2GRID_" + {{ grid_id | tojson }}];
        if (model) {
            model.set("callback_kwargs", JSON.stringify(data));
            model.save_changes();
        } else {
            // no kernel detected for callback
        }
        {% else %}
        // call custom js callback
        {{ callback }}
        {% endif %}
    });
});