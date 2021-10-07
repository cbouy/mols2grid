listObj.on("updated", function (list) {
    $('#mols2grid .data-img').click(function(e) {
        var data = {}
        data["mols2grid-id"] = parseInt($(this).closest(".cell")
                                               .attr("data-mols2grid-id"));
        $(this).siblings(".data").each(function() {
            let name = this.className.split(" ")
                                     .filter(cls => cls.startsWith("data-"))[0]
                                     .substring(5);
            data[name] = this.innerHTML;
        });
        {% if callback_type == "python" %}
        // call python func
        if (kernel_env === "jupyter") {
            kernel.execute("{{ callback }}("+JSON.stringify(data)+")");
        } else if (kernel_env === "colab") {
            (async function() {
                const result = await kernel.invokeFunction("{{ callback }}",
                                                           [data], {});
            })();
        }
        {% else %}
        // custom js callback
        {{ callback }}
        {% endif %}
    });
});