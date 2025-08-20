import type { AnyModel } from "@anywidget/types"
import $ from "jquery"
import type { WidgetModel } from "../widget"
import { RDKit } from "../initialize"

export interface Callback {
    callbackFn: string
    callbackType: string
}

export function onCallbackButtonClick(
    target: HTMLElement,
    model: AnyModel<WidgetModel>,
    callback: Callback
) {
    let data = {}
    // @ts-expect-error
    data["mols2grid-id"] = parseInt(
        // @ts-expect-error
        $(target).closest(".m2g-cell").attr("data-mols2grid-id")
    )
    // @ts-expect-error
    data["img"] = $(target).parent().siblings(".data-img").eq(0).get(0).innerHTML
    $(target)
        .parent()
        .siblings(".data")
        .not(".data-img")
        .each(function (_: number, el: HTMLElement) {
            let name = el.className
                .split(" ")
                .filter(cls => cls.startsWith("data-"))[0]
                .substring(5)
            // @ts-expect-error
            data[name] = el.innerHTML
        })

    if (callback.callbackType == "python") {
        // Trigger custom python callback.
        model.set("callback_kwargs", JSON.stringify(data))
        model.save_changes()
    } else {
        // Call custom js callback.
        const callbackFunction = new Function("data", "RDKit", "$", callback.callbackFn)
        callbackFunction(data, RDKit, $)
    }
}
