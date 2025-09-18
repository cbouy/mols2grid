import type { AnyModel } from "@anywidget/types"
import {$} from "../query"
import type { WidgetModel } from "../widget"
import { RDKit } from "../initialize"

export interface Callback {
    callbackFn: string
    callbackType: string
}

interface Data {
    "mols2grid-id": Number
    "img": string
}

export function onCallbackButtonClick(
    target: HTMLElement,
    model: AnyModel<WidgetModel>,
    callback: Callback
) {
    let data = <Data>{}
    data["mols2grid-id"] = parseInt(
        <string>$(target).closest(".m2g-cell").getAttr("data-mols2grid-id")[0]
    )
    data["img"] = $(target).parent.siblings(".data-img").elements[0].innerHTML
    $(target)
        .parent
        .siblings(".data")
        .filter(e => !e.classList.contains(".data-img"))
        .each(el => {
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
