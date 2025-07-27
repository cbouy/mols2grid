import type { AnyModel } from "@anywidget/types"
import $ from "jquery"
import type { WidgetModel } from "../widget"
import { RDKit } from "../initialize"

// Store an element's text content in the clipboard.
export function copyOnClick(target: HTMLElement) {
    let text = $(target).text()
    navigator.clipboard.writeText(text)

    // Blink the cell to indicate that the text was copied.
    $(target).addClass("m2g-copy-blink")
    setTimeout(() => {
        $(target).removeClass("m2g-copy-blink")
    }, 450)
}

// Cell click handler.
export function initCellClick(
    model: AnyModel<WidgetModel>,
    supportSelection: boolean,
    callbackType: string,
    callbackFn: string
) {
    $("#mols2grid .m2g-cell")
        .off("click")
        .on("click", function (e: JQuery.ClickEvent) {
            if ($(e.target).hasClass("m2g-info") || $(e.target).is(":checkbox")) {
                // Info button / Checkbox --> do nothing.
            } else if ($(e.target).is("div") && $(e.target).hasClass("data")) {
                // Data string --> copy text.
                copyOnClick(e.target)
            } else if ($(e.target).hasClass("m2g-callback")) {
                // Callback button.
                onCallbackButtonClick(e.target, model, callbackType, callbackFn)
            } else {
                // Outside checkbox --> toggle the checkbox.
                if (supportSelection) {
                    var chkbox: any = $(e.target)
                        .closest(".m2g-cell")
                        .find("input:checkbox")[0]
                    chkbox.checked = !chkbox.checked
                    $(chkbox).trigger("change")
                }
            }
        })
}

// Callback button
function onCallbackButtonClick(
    target: HTMLElement,
    model: AnyModel<WidgetModel>,
    callbackType: string,
    callbackFn: string
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

    if (callbackType == "python") {
        // Trigger custom python callback.
        model.set("callback_kwargs", JSON.stringify(data))
        model.save_changes()
    } else {
        // Call custom js callback.
        const callback = new Function("data", "RDKit", callbackFn)
        callback(data, RDKit)
    }
}
