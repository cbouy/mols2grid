import type { AnyModel } from "@anywidget/types"
import $ from "jquery"
import type { Callback, WidgetModel } from "../widget"
import { onCallbackButtonClick } from "./callback"

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
    callback: Callback,
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
                onCallbackButtonClick(e.target, model, callback)
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
