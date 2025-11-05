import type { AnyModel } from "@anywidget/types"
import { $ } from "../query"
import type { WidgetModel } from "../widget"
import { onCallbackButtonClick, type Callback } from "./callback"

// Store an element's text content in the clipboard.
export function copyOnClick(target: HTMLElement) {
    let text = $(target).text
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
    callback: Callback
) {
    const identifier = model.get("identifier")
    $(`#${identifier} .m2g-cell`)
        .off("click")
        .on("click", ev => {
            let $t = $(<HTMLElement>ev.target)
            if (
                $t.hasClass("m2g-info") ||
                (<HTMLElement>ev.target).matches("[type=checkbox]")
            ) {
                // Info button / Checkbox --> do nothing.
            } else if ((<HTMLElement>ev.target).matches("div") && $t.hasClass("data")) {
                // Data string --> copy text.
                copyOnClick(<HTMLElement>ev.target)
            } else if ($t.hasClass("m2g-callback")) {
                // Callback button.
                onCallbackButtonClick(<HTMLElement>ev.target, model, callback)
            } else {
                // Outside checkbox --> toggle the checkbox.
                if (supportSelection) {
                    var chkbox = <HTMLInputElement>(
                        $t.closest(".m2g-cell").find("input[type=checkbox]").elements[0]
                    )
                    chkbox.checked = !chkbox.checked
                    $(chkbox).change()
                }
            }
        })
}
