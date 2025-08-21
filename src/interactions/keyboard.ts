import $ from "jquery"

// Keyboard actions.
export function initKeyboard(identifier: string) {
    // Disable scroll when pressing UP/DOWN arrows
    $(`#${identifier} .m2g-cell`)
        .off("keydown")
        .on("keydown", function (e) {
            if (e.which == 38 || e.which == 40) {
                e.preventDefault()
            }
        })

    $(`#${identifier} .m2g-cell`)
        .off("keyup")
        .on("keyup", function (e: JQuery.KeyUpEvent) {
            let $t = $(e.target).closest(".m2g-cell")
            let chkbox: any = $t.find("input:checkbox")[0]
            if (e.which == 13) {
                // ENTER: toggle
                chkbox.checked = !chkbox.checked
                $(chkbox).trigger("change")
            } else if (e.which == 27 || e.which == 8) {
                // ESC/BACKSPACE: unselect
                chkbox.checked = false
                $(chkbox).trigger("change")
            } else if (e.which == 37) {
                // LEFT
                $t.prev().trigger("focus")
            } else if (e.which == 39) {
                // RIGHT
                $t.next().trigger("focus")
            } else if (e.which == 38 || e.which == 40) {
                var containerWidth = $t.parent().outerWidth()
                var cellWidth =
                    // @ts-expect-error
                    $t.outerWidth() + parseInt($t.css("marginLeft")) * 2
                // @ts-expect-error
                var columns = Math.round(containerWidth / cellWidth)
                var index = $t.index()
                if (e.which == 38) {
                    // UP
                    var indexAbove = Math.max(index - columns, 0)
                    $t.parent().children().eq(indexAbove).trigger("focus")
                } else if (e.which == 40) {
                    // DOWN
                    var total = $t.parent().children().length
                    var indexBelow = Math.min(index + columns, total)
                    $t.parent().children().eq(indexBelow).trigger("focus")
                }
            }
        })
}
