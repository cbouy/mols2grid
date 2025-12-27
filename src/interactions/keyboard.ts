import { $ } from "../query"

// Keyboard actions.
export function initKeyboard(el: HTMLElement) {
    // Disable scroll when pressing UP/DOWN arrows
    $(".m2g-cell", el)
        .off("keydown")
        .on("keydown", ev => {
            if (ev.code === "ArrowUp" || ev.code === "ArrowDown") {
                ev.preventDefault()
            }
        })

    $(".m2g-cell", el)
        .off("keyup")
        .on("keyup", ev => {
            let $t = $(<HTMLElement>ev.target).closest(".m2g-cell")
            let chkbox = <HTMLInputElement>$t.find("input[type=checkbox]").elements[0]
            switch (ev.code) {
                case "Enter":
                    // ENTER: toggle
                    chkbox.checked = !chkbox.checked
                    $(chkbox).change()
                    break
                case "ArrowLeft":
                    $t.previous().elements[0].focus()
                    break
                case "ArrowRight":
                    $t.next().elements[0].focus()
                    break
                case "ArrowUp":
                case "ArrowDown":
                    let containerWidth = $t.parent.elements[0].offsetWidth
                    let other = $t.elements[0]
                    let cellWidth = other.offsetWidth + parseInt($t.css("marginLeft")) * 2
                    let columns = Math.round(containerWidth / cellWidth)
                    let $cells = $t.parent.children
                    let index = $cells.elements.indexOf(other)
                    let targetIndex =
                        ev.code === "ArrowUp"
                            ? Math.max(index - columns, 0)
                            : Math.min(index + columns, $cells.elements.length)
                    $cells.elements[targetIndex].focus()
                    break
                case "Escape":
                case "Backspace":
                    // ESC/BACKSPACE: unselect
                    chkbox.checked = false
                    $(chkbox).change()
                    break
            }
        })
}
