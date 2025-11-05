import type { CSSOptions } from "../widget"
import { debounce } from "../utils"

export function addResizeHandler(el: HTMLElement, css: CSSOptions) {
    function resize() {
        const width = window.innerWidth
        function fits(colsMin: number, colsMax: number | null = null): boolean {
            return colsMax
                ? css.cellWidth * colsMin - 1 < width &&
                      width < css.cellWidth * colsMax - 1
                : css.cellWidth * colsMin - 1 < width
        }
        var basis: string | number = "auto"
        if (fits(4, 6)) {
            basis = "calc(100% / 4)"
        } else if (fits(6, 12)) {
            basis = "calc(100% / 6)"
        } else if (fits(12)) {
            basis = "calc(100% / 12)"
        }
        el.style.setProperty("--m2g-flex-basis", basis)
    }
    window.addEventListener("resize", debounce(resize, 100))
    // trigger auto sizing
    resize()
}
