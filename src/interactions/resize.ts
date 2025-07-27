import type { CSSOptions } from "../widget"

export function addResizeHandler(el: HTMLElement, css: CSSOptions) {
    // resize events
    // need to be updated JS side as CSS doesn't allow using `calc` in container queries
    // @container (min-width: {{ cell_width * 4 - 1 }}px) and (max-width: {{ cell_width * 6 - 1 }}px) {
    //     #mols2grid .m2g-cell {
    //         flex-basis: calc(100% / 4);
    //     }
    // }
    window.addEventListener("resize", () => {
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
    })
}
