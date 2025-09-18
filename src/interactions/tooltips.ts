import {
    computePosition,
    offset,
    autoPlacement,
    shift,
    arrow,
    type Placement,
    type ComputePositionConfig,
} from "@floating-ui/dom"
import { copyOnClick } from "./click"
import { $ } from "../query"


export interface TooltipOptions {
    tooltipPlacement: Placement | null
}

// Show tooltip when hovering the info icon.
export function initToolTip(identifier: string, options: TooltipOptions) {
    $(`#${identifier} .m2g-info`)
        .off("mouseenter")
        .on("mouseenter", ev => {
            const $t = $(<HTMLElement>ev.target).closest(".m2g-cell")
            if ($t.getAttr("m2g-tooltip-active")[0] === "1") {
                // already an existing tooltip
                return
            }
            const referenceEl = $t.elements[0]
            const contentEl = $t.find(".m2g-tooltip").elements[0]
            const tooltip = new Tooltip(
                identifier, <HTMLElement>ev.target, referenceEl, contentEl, options
            )
            $t.setAttr({"m2g-tooltip-active": "1"})
            tooltip.start()
            tooltip.show()
        })
}

const placementMap = new Map([
    ["top", ["bottom", -45]],
    ["right", ["left", 45]],
    ["bottom", ["top", 135]],
    ["left", ["right", -135]],
])

class Tooltip {
    private identifier: string
    private triggerEl: HTMLElement
    private referenceEl: HTMLElement
    private floatingEl: HTMLElement
    private arrowEl: HTMLElement
    private options: TooltipOptions
    private listeners: [string, () => void][]

    public constructor(
        identifier: string,
        triggerEl: HTMLElement,
        referenceEl: HTMLElement,
        contentEl: HTMLElement,
        options: TooltipOptions
    ) {
        this.identifier = identifier
        this.triggerEl = triggerEl
        this.referenceEl = referenceEl
        this.options = options || {}
        this.listeners = []
        this.floatingEl = document.createElement("div")
        this.floatingEl.classList.add("m2g-popover")
        this.floatingEl.innerHTML = <string>contentEl.getAttribute("data-content")
        $(this.referenceEl).closest(`#${identifier}`).append(this.floatingEl)
        this.arrowEl = $(this.floatingEl).find(".arrow").elements[0]
    }

    public show() {
        this.floatingEl.style.display = "block"
        this.update()
    }

    public hide(force: boolean = false) {
        const $t = $(this.referenceEl)
        if ($t.hasClass("m2g-keep-tooltip") && !force) {
            return
        }
        this.floatingEl.style.display = ""
        $t.removeClass("m2g-keep-tooltip")
        this.listeners.forEach(eventListener => {
            this.triggerEl.removeEventListener(...eventListener)
        })
        $(this.floatingEl).find(".copy-me").off("click")
        this.floatingEl.remove()
        $t.setAttr({"m2g-tooltip-active": "0"})
    }

    public click() {
        const $t = $(this.referenceEl)
        $t.toggleClass("m2g-keep-tooltip")
        if ($t.hasClass("m2g-keep-tooltip")) {
            this.show()
        } else {
            this.hide()
        }
    }

    public update() {
        const arrowLen = this.arrowEl.offsetWidth
        const floatingOffset = Math.sqrt(2 * arrowLen ** 2) / 2
        const options: ComputePositionConfig = {
            middleware: [
                offset(floatingOffset),
                // autoPlacement may be inserted here
                shift({ padding: 5 }),
                arrow({ element: this.arrowEl }),
            ],
        }
        if (this.options.tooltipPlacement) {
            options.placement = this.options.tooltipPlacement
        } else {
            // @ts-expect-error
            options.middleware.splice(1, 0, autoPlacement())
        }
        computePosition(this.referenceEl, this.floatingEl, options).then(
            ({ x, y, placement, middlewareData }) => {
                Object.assign(this.floatingEl.style, {
                    left: `${x}px`,
                    top: `${y}px`,
                })
                if (middlewareData.arrow) {
                    const { x: arrowX, y: arrowY } = middlewareData.arrow
                    const [staticSide, rotate] = <[string, number]>(
                        placementMap.get(placement.split("-")[0])
                    )
                    Object.assign(this.arrowEl.style, {
                        left: arrowX != null ? `${arrowX}px` : "",
                        top: arrowY != null ? `${arrowY}px` : "",
                        right: "",
                        bottom: "",
                        [staticSide]: `${-arrowLen / 2}px`,
                        transform: `rotate(${rotate}deg)`,
                    })
                }
            }
        )
    }

    public start() {
        // copy event
        $(this.floatingEl)
            .find(".copy-me")
            .on("click", ev => {
                copyOnClick(<HTMLElement>ev.target)
            })
        // add main events
        let events: [string, () => void][] = [
            ["mouseenter", () => this.show()],
            ["mouseleave", () => this.hide()],
            ["click", () => this.click()],
        ]
        events.forEach(eventListener => {
            this.triggerEl.addEventListener(...eventListener)
            this.listeners.push(eventListener)
        })
        $(`#${this.identifier} .m2g-functions`).on("click", _ => {
            this.hide(true)
        })
    }
}
