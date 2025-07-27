// Set iframe height to fit content.
export function fitIframe(iframe: HTMLIFrameElement, iframePadding: number) {
    // Only fit height when no specific height was given.
    if (iframe.getAttribute("height")) return

    // Initial fit + refit whenever the window width changes.
    var prevInnerWidth: number | null = null
    fitToHeight()
    $(window).on("resize", () => {
        if (window.innerWidth != prevInnerWidth) {
            prevInnerWidth = window.innerWidth
            fitToHeight()
        }
    })

    // Fit iframe height to content height.
    function fitToHeight() {
        var height = iframe?.contentDocument?.body?.scrollHeight
        iframe.style.height = (height || 0) + iframePadding + "px"
    }
}

export interface IFrameOptions {
    enabled: boolean
    width: number
    height: number
    padding: number
    allow: string
    sandbox: string
}

export function wrapIframe(
    el: HTMLElement,
    doc: string,
    options: IFrameOptions
): HTMLElement {
    let style = document.createElement("style")
    let css = document.createTextNode(`<style>
    /* Some CSS to integrate with Jupyter more cleanly */
    div.output_subarea {
        /* Undo an unfortunate max-width parameter
        that causes the output area to be too narrow
        on smaller screens. */
        max-width: none;
        /* Align the table with the content */
        padding: 0;
        /* Let it breathe */
        margin-top: 20px;
    }
    </style>`)
    style.appendChild(css)
    el.appendChild(style)
    let iframe = document.createElement("iframe")
    iframe.classList.add("mols2grid-iframe")
    iframe.width = options.width.toString()
    if (options.height) iframe.height = options.height.toString()
    if (options.allow) iframe.allow = options.allow
    if (options.sandbox) iframe.sandbox = options.sandbox
    iframe.srcdoc = doc
    el.appendChild(iframe)
    return iframe
}
