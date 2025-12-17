export async function waitForElement(
    root: HTMLElement,
    querySelector: string
): Promise<void> {
    if (root.querySelector(querySelector)) return
    const observer = new MutationObserver(() => {
        if (root.querySelector(querySelector)) {
            observer.disconnect()
        }
    })
    observer.observe(document.body, {
        childList: true,
        subtree: true,
    })
}

export function debounce<F extends (...args: Parameters<F>) => void>(
    func: F,
    waitFor: number
): (...args: Parameters<F>) => void {
    let timeout: number
    return (...args: Parameters<F>): void => {
        clearTimeout(timeout)
        timeout = setTimeout(() => func(...args), waitFor)
    }
}
