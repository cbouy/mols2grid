export function waitForElement<T = HTMLElement>(querySelector: string, timeout: number | null = null): Promise<T> {
    return new Promise((resolve, reject) => {
        var timer: number | null = null
        var el: T | null
        if (el = <T>document.querySelector(querySelector)) return resolve(el)
        const observer = new MutationObserver(() => {
            if (el = <T>document.querySelector(querySelector)) {
                observer.disconnect()
                if (timer !== null) clearTimeout(timer)
                return resolve(el)
            }
        })
        observer.observe(document.body, {
            childList: true,
            subtree: true,
        })
        if (timeout) {
            timer = setTimeout(() => {
                observer.disconnect()
                reject()
            }, timeout)
        }
    })
}

export function debounce<F extends (...args: Parameters<F>) => void>(
  func: F,
  waitFor: number,
): (...args: Parameters<F>) => void {
  let timeout: number;
  return (...args: Parameters<F>): void => {
    clearTimeout(timeout);
    timeout = setTimeout(() => func(...args), waitFor);
  };
}
