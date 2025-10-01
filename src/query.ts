type Root = HTMLElement | string | null
type Selector<T extends HTMLElement> = string | T | undefined | (T | null)[]

function findAll(selector: string, root: Root = null): HTMLElement[] {
    let src: Document | HTMLElement
    if (!root) {
        src = document
    }
    else if (typeof root === "string") {
        src = findAll(root)[0]
    } else {
        src = root
    }
    return Array.from(src.querySelectorAll(selector))
}

const eventsCache: Map<[HTMLElement, keyof GlobalEventHandlersEventMap], EventListener> = new Map()

export class Query<T extends HTMLElement = HTMLElement> {
    public selector: Selector<T>
    public elements: T[]

    public constructor(
        selector: Selector<T>, root: Root = null
    ) {
        if (typeof selector === "string") {
            this.elements = <T[]>findAll(selector, root)
        } else if (Array.isArray(selector)) {
            this.elements = <T[]>selector.filter(x => x)
        } else {
            this.elements = [<T>selector]
        }
        if (this.elements.length === 0) {
            console.error("Empty selection for query", selector)
        }
    }
    
    public each(callbackFn: (el: T) => void) {
        this.elements.forEach(callbackFn)
    }
    
    public get parent(): Query<HTMLElement> {
        return new Query(<HTMLElement>this.elements[0].parentElement)
    }
    
    public closest(selector: string): Query<HTMLElement> {
        let v = this.elements[0]
        return new Query(<HTMLElement>v.closest(selector))
    }
    
    public filter(filterFn: (e: T) => boolean) {
        return new Query(this.elements.filter(filterFn))
    }
    
    public index(ix: number): Query<T> {
        return new Query(this.elements[ix])
    }
    
    
    public get children(): Query<HTMLElement> {
        let children: HTMLElement[] = []
        this.each(el => children.push(...<HTMLElement[]>Array.from(el.children)))
        return new Query(children)
    }

    public find(selector: string | null = null): Query<HTMLElement> {
        if (!selector) {
            return this.children
        }
        let matches: HTMLElement[] = []
        this.each(el => matches.push(...<HTMLElement[]>Array.from(el.querySelectorAll(selector))))
        return new Query(matches)
    }

    public siblings(selector: string | null = null): Query<HTMLElement> {
        return this.parent.filter(child => this.elements.some(el => child !== el)).find(selector)
    }

    public append(html: string | HTMLElement) {
        if (typeof html === "string") {
            const template = document.createElement('template')
            template.innerHTML = html
            template.content.childNodes.forEach(
                node => this.elements[0].append(node)
            )
        } else {
            this.elements[0].append(html)
        }
    }

    public get text(): string {
        return <string>this.elements[0].textContent
    }

    public set text(text: string) {
        this.each(el => el.textContent = text)
    }

    public addClass(name: string) {
        this.each(el => el.classList.add(name))
    }

    public removeClass(name: string) {
        this.each(el => el.classList.remove(name))
    }

    public toggleClass(name: string) {
        this.each(el => el.classList.toggle(name))
    }

    public hasClass(name: string): boolean {
        return this.elements.some(el => el.classList.contains(name))
    }

    public getAttr(name: string): (string | null)[] {
        return this.elements.map(el => el.getAttribute(name))
    }

    public setAttr(items: object): Query<T> {
        this.each(el => {
            for (const [key, value] of Object.entries(items)) {
                el.setAttribute(key, value)
            }
        })
        return this
    }
    
    public delAttr(name: string): Query<T> {
        this.each(el => el.removeAttribute(name))
        return this
    }

    public hide(): Query<T> {
        this.each(el => {
            el.style.display = "none"
        })
        return this
    }
    
    public show(): Query<T> {
        this.each(el => {
            el.style.display = ""
        })
        return this
    }

    public off(event: keyof GlobalEventHandlersEventMap): Query<T> {
        navigator.locks.request("eventsCache", (_) => {
            this.each(el => {
                let callbackFn = eventsCache.get([el, event])
                if (callbackFn) {
                    el.removeEventListener(event, callbackFn)
                    eventsCache.delete([el, event])
                }
            })
        })
        return this
    }

    public on<K extends keyof GlobalEventHandlersEventMap, V = GlobalEventHandlersEventMap[K]>(event: K, callbackFn: (ev: V) => void): Query<T> {
        navigator.locks.request("eventsCache", (_) => {
            this.each(el => {
                el.addEventListener(event, <EventListener>callbackFn)
                eventsCache.set([el, event], <EventListener>callbackFn)
            })
        })
        return this
    }

    public change(): Query<T> {
        let ev = new Event("change", {bubbles: true})
        this.each(el => el.dispatchEvent(ev))
        return this
    }

    public previous(): Query<T> {
        return new Query(<T>this.elements[0].previousElementSibling)
    }

    public next(): Query<T> {
        return new Query(<T>this.elements[0].nextElementSibling)
    }

    public css<K extends keyof CSSStyleDeclaration>(name: K): CSSStyleDeclaration[K] {
        return getComputedStyle(this.elements[0])[name]
    }
}

export function $<T extends HTMLElement = HTMLElement>(selector: Selector<T>): Query<T> {
    return new Query(selector)
} 