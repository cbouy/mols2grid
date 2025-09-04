import { type CSSOptions } from "./widget"
import { initStyling } from "./initialize"

export function setupHTML(
    container: HTMLElement,
    identifier: string,
    sortBy: string,
    sortCols: string[],
    supportSelection: boolean,
    css: CSSOptions,
    customHeader: string,
): HTMLElement {
    let defaultSort = sortBy.replace("data-", "")
    let sortDisplayValue = defaultSort == "mols2grid-id" ? "Index" : defaultSort
    let chunks: string[] = []
    chunks.push(`<div id="${identifier}">
    <!-- Pagination & search -->
    <div class="m2g-functions">
        <!-- Rows are used to collapse functions into two rows on smaller screens -->
        <div class="m2g-row">
            <!-- Pagination -->
            <ul class="m2g-pagination" class="d-flex"></ul>
            <div class="m2g-gap"></div>

            <!-- Sort dropdown -->
            <div class="m2g-dropdown m2g-sort">
                <select>`)
    let i = 0
    sortCols.forEach((name: string) => {
        let displayName = name.replace("data-", "")
        let selected = defaultSort === displayName ? " selected" : ""
        if (i === 0) {
            chunks.push(`<option value="mols2grid-id"${selected}>Index</option>`)
        } else {
            chunks.push(`<option value="${name}"${selected}>${displayName}</option>`)
        }
        i += 1
    })
    if (supportSelection) {
        chunks.push(`<option value="checkbox">Selected</option>`)
    }
    chunks.push(`                </select>
                <div class="m2g-order"></div>
                <div class="m2g-display">${sortDisplayValue}</div>
            </div>
        </div>
        <div class="m2g-row">
            <!-- Search bar -->
            <div class="m2g-search-wrap">
                <input
                    type="text"
                    class="m2g-searchbar form-control"
                    placeholder="Search"
                    aria-label="Search"
                    aria-describedby="basic-addon1"
                />
                <div class="m2g-search-options">
                    <div class="m2g-option m2g-search-text sel">Text</div>
                    <div class="m2g-option m2g-search-smarts">SMARTS</div>
                </div>
            </div>

            <!-- Action dropdown -->
            <div class="m2g-dropdown m2g-actions">
                <select>
                    <option hidden>-</option>
                    <option value="select-all">Select all</option>
                    <option value="select-matching">Select matching</option>
                    <option value="unselect-all">Unselect all</option>
                    <option value="invert">Invert</option>
                    <option value="copy">Copy to clipboard</option>
                    <option value="save-smiles">Save SMILES</option>
                    <option value="save-csv">Save CSV</option>
                </select>
                <div class="m2g-icon">
                    <svg width="20" height="20" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg">
                        <path d="M11.5 4C11.5 4.82843 10.8284 5.5 10 5.5C9.17157 5.5 8.5 4.82843 8.5 4C8.5 3.17157 9.17157 2.5 10 2.5C10.8284 2.5 11.5 3.17157 11.5 4ZM11.5 10C11.5 10.8284 10.8284 11.5 10 11.5C9.17157 11.5 8.5 10.8284 8.5 10C8.5 9.17157 9.17157 8.5 10 8.5C10.8284 8.5 11.5 9.17157 11.5 10ZM10 17.5C10.8284 17.5 11.5 16.8284 11.5 16C11.5 15.1716 10.8284 14.5 10 14.5C9.17157 14.5 8.5 15.1716 8.5 16C8.5 16.8284 9.17157 17.5 10 17.5Z"/>
                    </svg>
                </div>
            </div>
        </div>
    </div>

    <!-- Grid -->
    <div class="m2g-list"></div>
</div>`)
    initStyling(container, css)
    let content = chunks.join("\n")
    container.innerHTML = customHeader ? customHeader + content : content
    return container
}
