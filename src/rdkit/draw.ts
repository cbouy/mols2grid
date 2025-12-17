import type { JSMol } from "@rdkit/rdkit"
import type { SmartsMatches } from "./smarts"
import { $ } from "../query"
import { RDKit } from "./loader"

export interface DrawOptions {
    width: number
    height: number
}

export interface MolOptions {
    removeHs: boolean
    preferCoordGen: boolean
}

export function getEmptySvg(drawOptions: DrawOptions): string {
    return `<svg width="${drawOptions.width}" height="${drawOptions.height}" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 ${drawOptions.width} ${drawOptions.height}"></svg>`
}

// Generate images for the currently displayed molecules.
export async function drawMol(
    smiles: string | null,
    index: number,
    templateMol: JSMol | null | undefined,
    drawOptions: DrawOptions,
    molOptions: MolOptions,
    smartsMatches: SmartsMatches,
    placeholder: string
): Promise<string> {
    if (!smiles) {
        return placeholder
    }
    var mol: JSMol | null
    mol = RDKit.get_mol(smiles, `{"removeHs": ${molOptions.removeHs}}`)
    if (!mol || !mol.is_valid()) {
        return placeholder
    }
    var highlights = smartsMatches.get(index)
    if (highlights && templateMol) {
        var details = Object.assign({}, drawOptions, highlights)
        var jsonDetails = JSON.stringify(details)
        mol.generate_aligned_coords(
            templateMol,
            `{"useCoordGen": ${molOptions.preferCoordGen}}`
        )
    } else {
        var jsonDetails = JSON.stringify(drawOptions)
    }
    const svg = mol.get_svg_with_highlights(jsonDetails)
    mol?.delete()
    return svg
}

// Update images when the list is updated.
export function initMolDrawing(
    el: HTMLElement,
    smilesCol: string,
    drawOptions: DrawOptions,
    molOptions: MolOptions,
    smartsMatches: SmartsMatches,
    placeholder: string
) {
    var query = $<HTMLInputElement>(".m2g-searchbar", el).elements[0].value
    var templateMol: JSMol | null | undefined = null
    if (!query || typeof query !== "string") {
        smartsMatches.clear()
    } else {
        templateMol = RDKit.get_qmol(query)
        if (templateMol && templateMol.is_valid()) {
            templateMol.set_new_coords(molOptions.preferCoordGen)
        } else {
            templateMol = null
            smartsMatches.clear()
        }
    }
    $(".m2g-cell:not(.m2g-phantom)", el).each(cell => {
        const $t = $(cell)
        const imgEl = $t.find(".data-img").elements[0]
        imgEl.innerHTML = placeholder
        const smiles = $t.find(`.data-${smilesCol}`).index(0).text
        const index = parseInt(<string>cell.getAttribute("data-mols2grid-id"))
        drawMol(
            smiles,
            index,
            templateMol,
            drawOptions,
            molOptions,
            smartsMatches,
            placeholder
        ).then(svg => {
            imgEl.innerHTML = svg
        })
    })
    if (templateMol) {
        ;(<JSMol>templateMol).delete()
    }
}
