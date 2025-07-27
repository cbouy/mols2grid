import { RDKit } from "../initialize"
import type { JSMol } from "@rdkit/rdkit"
import type { SmartsMatches } from "./smarts"

export interface DrawOptions {
    width: number
    height: number
}

export interface MolOptions {
    removeHs: boolean
    preferCoordGen: boolean
}

// Generate images for the currently displayed molecules.
export function drawMol(
    smiles: string,
    index: number,
    templateMol: JSMol | null | undefined,
    drawOptions: DrawOptions,
    molOptions: MolOptions,
    smartsMatches: SmartsMatches
): string {
    var mol = RDKit?.get_mol(smiles, `{"removeHs": ${molOptions.removeHs}}`)
    if (!mol || !mol.is_valid()) {
        var svg = `<svg width="${drawOptions.width}" height="${drawOptions.height}" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 ${drawOptions.width} ${drawOptions.height}"></svg>`
    } else {
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
        svg = mol.get_svg_with_highlights(jsonDetails)
    }
    mol?.delete()
    return svg
}

// Update images when the list is updated.
export function initMolDrawing(
    smilesCol: string,
    drawOptions: DrawOptions,
    molOptions: MolOptions,
    smartsMatches: SmartsMatches
) {
    var query = $("#mols2grid .m2g-searchbar").val()
    var templateMol: JSMol | null | undefined = null
    if (!query || typeof query !== "string") {
        smartsMatches.clear()
    } else {
        templateMol = RDKit?.get_qmol(query)
        if (templateMol && templateMol.is_valid()) {
            templateMol.set_new_coords(molOptions.preferCoordGen)
        }
    }
    $("#mols2grid .m2g-cell").each(function (_: number, el: HTMLElement) {
        var $t = $(el)
        var smiles = $t.children(`.data-${smilesCol}`).first().text()
        var index = parseInt(<string>el.getAttribute("data-mols2grid-id"))
        var svg = drawMol(
            smiles,
            index,
            templateMol,
            drawOptions,
            molOptions,
            smartsMatches
        )
        $t.children(".data-img").html(svg)
    })
    if (templateMol) {
        ;(<JSMol>templateMol).delete()
    }
}
