import { RDKit } from "../initialize"
import { MolGrid } from "../molgrid"
import { JSMol } from "@rdkit/rdkit"
import { $ } from "../query"

export interface SmartsOptions {
    removeHs: boolean
    onTheFlyRendering: boolean
    substructHighlight: boolean
    singleHighlight: boolean
}

export interface QueryResult {
    atoms: number[]
    bonds: number[]
}

export type SmartsMatches = Map<number, QueryResult>

export function smartsSearchFactory(
    molgrid: MolGrid,
    options: SmartsOptions,
    smartsMatches: SmartsMatches
) {
    return (_: string, columns: Array<any>) => {
        var smilesCol: string = columns[0]
        var query = $<HTMLInputElement>(`#${molgrid.listObj.listContainer.id} .m2g-searchbar`).elements[0].value
        if (typeof query !== "string") {
            return
        }
        var qmol = RDKit?.get_qmol(query)
        if (!qmol) {
            return
        }
        if (qmol.is_valid()) {
            molgrid.listObj.items.forEach((item: any) => {
                var smiles = item.values()[smilesCol]
                var mol = RDKit?.get_mol(smiles, `{"removeHs": ${options.removeHs}}`)
                if (!mol) {
                    item.found = false
                    return
                }
                if (mol.is_valid()) {
                    var jsonResults = mol.get_substruct_matches(<JSMol>qmol)
                    if (jsonResults === "{}") {
                        item.found = false
                    } else {
                        item.found = true
                        if (options.onTheFlyRendering && options.substructHighlight) {
                            let results: QueryResult[] = JSON.parse(jsonResults)
                            if (options.singleHighlight) {
                                var highlights = results[0]
                            } else {
                                var highlights = <QueryResult>{ atoms: [], bonds: [] }
                                results.forEach(function (match) {
                                    highlights["atoms"].push(...match.atoms)
                                    highlights["bonds"].push(...match.bonds)
                                })
                            }
                            var index: number = item.values()["mols2grid-id"]
                            smartsMatches.set(index, highlights)
                        }
                    }
                } else {
                    item.found = false
                }
                mol.delete()
            })
        }
        qmol.delete()
    }
}
