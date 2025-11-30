import { JSMol, RDKitModule } from "@rdkit/rdkit"
import { MolGrid } from "../molgrid"
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
    el: HTMLElement,
    molgrid: MolGrid,
    options: SmartsOptions,
    smartsMatches: SmartsMatches
) {
    function search(smilesCol: string): void {
        var query = $<HTMLInputElement>(".m2g-searchbar", el).elements[0].value
        if (typeof query !== "string") {
            return
        }
        // a bit dodgy but we can't use the async loader here
        // @ts-expect-error
        const RDKit: RDKitModule = window.__mol2gridRDKitModule__

        const qmol = RDKit.get_qmol(query)
        if (!qmol) {
            return
        }
        if (qmol.is_valid()) {
            molgrid.listObj.items.forEach((item: any) => {
                const smiles = item.values()[smilesCol]
                const mol = RDKit.get_mol(smiles, `{"removeHs": ${options.removeHs}}`)
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
                                var highlights = <QueryResult>{
                                    atoms: [],
                                    bonds: [],
                                }
                                results.forEach(function (match) {
                                    highlights["atoms"].push(...match.atoms)
                                    highlights["bonds"].push(...match.bonds)
                                })
                            }
                            const index: number = item.values()["mols2grid-id"]
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

    // wrapper for list.js
    return (_: string, columns: Array<any>) => {
        search(columns[0])
    }
}
