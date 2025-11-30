import initRDKitModule from "@rdkit/rdkit"
import type { RDKitModule } from "@rdkit/rdkit"
import { sleep } from "../utils"

export async function loadRDKit(initialize: boolean = false): Promise<RDKitModule> {
    if (window) {
        if (initialize) {
            console.log("Initializing RDKit")
            // @ts-expect-error
            let RDKit: RDKitModule = await initRDKitModule({
                locateFile: () =>
                    "https://unpkg.com/@rdkit/rdkit@2025.3.4-1.0.0/dist/RDKit_minimal.wasm",
            })
            console.log("RDKit version:", RDKit.version())
            // @ts-expect-error
            window.__mol2gridRDKitModule__ = RDKit
            return RDKit
        }
        // @ts-expect-error
        while (typeof window.__mol2gridRDKitModule__ === "undefined") {
            await sleep(25)
        }
        // @ts-expect-error
        return window.__mol2gridRDKitModule__
    }
    console.log(
        "No accessible window to cache RDKit module into, this will probably cause slowdowns."
    )
    // @ts-expect-error
    return await initRDKitModule({
        locateFile: () =>
            "https://unpkg.com/@rdkit/rdkit@2025.3.4-1.0.0/dist/RDKit_minimal.wasm",
    })
}
