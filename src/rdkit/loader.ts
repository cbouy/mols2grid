import initRDKitModule from "@rdkit/rdkit"
import type { RDKitModule } from "@rdkit/rdkit"

console.log("Initializing RDKit")
// @ts-expect-error
export const RDKit: RDKitModule = await initRDKitModule({
    locateFile: () =>
        "https://unpkg.com/@rdkit/rdkit@2025.3.4-1.0.0/dist/RDKit_minimal.wasm",
})
console.log("RDKit version:", RDKit.version())
