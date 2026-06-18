# MolPatcher: Protein-Ligand Attachment Tool

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

MolPatcher is a Python-based computational chemistry tool designed to automate the geometric alignment, topological merging, and conformational optimization of molecular patches onto target protein residues.

It handles the heavy lifting of bridging complex topologies (such as attaching large, branched N-glycans or flexible custom small molecules) by calculating precise spatial superpositions, resolving steric clashes, and extracting cross-boundary forcefield parameters without relying on external visualization GUIs.

## Core Features

* **Automated Topology Stitching:** Splices GROMACS `.itp` files, re-indexes atomic numbering, writes junction bonds, and balances local electrostatics.
* **Steric Resolution:** Features a built-in `RotamerSweeper` that uses a localized graph-traversal search to identify van der Waals overlaps. It sweeps sidechain dihedral angles ( $\chi_1$, $\chi_2$, etc. ) to isolate an optimal pose with minimal clashes.
* **Dynamic Forcefield Extraction:** Automatically parses raw CHARMM `.prm` and `.str` databases to extract the missing physical parameters (force constants, equilibrium distances) for newly formed covalent junctions and writes them to a GROMACS-ready `#include` file.
* **Configuration-First Design:** Chemistry rules are defined in internal templates, meaning daily runs require only a simple, highly readable JSON input file.

---
## Installation

MolPatcher is structured for modern Python environments and can be installed directly from source.

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/hillel-lerner/MolPatcher.git](https://github.com/hillel-lerner/MolPatcher.git)
   cd MolPatcher
   ```

2. **Create and activate the environment:**

    ```bash
    conda env create -f environment.yml
    conda activate mol_patcher
    ```
3. **Install the package locally:**

    ```bash
    pip install -e .
    ```

*(Note: The `-e` flag installs the package in editable mode, allowing you to modify the source code in `src/` without needing to reinstall.)*

---

## Quickstart: PFP-Lysine Patching
To run a patch, simply initialize a run file, map your file paths, and execute. You can find the required input files for this tutorial in the `examples/` directory.

**1. Initialize the Template** Generate a blank run file for the PFP-Lysine system in your current directory:

```bash
molpatcher --init pfp_lys

```

**2. Configure the Run** Open the generated `run.json` and map the file paths to point to the examples directory:

```json
{
    "_comment": "Mutating Lysine 188 on Chain A to PFP.",
    "template": "pfp_lys",
    "base": {
        "pdb": "examples/protein.pdb",
        "itp": "examples/PROB.itp",
        "resid": 188,
        "chain": "B"
    },
    "patch": {
        "pdb": "examples/pfp_patch_new.pdb",
        "itp": "examples/pfp_patch_new.itp"
    },
    "forcefields": []
}

```

**3. Execute the Patch** Run the engine. MolPatcher will automatically resolve the topological junction, optimize the sidechain rotamers to prevent steric clashes, and output the stitched files:

```bash
molpatcher run.json

```

---

## Output Structure & Logging

MolPatcher strictly separates your raw inputs from its generated topologies to preserve reproducibility. It outputs a dedicated directory containing the stitched coordinates, the merged topology, and a master `.log` file.

```text
patched_lyx_188_PROB/
├── toppar/
│   ├── patched_lyx_188.itp
│   └── junction_ff_lyx_188.itp
├── patched_lyx_188.pdb
└── patched_lyx_188.log

```

### The Execution Log

The generated `.log` file provides exact instructions for integrating the new system into your MD workflow, including Periodic Boundary Condition (PBC) sizing and the explicit parameters extracted from the forcefield databases.

```text
========== MOLPATCHER EXECUTION LOG ==========
Target : LYS 188 (Chain A)
Patch  : pfp_patch_new.pdb
==============================================

--- GROMACS TOPOLOGY INSTRUCTIONS ---
1. Add the following to your master .top file (after the forcefield #includes):
   #include "patched_lyx_188.itp"
   #include "junction_ff_lyx_188.itp"

2. Add the following to the [ molecules ] directive:
   PROA            1

--- PERIODIC BOUNDARY CONDITIONS (PBC) ---
Molecule Max Diagonal : 12.450 nm
Applied PBC Buffer    : 4.108 nm
Optimal Cubic Box     : 16.558 nm

--- JUNCTION FORCEFIELD PARAMETERS ---
[ bondtypes ]
CC3162   NC2D1    1   0.14300    267776.0

```

---

## Advanced Usage: Branched N-Glycans & Pre-Docked Complexes

For large patches, such as a large N-glycan with a pre-bound Synthetic Carbohydrate Receptor (SCR), MolPatcher can handle multi-topology integration.

You can initialize a glycan-specific run file:

```bash
molpatcher --init nglycan_asn

```

The resulting `run.json` supports additional keys for secondary topologies (`scr_itp`) and raw CHARMM forcefield files (`forcefields`):

```json
{
    "_comment": "Patching a Glycan/SCR complex.",
    "template": "nglycan_asn",
    "base": {
        "pdb": "protein.pdb",
        "itp": "PROB.itp",
        "resid": 67,
        "chain": "B"
    },
    "patch": {
        "pdb": "glycan_rec.pdb",
        "itp": "CARB.itp"
    },
    "scr_itp": "LIG.itp",
    "forcefields": [
        "/path/to/charmm36/par_all36_carb.prm",
        "/path/to/charmm36/toppar_all36_carb_glycopeptide.str"
    ]
}

```

*Note: For extremely crowded pockets, the rigid-body `RotamerSweeper` will calculate and apply the lowest possible penalty pose. It is recommended to perform a standard steepest-descent energy minimization in GROMACS to allow the surrounding protein sidechains to relax around the new attachment before beginning production MD.*

---

## License

MolPatcher is distributed under the GNU General Public License v3.0 (GPLv3). See `LICENSE` for more information.
```

```
## Development & Contributing

If you want to modify the source code or contribute to MolPatcher, it is recommended to set up an isolated development environment using Conda.

```bash
git clone [https://github.com/hillel-lerner/MolPatcher.git](https://github.com/hillel-lerner/MolPatcher.git)
cd MolPatcher
conda env create -f environment.yml
conda activate mol_patcher
pip install -e .
