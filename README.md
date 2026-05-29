# MolPatcher: Protein-Ligand Attachment Tool

MolPatcher is a Python-based computational chemistry tool designed to automate the geometric alignment, topological merging, and conformational optimization of molecular patches onto target protein residues.

## Core Architecture

The package consists of modular components designed to parse standard topologies, calculate geometric superpositions, resolve steric clashes, and dynamically generate junction forcefield parameters without relying on external visualization GUIs.

* **`pdb_io.py` & `topology_io.py`**: Handles parsing and writing of PDB coordinate files and GROMACS `.itp` files.
* **`mol_record.py`**: Standardizes molecular data into Python dataclasses (`Mol`, `PdbRecord`, `ItpAtom`).
* **`geometry.py`**: The spatial engine of the package.
  * **`PatchAligner`**: Uses Kabsch-style vector alignment for optimal 3D superposition and dihedral manipulation.
  * **`MolGraph`**: Builds a NetworkX-based connectivity graph for distance-based bond mapping.
  * **`StericChecker`**: Evaluates spatial environments using a localized graph-traversal search to identify overlaps based on VdW radii (Charry, Tkatchenko, *J. Chem. Theory Comput.* **2024**, 20, 7469–7478.)
* **`sweeper.py`**: Orchestrates conformational optimization via `RotamerSweeper` and `SweepConductor` to isolate clash-free states.
* **`stitcher.py`**: Executes the topological surgery. Splices aligned patch coordinates, writes junction bonds, balances electrostatics, and maps junction interactions for parameter extraction.
* **`topology_tools.py`**: Re-indexes the merged topology output by `stitcher.py` to ensure sequential GROMACS numbering.
* **`combine_ff.py`**: Reads raw CHARMM `.prm` and `.str`  files to extract and translate forcefield parameters for the newly created junction.
* **`utilities.py`**: Helper functions for vector math, dihedral rotations, and dynamic MD box sizing.
* **`main.py`**: CLI entry point that executes the MolPatcher pipeline.

---

## Installation

MolPatcher uses a conda environment to manage dependencies and install the command-line tool.

1. **Clone the repository:**

    ```bash
    git clone https://github.com/hillel-lerner/MolPatcher.git
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

---

## Usage

Once installed, you can execute the pipeline using the `molpatcher` command.

### Command-Line Arguments

| Argument | Description | Required | Format |
| :--- | :--- | :--- | :--- |
| `-base` | Base protein inputs. | **Yes** | `[PDB] [ITP] [RESID] [CHAIN]` |
| `-patch` | Patch molecule inputs. | **Yes** | `[PDB] [ITP] [Optional RESID] [Optional CHAIN]` |
| `-config` | Path to the JSON configuration template. | **Yes** | String (filepath) |
| `-o`, `--outdir` | Custom directory for outputs. | No | String (filepath, Default: `./`) |
| `-ff`, `--forcefield` | Additional forcefield files to copy into the output `toppar/` directory. | No | `nargs=*` (multiple filepaths separated by spaces)|

### Basic Command Structure

```bash
molpatcher -config [JSON] -base [BASE_PDB] [BASE_ITP] [BASE_RESIDUE_ID] [BASE_CHAIN] -patch [PATCH_PDB] [PATCH_ITP] -ff [FORCEFIELD_1] [FORCEFIELD_2]

```

### Example Commands

**1. Flexible Small Molecule (PFP-Lysine)**

```bash
molpatcher -config pfp_lysine.json \
           -base protein.pdb protein_PROB.itp 188 B \
           -patch pfp_ligand.pdb pfp_ligand.itp
```

**1. Rigid Molecule Appending (N-Glycan)**

```bash
molpatcher -config glycan_aspargine.json \
           -base protein.pdb protein_PROB.itp 67 B \
           -patch m9_glycan.pdb m9_glycan.itp 1 A \
           -ff forcefield.itp ligand_forcfield.itp
```

> **Note**: MolPatcher automatically checks the residue numbering and labeling on your input files. It will renumber residues sequentially starting from 1 and clean up `.pdb` and `.itp` files to ensure compatibility with GROMACS processing tools. Always verify that you are modifying the correct residue when your input files contain non-integer residue numbers (e.g., 94B).

---

## JSON Configuration Templates

The specific method used to attach the base and patch molecule is determined by the JSON file specified. Each configuration uses a specific geometry engine.

**1. Flexible/Merged Linkages**
The `chi_definitions` here are used to define the dihedral angles that the `RotamerSweeper` rotates to resolve steric clashes. For completely custom interactions, parameters can be explicitly defined in an ff_parameters block..

```json
{
    "system_name": "pfp_lysine",
    "target_res_name": "LYS",
    "new_res_name": "LYX",
    "base_bridge": "NZ",
    "patch_bridge": "C7",
    "patch_merged_atom": "N",
    "chi_definitions": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"]
    ]
}
```

**2. Rigid Geometry Linkages**
This configuration uses a `"geometry"` block to enforce strict bond angles/dihedrals, bypassing the RotamerSweeper.

```json
{
    "system_name": "glycan_aspargine",
    "target_res_name": "ASN",
    "new_res_name": "ASX",
    "anchors": {
        "base": {"res_name": "ASN", "parent": "CG", "anchor": "ND2", "ref": "CB", "carbonyl": "OD1"},
        "patch": {"res_name": "BGL", "anchor": "C1", "leaving": "O1", "child": "O5"}
    },
    "geometry": {
        "bond_length": 1.43,
        "angle_target": 120.0,
        "omega_target": 180.0,
        "phi_target": 180.0
    },
  "charmm_forcefield_files": [
        "~/charmm36/par_all36_carb.prm",
        "~/charmm36/toppar_all36_carb_glycopeptide.str"
    ]
}
```

## Output Structure & Logging

MolPatcher separates your raw inputs from its generated topologies. It preserves the original input files for reproducibility and generates a master .log file containing explicit topology `#include` instructions and extracted forcefield parameters.

```text
patched_asx_67_PROB/
├── infiles/
│   ├── base.pdb
│   ├── base.itp
│   ├── patch.pdb
│   └── patch.itp
├── toppar/
│   ├── patched_asx_67.itp
│   ├── junction_ff_asx_67.itp
│   └── [any files passed via -ff]
├── patched_asx_67.pdb
└── patched_asx_67.log
```

**Example `patcher.log` Output:**
This file provides the exact GROMACS `.top` syntax.

```text
========== MOLPATCHER EXECUTION LOG ==========
Target : ASN 67 (Chain B)
Patch  : m9_glycan.pdb
==============================================

--- GROMACS TOPOLOGY INSTRUCTIONS ---
1. Add the following to your master .top file (after the forcefield #includes):
   #include "patched_asx_67.itp"

2. Add the following to the [ molecules ] directive at the bottom of your .top file:
   PROB            1

--- PERIODIC BOUNDARY CONDITIONS (PBC) ---
Molecule Max Diagonal : 12.450 nm
Applied PBC Buffer    : 4.108 nm
Optimal Cubic Box     : 16.558 nm

--- JUNCTION FORCEFIELD PARAMETERS ADDED ---
[ bondtypes ]
CC3162   NC2D1    1   0.14300    267776.0

[ angletypes ]
CC3162   NC2D1    CC2O1    5  120.00000   418.4   0.00000  0.0
```

## Automated Conformational Optimization

MolPatcher includes a built-in **RotamerSweeper** that automatically resolves steric clashes introduced by the patch. If the initial alignment results in a collision, the tool evaluates a canonical set of staggered rotameric states for the target residue to isolate a clash-free geometry.the tool executes three tiers of search:

The tool provides real-time terminal feedback via an animation spinner and progress indicators during the $O(N^2)$ graph-building phase.

> **Note**: The `RotamerSweeper` is intelligently bypassed if the input JSON configures a rigid geometric template (e.g., specific N-glycan torsions). This ensures explicitly calculated junction geometries are not overwritten. While MolPatcher optimizes for sterics, it is still recommended to visually inspect the resulting PDB in PyMOL or VMD before beginning long-production MD simulations.

---

## Dynamic Forcefield Extraction

When a new covalent bond is formed between two distinct forcefields (e.g., a CHARMM36 protein and a CHARMM36 carbohydrate), the standard GROMACS `grompp` compiler will fail due to missing cross-boundary parameters.

MolPatcher natively solves this by generating a `junction_ff_[res].itp` file.

* By providing a list of raw CHARMM `.prm` and `.str` files in your JSON config (`charmm_forcefield_files`), MolPatcher will automatically map the newly formed topology.
* It parses the raw CHARMM databases, extracts the required physical parameters (force constants, equilibrium distances), translates them into GROMACS units, and writes them into a formatted global `#include` file.
* Aliasing logic natively handles boundary crossings (e.g., treating a modified glycan amide carbon like a standard protein backbone carbon when looking up historical backbone parameters).
