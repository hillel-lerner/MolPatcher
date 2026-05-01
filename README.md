# MolPatcher: Protein-Ligand Attachment Tool

MolPatcher is a Python-based computational chemistry tool designed to automate the geometric alignment, topological merging, and conformational optimization of molecular patches onto target protein residues.

## Core Architecture

The package consists of modular components designed to parse standard topologies, calculate geometric superpositions, and resolve steric clashes without relying on external visualization GUIs.

* **`pdb_io.py` & `topology_io.py`**: Handles parsing and writing of PDB coordinate files and GROMACS `.itp` files.
* **`mol_record.py`**: Standardizes molecular data into Python dataclasses (`Mol`, `PdbRecord`, `ItpAtom`).
* **`geometry.py`**: The engine of the package.
  * **`PatchAligner`**: Uses Kabsch-style vector alignment for optimal 3D superposition.
  * **`MolGraph`**: Builds a NetworkX-based connectivity graph for distance-based bond mapping.
  * **`StericChecker`**: Evaluates spatial environments using a localized graph-traversal search to identify overlaps based on VdW radii (Charry, Tkatchenko, *J. Chem. Theory Comput.* **2024**, 20, 7469–7478.)
  * **`RotamerSweeper`**: Executes a three-tiered conformational search to find clash-free states.
* **`stitcher.py`**: Executes the topological surgery. Splices aligned patch coordinates, writes junction bonds, and balances electrostatics.
* **`topology_tools.py`**: Reindexes the topology modified in `stitcher.py`.
* **`utilities.py`**: Helper functions for vector math, dihedral rotations, and dynamic MD box sizing.
* **`main.py`**: Excecutes the MolPatcher pipeline.

## Standalone Utility Functions

* **`combine_ff.py`**: Use to build combined master forcefield files with custom parameters.

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
| `-ff` | **Toggle:** Copy `forcefield_master.itp` to output. | No | Flag |

### Basic Command Structure

```bash
molpatcher -config [JSON] -base [BASE_PDB] [BASE_ITP] [BASE_RESIDUE_ID] [BASE_CHAIN] -patch [PATCH_PDB] [PATCH_ITP]

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
           -patch m9_glycan.pdb m9_glycan.itp 1 A
```

> **Note**: MolPatcher automatically checks the residue numbering and labeling on your input files. It will renumber residues sequentially starting from 1 and clean up `.pdb` and `.itp` files to ensure compatibility with GROMACS processing tools. Always verify that you are modifying the correct residue when your input files contain non-integer residue numbers (e.g., 94B).

---

## JSON Configuration Templates

The specific method used to attach the base and patch molecule is determined by the JSON file specified. Each configuration uses a specific geometry engine.

**1. Flexible/Merged Linkages**
The `chi_definitions` here are used to define the dihedral angles that the `RotamerSweeper` rotates

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
    }
}
```

## Output Structure & Logging

MolPatcher preserves the original input files for reproducibility and generates a master log with topology `#include` instructions.

```text
patched_asx_67_PROB/
├── infiles/
│   ├── base_input.pdb
│   ├── base_input.itp
│   ├── patch_input.pdb
│   └── patch_input.itp
├── patched_asx_67.pdb
├── patched_asx_67.itp
└── patcher.log
```

**Example `patcher.log` Output:**
This file provides the exact GROMACS `.top` syntax.

```text
========== MOLPATCHER EXECUTION LOG ==========
Target : ASN 67 (Chain B)
Patch  : m9_glycan.pdb
==============================================

--- GROMACS TOPOLOGY INSTRUCTIONS ---
1. Add the following to your master .top file:
   #include "patched_asx_67.itp"

2. Add the following to the [ molecules ] directive:
   PROB            1

--- PERIODIC BOUNDARY CONDITIONS (PBC) ---
Molecule Max Diagonal : 12.450 nm
Applied PBC Buffer    : 4.108 nm
Optimal Cubic Box     : 16.558 nm
```

## Automated Conformational Optimization

MolPatcher includes a built-in **RotamerSweeper** that automatically resolves steric clashes introduced by the patch. If the initial alignment results in a collision, the tool executes three tiers of search:

1. **Tier 1 (Canonical)**: Tests the most common staggered rotameric states for the target residue.
2. **Tier 2 (Systematic Wiggle)**: Performs a fine-grained 30° systematic sweep of the $\chi_4$ dihedral.
3. **Tier 3 (Patch Twist)**: If there are still steric clashes, the tool performs 15° rotations around the new junction bond to clear local environment conflicts.

The tool provides real-time terminal feedback via an animation spinner and progress indicators during the $O(N^2)$ graph-building phase.

> **Note**: While MolPatcher optimizes for sterics, it is still recommended to visually inspect the resulting PDB in PyMOL or VMD before beginning long-production MD simulations.

> **Note**: The `RotamerSweeper` is intelligently bypassed if the input JSON configures a rigid geometric template (e.g., specific N-glycan torsions). This ensures explicitly calculated junction geometries are not overwriten

---

## Forcefield Staging (`-ff`)

The `-ff` flag is a project-specific utility designed to streamline workflow for the **6oge** protein system.

* **Behavior**: When enabled, MolPatcher automatically locates `forcefields/forcefield_master.itp` within the project root and copies it directly into your specific output directory alongside your patched PDB and ITP files.
* **Note**: This master file is pre-configured for specific parameters. For other protein systems or custom parameters, you can use the provided `combine_ff.py` tool in the `forcefields/` directory to generate a compatible master file before running the patcher.
