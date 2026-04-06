# MolPatcher: Protein-Ligand Attachment Tool

MolPatcher is a Python-based computational chemistry tool designed to automate the geometric alignment, topological merging, and conformational optimization of molecular patches onto target protein residues.

## Core Architecture

The package consists of modular components designed to parse standard topologies, calculate geometric superpositions, and resolve steric clashes without relying on external visualization GUIs.

* **`pdb_io.py` & `topology_io.py`**: Handles parsing and writing of PDB coordinate files and GROMACS `.itp` files. 
* **`mol_record.py`**: Standardizes molecular data into Python dataclasses (`Mol`, `PdbRecord`, `ItpAtom`).
* **`geometry.py`**: The engine of the package.
    * **`PatchAligner`**: Uses Kabsch-style vector alignment for optimal 3D superposition.
    * **`MolGraph`**: Builds a NetworkX-based connectivity graph for distance-based bond mapping.
    * **`StericChecker`**: Evaluates spatial environments using a localized graph-traversal search to identify overlaps based on 2024 Charry/Tkatchenko VdW radii.
    * **`RotamerSweeper`**: Executes a three-tiered conformational search to find clash-free states.
* **`stitcher.py`**: Executes the topological surgery. Splices aligned patch coordinates, writes junction bonds, and balances electrostatics.
* **`topology_tools.py`**: Reindexes the topology modified in `stitcher.py`.
* **`utilities.py`**: Helper functions for vector math, dihedral rotations, and dynamic MD box sizing.
* **`main.py`**: Excecutes the MolPatcher pipeline.

## Standalone Utility Functions

* **`combine_ff.py`**: Use to build combined master forcefield files with custom parameters.
* **`fix_itp_resnr.py`**: Use to fix and reorder any topology files with non-integer residue numbers (e.g., 94B).

---

## Installation

MolPatcher uses a conda environment to manage dependencies and install the command-line tool.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/hillel-lerner/MolPatcher.git
    ```
2.  **Create and activate the environment:**
    ```bash
    conda env create -f environment.yml
    conda activate mol_patcher
    ```
---

## Usage

Once installed, you can execute the pipeline using the `molpatcher` command.

**Basic Command Structure:**
```bash
molpatcher --pdb [INPUT_PDB] --itp [INPUT_ITP] --res [TARGET_RESIDUE_ID] --chain [TARGET_CHAIN]
```

**Optional Arguments:**
* `-o`, `--outdir`: Specify a custom directory to save outputs. Defaults to the current working directory.

**Example:**
```bash
molpatcher --pdb protein.pdb --itp protein.itp --res 145 --chain B -o ./results/
```

---

## Automated Conformational Optimization

MolPatcher includes a built-in **RotamerSweeper** that automatically resolves steric clashes introduced by the patch. If the initial alignment results in a collision, the tool executes three tiers of search:

1.  **Tier 1 (Canonical)**: Tests the most common staggered rotameric states for the target residue.
2.  **Tier 2 (Systematic Wiggle)**: Performs a fine-grained 30° systematic sweep of the $\chi_4$ dihedral.
3.  **Tier 3 (Patch Twist)**: If the residue is "stuck," the tool performs 15° rotations around the new junction bond to clear local environment clashes.

The tool provides real-time terminal feedback via an animation spinner and progress indicators during the $O(N^2)$ graph-building phase.

> **Note**: While MolPatcher optimizes for sterics, it is still recommended to visually inspect the resulting PDB in PyMOL or VMD before beginning long-production MD simulations.