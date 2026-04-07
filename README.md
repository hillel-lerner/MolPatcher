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

**Command-Line Arguments**

| Argument | Shorthand | Description | Required |
| :--- | :--- | :--- | :--- |
| `-pdb`, `--pdb` | | Path to the input PDB coordinate file. | **Yes** |
| `-itp`, `--itp` | | Path to the input GROMACS `.itp` topology file. | **Yes** |
| `-res`, `--res` | | The target residue ID (sequential integer). | **Yes** |
| `-chain`, `--chain` | `-c` | The target chain ID (e.g., A, B, C). | No (Default: `" "`) |
| `-outdir`, `--outdir` | `-o` | Custom directory for outputs. | No (Default: `./`) |
| `-ff`, `--ff` | | **Toggle:** Copy `forcefield_master.itp` to output. | No (Flag) |

**Basic Command Structure:**
```bash
molpatcher -pdb [INPUT_PDB] -itp [INPUT_ITP] -res [TARGET_RESIDUE_ID] -chain [TARGET_CHAIN]
```

**Example Commands**
```bash
molpatcher -pdb protein.pdb -itp protein.itp -res 145 -chain B -o ./patched_results -ff
```

---

## Forcefield Staging (`-ff`)

The `-ff` flag is a project-specific utility designed to streamline workflow for the **6oge** protein system. 

* **Behavior**: When enabled, MolPatcher automatically locates `forcefields/forcefield_master.itp` within the project root and copies it directly into your specific output directory alongside your patched PDB and ITP files.
* **Note**: This master file is pre-configured for specific parameters. For other protein systems or custom parameters, you can use the provided `combine_ff.py` tool in the `forcefields/` directory to generate a compatible master file before running the patcher.

---

## Automated Conformational Optimization

MolPatcher includes a built-in **RotamerSweeper** that automatically resolves steric clashes introduced by the patch. If the initial alignment results in a collision, the tool executes three tiers of search:

1.  **Tier 1 (Canonical)**: Tests the most common staggered rotameric states for the target residue.
2.  **Tier 2 (Systematic Wiggle)**: Performs a fine-grained 30° systematic sweep of the $\chi_4$ dihedral.
3.  **Tier 3 (Patch Twist)**: If there are still steric clashes, the tool performs 15° rotations around the new junction bond to clear local environment conflicts.

The tool provides real-time terminal feedback via an animation spinner and progress indicators during the $O(N^2)$ graph-building phase.

> **Note 1**: While MolPatcher optimizes for sterics, it is still recommended to visually inspect the resulting PDB in PyMOL or VMD before beginning long-production MD simulations.

> **Note 2**: MolPatcher automatically checks the residue numbering and lableing on your input files. It will renumber residues sequentially starting from 1 and clean up `.pdb` and `.itp` files to ensure compatibility with GROMACS processing tools. Always verify that you are modifying the correct residue when your input files contain non-integer residue numbers (e.g., 94B).