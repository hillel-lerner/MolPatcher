# MolPatcher: Protein-Ligand Attachment Tool

MolPatcher is a Python-based computational chemistry tool designed to automate the geometric alignment and topological merging of molecular patches onto target protein residues. 

## Core Architecture

The package consists of Python scripts designed to parse standard topologies, calculate geometric superpositions, and splice physical and topological arrays without relying on external visualization GUIs.

* **`pdb_io.py` & `topology_io.py`**: Handles the parsing and writing of PDB coordinate files and GROMACS `.itp` topology files, converting raw text into manageable Python objects.
* **`mol_record.py`**: Contains the dataclasses (e.g., `Mol`, `PdbRecord`, `ItpAtom`) that standardize attributes across the package.
* **`geometry.py`**: Houses the `PatchAligner` and `MolGraph` classes. Calculates the optimal rotation and translation matrices to superimpose anchor atoms using Kabsch-style vector alignment, and builds network connectivity matrices for distance-based bond mapping.
* **`stitcher.py`**: Executes the topological surgery. Analyzes the `MolGraph` to identify displaced atoms (e.g., leaving groups/hydrogens), splices the aligned patch coordinates into the base molecule, writes the new junction bonds, balances integer electrostatics, and renames the patched residue to `LYX`.
* **`topology_tools.py`**: Synchronizes and re-indexes the new global topology to ensure contiguous atom numbering for GROMACS compatibility.
* **`utilities.py`**: Provides mathematical helper functions for vector math, spatial distance calculations, and dynamic MD simulation box sizing.
* **`combine_ff.py`**: A standalone utility to merge base forcefields, modified parameters, and custom junction dihedrals into a single, deduplicated master forcefield.

## Requirements
* Python 3.11
* `numpy` (for matrix operations)
* `scipy` (for spatial alignments)
* `networkx` (for connectivity matrices)
---


## Installation & Usage

MolPatcher uses a conda environment to manage its dependencies and install its command-line tool.

### Setup (One Time Only)

MolPatcher uses a conda environment to manage its dependencies. 

1. Clone or download this repository to your local machine.
2. Navigate to the project directory in your terminal and create the environment. This will automatically install the `molpatcher` terminal command.
```bash
conda env create -f environment.yml
```
3. Activate the environment:
```bash
conda activate patcher
```
---

### Usage

Once installed, you can execute the pipeline from anywhere using the `molpatcher` command. Just activate the environment before running your command.

**Basic Command Structure:**
```bash
conda activate patcher

molpatcher --pdb [INPUT_PDB] --itp [INPUT_ITP] --res [TARGET_RESIDUE_ID] --chain [TARGET_CHAIN]
```

**Optional Arguments:**
* `-o`, `--outdir`: Specify a custom directory to save the patched outputs. Defaults to your current working directory. It is not recommended to run this package inside of the MolPatcher directory
* `--ff`, `--forcefield`: Provide the path to a master forcefield directory. If provided, the pipeline will copy this directory into your output folder.

**Example:**
```bash
molpatcher --pdb pdbs/fab_base.pdb --itp topology_files/PROB.itp --res 42 --chain B -o /home/user/simulations/ --ff /home/user/master_toppar/
```

## [!WARNING] 

### Important Warning: Geometric Verification

**Full steric clash detection is not yet implemented in this version of MolPatcher.** While the `PatchAligner` mathematically superimposes the anchor atoms to optimize the junction geometry, it does not evaluate the surrounding spatial environment for collisions between the newly attached patch and the rest of the protein backbone or adjacent side chains.

**Mandatory Post-Processing:** You must manually open the resulting PDB file in a molecular visualization software (such as PyMOL, Jmol, VMD, or Chimera) to visually verify the patch geometry and ensure there are no severe steric clashes before proceeding to molecular dynamics simulations.

