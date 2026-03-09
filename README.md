# MolPatcher: Protein-Ligand Attachment Tool

**MolPatcher** is a computational chemistry tool designed to precisely insert ligand "patches" into "base" protein residues at a specified anchor site while maintaining topology synchronization and geometric alignment.

## Key Features

* **Surgical Residue Insertion**: Rather than append atoms to the end of a file, MolPatcher inserts patch atoms immediately after the anchor, maintaining physical residue contiguity.
* **Synchronized Topology & PDB**: Simultaneously generates valid PDB and ITP files where atom serials, charge groups (`cgnr`), and all bonded interaction pointers (bonds, angles, dihedrals) are automatically re-indexed.
* **Precision Geometry**: Uses a 3-point Kabsch algorithm to orient patches with sub-angstrom precision relative to the target residue.
* **CHARMM-GUI Optimized**: Specifically designed for all-atom PDBs and naming conventions generated via **CHARMM-GUI**.
* **Safety Layers**: Features an anchor protection layer to prevent backbone deletion and a "fail-fast" mechanism for missing anchors.

---

## Installation & Usage

### Dependencies

* **Python 3.10+**
* **NumPy** & **SciPy**

### Execution

MolPatcher supports both CLI for production and interactive mode for debugging.

**CLI Mode:**

```bash
python main.py --pdb protein.pdb --itp protein.itp --res 188 --chain "A"

```

**Interactive Mode (VS Code):**
Run `main.py` directly; it defaults to pre-configured research parameters if no arguments are provided.

---

## Adding New Patch Types

To integrate a new chemical modification (e.g., a halogenated group), follow these steps:

### 1. File Preparation

Place the new `.pdb` (with 3 anchors) and `.itp` files in the `pdbs/` and `itps/` directories.

### 2. Configure Anchors

In `main.py`, define the **Mobile Anchors** (patch-side) and **Target Anchors** (protein-side) to lock the geometry.

### 3. Update the Call

Update `run_patch` to load your new ligand files:

```python
pfp_atoms = mol_stitcher.get_new_patch_pdb("new_patch.pdb")
pfp_mol = Mol(name="NEW", records=pfp_atoms)
pfp_mol.load_itp("new_patch.itp")

```

### 4. Set the Junction Bridge

In `mol_stitcher.py`, verify the junction bond (e.g., base N to patch C) and associated angles/dihedrals reflect the new chemistry.

---

## Project Roadmap

* **Steric Clash Resolution**: Implementing automatic ligand rotation to prevent steric clashed from occuring between the patch and the base molecule (protein) as a whole.
* **Advanced Deletion**: Implementing network-connectivity logic to prune atoms more reliably than distance-based methods.