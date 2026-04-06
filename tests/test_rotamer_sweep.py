import os
import sys

# Ensure the project root is in the path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.geometry import PatchAligner, MolGraph, StericChecker, RotamerSweeper
from mol_patcher.mol_record import Mol

# --- 1. SETUP PARAMETERS ---
prot_path = os.path.join(project_root, 'pdbs', "PROBC_fixed.pdb")
top_path = os.path.join(project_root, 'topology_files', 'PROBC_fixed.itp')
target_res_seq = 145  
chain = 'B'
output_path = os.path.join(project_root, 'pdbs', 'PROBC_patched_optimized.pdb')

# --- 2. LOAD PROTEIN AND PATCH ---
parser = PdbParser()
headers, prot_records, t_name = parser.read_file(prot_path)
base_mol = Mol(name=t_name, records=prot_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
base_mol.load_itp(top_path)

loader = Patchloader()
patch_pdb_atoms = loader.get_pfp_pdb()
patch_mol = Mol("patch", patch_pdb_atoms, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
loader.get_pfp_itp(patch_mol)

# --- 3. LOCATE ANCHORS ---
print(f"Locating anchors for Lysine {target_res_seq}...")
chain_check = chain.strip()

target_anchors = []
for name in ["CE", "CD", "NZ"]:
    anchor = next(r for r in base_mol.records 
                if r.res_seq == target_res_seq 
                and r.name.strip() == name 
                and (not chain_check or r.chain.strip() == chain_check))
    target_anchors.append(anchor)

pfp_anchors = [
    next(r for r in patch_mol.records if r.name.strip() == "C10"),
    next(r for r in patch_mol.records if r.name.strip() == "C11"),
    next(r for r in patch_mol.records if r.name.strip() == "N")
]

# --- 4. ALIGN AND STITCH ---
print("Aligning and stitching patch...")
aligner = PatchAligner(patch_mol.records, pfp_anchors, target_anchors)
aligned_pfp_atoms = aligner.implement_align()

stitcher = Stitcher(base_mol, patch_mol, target_res_id=target_res_seq)
stitched_mol = stitcher.stitch_molecules(
    aligned_patch_atoms=aligned_pfp_atoms, 
    target_reference=target_anchors[-1], 
    target_anchors=target_anchors
)

# --- 5. OPTIMIZATION SETUP ---
graph = MolGraph(stitched_mol, distXX=1.90)
graph = MolGraph(stitched_mol)

backbone_names = ['N', 'H', 'HN', 'CA', 'HA', 'C', 'O']

moving_atoms = [a for a in stitched_mol.records 
                if a.res_seq == target_res_seq and a.name.strip() not in backbone_names]

static_atoms = [a for a in stitched_mol.records if a.res_seq != target_res_seq]

checker = StericChecker(graph, moving_atoms, static_atoms, tolerance=0.75)
sweeper = RotamerSweeper(stitched_mol, graph, target_res_seq)

# --- 6. EXECUTE SWEEP ---
print(f"Starting optimization for {len(moving_atoms)} moving atoms...")
success = sweeper.run_sweep(checker)

# --- 7. SAVE OUTPUT ---
if success:
    print(f"Success! Saving optimized PDB to {output_path}")
    builder = PdbBuilder(output_path, stitched_mol.records, headers)
    builder.write_pdb()
else:
    print("Optimization failed to find a clash-free conformation.")