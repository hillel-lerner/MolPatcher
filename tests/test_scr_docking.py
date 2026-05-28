import os
import sys

if os.path.dirname(os.path.dirname(os.path.abspath(__file__))) not in sys.path:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.geometry import MolGraph
from mol_patcher.sweeper import SweepConductor
from mol_patcher.mol_record import Mol

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

scr_pdb = os.path.join(project_root, "pdbs", "SCR007.pdb")
scr_top = os.path.join(project_root, "topology_files", "SCR007.itp")
glycan_pdb = os.path.join(project_root, "pdbs", "m9.pdb")
glycan_top = os.path.join(project_root, "topology_files", "CARB-m9.itp")

parser = PdbParser()
scr_headers, scr_records, scr_name = parser.read_file(scr_pdb)
glycan_headers, glycan_records, glycan_name = parser.read_file(glycan_pdb)

scr_mol = Mol(name=scr_name, records=scr_records)
scr_mol.load_itp(scr_top)
glycan_mol = Mol(name=glycan_name, records=glycan_records)
glycan_mol.load_itp(glycan_top)

combined_records = glycan_mol.records + scr_mol.records
combined_atoms = glycan_mol.atoms + scr_mol.atoms
combined_bonds = glycan_mol.bonds + scr_mol.bonds

complex_mol = Mol(
    "Complex", records=combined_records, atoms=combined_atoms, bonds=combined_bonds
)
complex_graph = MolGraph(complex_mol)

dummy_config = {"patch_merged_atom": "C7", "docking": {"max_dist_nm": 1.5}}

conductor = SweepConductor(
    mol=complex_mol,
    graph=complex_graph,
    target_resid=67,
    config=dummy_config,
    chain_id=" ",
)
success = conductor.rotate_and_dock(scr_mol=scr_mol, receptor_atoms=glycan_mol.records)
outfile = os.path.join(project_root, "pdbs", "scr007_m9_docked.pdb")
if success:
    builder = PdbBuilder(outfile, combined_records, glycan_headers)
    print(f"Outfile: {outfile}")
    builder.write_pdb()
