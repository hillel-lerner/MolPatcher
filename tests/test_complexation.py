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
complex_pdb = os.path.join(project_root, "pdbs", "nipah_m9_optimized_67.pdb")
complex_top = os.path.join(project_root, "topology_files", "patched_asn_67.itp")

parser = PdbParser()
scr_headers, scr_records, scr_name = parser.read_file(scr_pdb)
complex_headers, complex_records, complex_name = parser.read_file(complex_pdb)

scr_mol = Mol(name=scr_name, records=scr_records)
scr_mol.load_itp(scr_top)
complex_mol = Mol(name=complex_name, records=complex_records)
complex_mol.load_itp(complex_top)

combined_records = complex_mol.records + scr_mol.records
combined_atoms = complex_mol.atoms + scr_mol.atoms
combined_bonds = complex_mol.bonds + scr_mol.bonds

combined_mol = Mol(
    "Combined", records=combined_records, atoms=combined_atoms, bonds=combined_bonds
)
combined_graph = MolGraph(combined_mol)

dummy_config = {"patch_merged_atom": "C7", "docking": {"max_dist_nm": 1.5}}

conductor = SweepConductor(
    mol=combined_mol,
    graph=combined_graph,
    target_resid=68,
    config=dummy_config,
    chain_id="A",
)
glycan_target_atoms = [
    r for r in complex_mol.records if "MAN" in r.res_name or "GLC" in r.res_name
]

print(f"Debug: Isolated {len(glycan_target_atoms)} glycan atoms for centering")

success = conductor.rotate_and_dock(scr_mol=scr_mol, receptor_atoms=glycan_target_atoms)
outfile = os.path.join(project_root, "pdbs", "scr008_complex_docked.pdb")
if success:
    builder = PdbBuilder(outfile, combined_records, complex_headers)
    print(f"Outfile: {outfile}")
    builder.write_pdb()
