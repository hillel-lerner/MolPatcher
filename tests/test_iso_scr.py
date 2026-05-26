import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.geometry import random_rot
from mol_patcher.sweeper import SCRSweeper
from mol_patcher.mol_record import Mol

scr_pdb = os.path.join(project_root, "pdbs", "SCR007.pdb")
scr_top = os.path.join(project_root, "topology_files", "SCR007.itp")
glycan_pdb = os.path.join(project_root, "pdbs", "m9.pdb")
glycan_top = os.path.join(project_root, "topology_files", "CARB-m9.itp")

parser = PdbParser()
scr_headers, scr_records, scr_name = parser.read_file(scr_pdb)
glycan_headers, glycan_records, glycan_name = parser.read_file(glycan_pdb)

scr_mol = Mol(
    name=scr_name, records=scr_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[]
)
scr_mol.load_itp(scr_top)
glycan_mol = Mol(
    name=glycan_name,
    records=glycan_records,
    atoms=[],
    bonds=[],
    pairs=[],
    angles=[],
    dihs=[],
)
glycan_mol.load_itp(glycan_top)

scr_sweep = SCRSweeper()
test = scr_sweep.dock_near(glycan_atoms=glycan_mol.records, scr_atoms=scr_mol.records)


outfile = os.path.join(project_root, "pdbs", "scr007_m9_docked.pdb")

# builder = PdbBuilder(outfile, )
# builder.write_pdb()
