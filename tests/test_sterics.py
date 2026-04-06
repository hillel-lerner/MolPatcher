import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

# Adjust these imports based on what file you saved StericChecker into
from mol_patcher.pdb_io import PdbParser
from mol_patcher.geometry import MolGraph
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import StericChecker 

pdb_dir = os.path.join(project_root, 'pdbs')
test_pdb = os.path.join(pdb_dir, "PROBC_fixed.pdb") 

# 1. Parse the file
parser = PdbParser()
header_lines, final_atoms, file_name = parser.read_file(test_pdb)


test_mol = Mol(name="TestMol", records=final_atoms)
print("Building molecular graph...")
graph = MolGraph(test_mol)

moving_atoms = [atom for atom in final_atoms if atom.res_seq == 1]
static_atoms = [atom for atom in final_atoms if atom.res_seq != 1]

print(f"Testing {len(moving_atoms)} moving atoms against {len(static_atoms)} static atoms...")
checker = StericChecker(graph, moving_atoms, static_atoms)
total_clashes = checker.check_clashes()

print(f"Total clashes found: {total_clashes}")