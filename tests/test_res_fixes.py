import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser, PdbBuilder 

pdb_dir = os.path.join(project_root, 'pdbs')
test_pdb = os.path.join(pdb_dir, "nipah.pdb")
output_pdb = os.path.join(pdb_dir, "nipah_fixed.pdb")


parser = PdbParser()
headers, raw_records, t_name = parser.read_file(test_pdb) 
final_atoms = parser.fix_chain_id(raw_records)



builder = PdbBuilder(output_pdb, final_atoms, headers=headers)
builder.write_pdb()