import os
import sys
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyParser
from mol_patcher.stitcher import *
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import *
from mol_patcher.utilities import *

if __name__ == "__main__":
    print("Testing Anomer Definitons...\n")

    pdb_dir = os.path.join(project_root, 'pdbs')
    itp_dir = os.path.join(project_root, "topology_files")

    glycan_pdb = os.path.join(pdb_dir, "m9.pdb")
    glycan_itp = os.path.join(itp_dir, "CARB-m9.itp")
    config = "glycan_asparagine.json"
    library_path = os.path.join("configs", "glycosidic_angles.json" )

    parser = PdbParser()
    headers, records, mol_name = parser.read_file(glycan_pdb)

    itp_parser = TopologyParser() 
    _, itp_atoms, itp_bonds, _, _, _ = itp_parser.read_file(glycan_itp)

    mol = Mol(name=mol_name, records=records)
    mol.bonds = itp_bonds 

    itp_to_pdb_map = {}
    pdb_to_itp = Stitcher._build_pdb_to_itp_map(records=records, atoms=itp_atoms)

    for rec_id, itp_atom in pdb_to_itp.items():
        matching_record = next(r for r in records if id(r) == rec_id)
        itp_to_pdb_map[itp_atom.number] = matching_record
    
    graph = MolGraph(mol, itp_to_pdb_map=itp_to_pdb_map)
    
    
    sweeper = GlycanSweeper(mol, graph, config, library_path)
    
    sweeper.find_linkages(patch_chain=" ")

    moving_atoms = mol.records
    static_atoms = []
    
    checker = StericChecker(graph, moving_atoms, static_atoms, tolerance=0.75)
    
    # 2. Execute the sweep
    print("\nStarting Glycosidic Sweep...")
    sweeper.run_sweep(checker)
    
    out_pdb = os.path.join(pdb_dir,  'm9_swept.pdb')
    
    builder = PdbBuilder(out_pdb, mol.records, headers)
    builder.write_pdb()
    
    print(f"\nTest Complete. Wrote optimized glycan to {out_pdb}")
    