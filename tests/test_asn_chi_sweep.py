import os
import sys
import json

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph, RotamerSweeper

if __name__ == "__main__":
    print("Testing ASN Rotamer Physics...\n")

    pdb_dir = os.path.join(project_root, 'pdbs')
    
    base = os.path.join(pdb_dir, "nipah_fixed.pdb")
    glycan = os.path.join(pdb_dir, "m9.pdb") 
    target_res_seq = 67 
    chain_id = "A"

    config_path = os.path.join(project_root, 'configs', 'nglycan_asn.json')

    # Load Config
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Load PDB
    parser = PdbParser()
    headers, records, mol_name = parser.read_file(base)
    mol = Mol(name=mol_name, records=records)
    
    print("Building MolGraph (isolates the sidechain from the backbone)...")

    graph = MolGraph(mol) 


    sweeper = RotamerSweeper(mol, graph, target_res_seq, config, chain=chain_id)

    canonical_rotamers = config.get("canonical_rotamers", [])
    if not canonical_rotamers:
        print("Error: No canonical rotamers found in config.")
        sys.exit()

    print(f"\nFound {len(canonical_rotamers)} canonical poses. Generating files...")

    # Force the sweeper to apply each pose and save the result
    for i, pose in enumerate(canonical_rotamers):
        print(f"Applying Pose {i+1}: Chi1={pose[0]}, Chi2={pose[1]}")
        sweeper.apply_pose(pose)

        out_pdb = os.path.join(pdb_dir, f"asn_pose_{i+1}.pdb")
        builder = PdbBuilder(out_pdb, mol.records, headers)
        builder.write_pdb()

    print(f"\nTest Complete. Wrote {len(canonical_rotamers)} PDBs to {pdb_dir}")