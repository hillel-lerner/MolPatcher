import os
import sys
import json

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyParser
from mol_patcher.stitcher import Stitcher
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph, PatchAligner
from mol_patcher.utilities import identify_optimization_clusters
from mol_patcher.sweeper import SweepConductor


if __name__ == "__main__":
    print("Testing Global Orchestrator...\n")

    # --- SET TARGET VARIABLES ---
    pdb_dir = os.path.join(project_root, 'pdbs')
    itp_dir = os.path.join(project_root, "topology_files")
    
    base_pdb = os.path.join(pdb_dir, "nipah.pdb") # The protein
    base_itp = os.path.join(itp_dir, "nipah.itp") 
    
    patch_pdb = os.path.join(pdb_dir, "m9.pdb")         # The glycan
    patch_itp = os.path.join(itp_dir, "CARB-m9.itp")    # The glycan topology
    
    target_res_seq = 67
    chain_id = "A" 
    # -----------------------------------------


    config_path = os.path.join(project_root, 'configs', 'nglycan_asn.json')
    with open(config_path, 'r') as f:
        raw_config = json.load(f)
        junction_config = raw_config[0] if isinstance(raw_config, list) else raw_config

# --- LOAD FILES ---
    parser = PdbParser()
    base_headers, raw_base, base_name = parser.read_file(base_pdb)
    base_records = parser.fix_chain_id(raw_base) 
    _, patch_records, patch_name = parser.read_file(patch_pdb)
    
    # Initialize Mol objects
    base_mol = Mol(name=base_name, records=base_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
    patch_mol = Mol(name=patch_name, records=patch_records)
    base_mol.load_itp(base_itp)
    # Load the full topology directly into the patch molecule (just like main.py!)
    patch_mol.load_itp(patch_itp)

    # --- ALIGNMENT ---
    print(f"Aligning {patch_name} to {base_name}...")
    chain_check = chain_id.strip() if chain_id.strip() else ""
    
    target_anchors = []
    for name in junction_config["base_anchors"]:
        anchor = next(r for r in base_mol.records 
                    if r.res_seq == target_res_seq 
                    and r.name.strip() == name 
                    and (not chain_check or r.chain.strip() == chain_check))
        target_anchors.append(anchor)

    aligner = PatchAligner(patch_records, [], target_anchors)
    aligned_patch_atoms = aligner.align_from_template(base_records, patch_records, junction_config, target_res_seq)

    # --- STITCHING ---
    print(f"Stitching {patch_name} to {base_name} at residue {target_res_seq}...")
    stitcher = Stitcher(base_mol, patch_mol, target_res_id=target_res_seq, config=junction_config, itp_chains=[chain_check] if chain_check else [])
    
    # Notice that your stitcher ALREADY returns the graph! We don't need to rebuild it manually.
    stitched_mol, graph, junction_log = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_patch_atoms, 
        target_reference=target_anchors[0], 
        target_anchors=target_anchors
    )

    # --- INITIALIZE THE ORCHESTRATOR ---
    print("Identifying moving vs static clusters...")
    moving_atoms, static_atoms = identify_optimization_clusters(
        stitched_mol, target_res_seq, chain_id, junction_config, base_records
    )
    
    conductor = SweepConductor(
        stitched_mol, graph, target_res_seq, junction_config, project_root, chain_id
    )

    # --- EXECUTE THE GLOBAL SWEEP ---
    success = conductor.optimize(moving_atoms, static_atoms)

    # --- SAVE THE RESULT ---
    if success:
        out_name = f"nipah_m9_optimized_{target_res_seq}.pdb"
    else:
        out_name = f"nipah_m9_failed_{target_res_seq}.pdb"
        
    out_pdb = os.path.join(pdb_dir, out_name)
    builder = PdbBuilder(out_pdb, stitched_mol.records, base_headers)
    builder.write_pdb()
    
    print(f"\nTest Complete. Output saved to {out_pdb}")