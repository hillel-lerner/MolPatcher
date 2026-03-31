import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.geometry import PatchAligner
from mol_patcher.pdb_io import PdbParser
from mol_patcher.mol_record import Mol
from mol_patcher.topology_io import TopologyBuilder

if __name__ == "__main__":
    print("Running Refactored Stitcher Test...")

    cdir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(cdir)
    pdb_dir = os.path.join(project_root, 'pdbs')
    top_dir = os.path.join(project_root, 'topology_files')
    
    pdb_file = os.path.join(pdb_dir, "step3_input.pdb")
    itp_file = os.path.join(top_dir, "PROE.itp")

    loader = Patchloader()
    pfp_atoms = loader.get_pfp_pdb()
    patch_mol_obj = Mol("patch", pfp_atoms, [], [], [], [], [])
    loader.get_pfp_itp(patch_mol_obj) # Loads the real ITP data into the empty lists

    loader.get_pfp_itp(patch_mol_obj) # Loads the real ITP data into the empty lists

    _, base_records, _, _ = PdbParser.read_file(pdb_file)
    base_mol_obj = Mol("protein", base_records, [], [], [], [], [])
    base_mol_obj.load_itp(itp_file)

    # Define the target residue dynamically for the test
    test_res = 136

    stitcher = Stitcher(base_mol_obj, patch_mol_obj, target_res_id=test_res)

    # Use the dynamic variable in the exact-order extraction
    target_anchors = [
        next(r for r in base_mol_obj.records if r.res_seq == test_res and r.name.strip() == "CE"),
        next(r for r in base_mol_obj.records if r.res_seq == test_res and r.name.strip() == "CD"),
        next(r for r in base_mol_obj.records if r.res_seq == test_res and r.name.strip() == "NZ")
    ]

    pfp_anchors = [
        next(r for r in patch_mol_obj.records if r.name.strip() == "C10"),
        next(r for r in patch_mol_obj.records if r.name.strip() == "C11"),
        next(r for r in patch_mol_obj.records if r.name.strip() == "N")
    ]

    aligner = PatchAligner(patch_mol_obj.records, pfp_anchors, target_anchors)
    aligned_patch_atoms = aligner.implement_align()

    target_ref = target_anchors[0]
    stitched_mol = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_patch_atoms, 
        target_reference=target_ref, 
        target_anchors=target_anchors
    )
    print(f"SUCCESS: Stitched molecule has {len(stitched_mol.records)} atoms.")
    
    out_pdb_path = os.path.join(pdb_dir, "test_stitched.pdb")
    
    with open(out_pdb_path, 'w') as f:
        for atom in stitched_mol.records:
            occ = getattr(atom, 'occupancy', 1.00)
            temp = getattr(atom, 'temp_factor', 0.00)
            
            line = "{:<6s}{:>5d} {:^4s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}\n".format(
                str(atom.record_type),
                int(atom.serial),   
                str(atom.name),
                " ", 
                str(atom.res_name),
                str(atom.chain),
                int(atom.res_seq),  
                " ", 
                float(atom.x), float(atom.y), float(atom.z),
                float(occ),
                float(temp),
                str(atom.seg_id)
            )
            f.write(line)
            
    print(f"Saved test output to: {out_pdb_path}")

    out_itp_path = os.path.join(top_dir, "test_stitched.itp")
    
    builder = TopologyBuilder(stitched_mol, out_itp_path, "PROE")
    builder.write_itp()
    
    print(f"Saved topology test output to: {out_itp_path}")