import argparse
import sys, os

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyBuilder
from mol_patcher.geometry import PatchAligner
from mol_patcher import utilities

def run_patch(pdb_file, res_id, chain, itp_file):

    # Paths (adjust to your structure)
    cdir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(cdir, 'pdbs', pdb_file)
    itp_path = os.path.join(cdir, 'topology_files', itp_file)
    
    # Extract the molecule name for the ITP file (e.g., "PROE" from "PROE.itp")
    base_mol_name = itp_file.split('.')[0]
    
    # 1. Load Patch Data via Patchloader
    loader = Patchloader()
    pfp_atoms = loader.get_pfp_pdb()
    patch_mol = Mol("patch", pfp_atoms, [], [], [], [], [])
    loader.get_pfp_itp(patch_mol)

    # 2. Load Base Protein Data
    headers, base_records, ters, t_name = PdbParser.read_file(pdb_path)
    base_mol = Mol(name=t_name, records=base_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
    base_mol.load_itp(itp_path)

    # 3. Locate Anchors
    # Handle empty chains to ensure precise matching
    chain_check = chain.strip() if chain.strip() else ""

    target_anchors = []
    for name in ["CE", "CD", "NZ"]:
        anchor = next(r for r in base_mol.records 
                    if r.res_seq == res_id 
                    and r.name.strip() == name 
                    and (not chain_check or r.chain.strip() == chain_check))
        target_anchors.append(anchor)

    pfp_anchors = [
        next(r for r in patch_mol.records if r.name.strip() == "C10"),
        next(r for r in patch_mol.records if r.name.strip() == "C11"),
        next(r for r in patch_mol.records if r.name.strip() == "N")
    ]

    # 4. Align the Patch
    aligner = PatchAligner(patch_mol.records, pfp_anchors, target_anchors)
    aligned_pfp_atoms = aligner.implement_align()

    # 5. Execute the Stitcher
    # This automatically handles topology mapping, atom typing, and electrostatic balancing
    stitcher = Stitcher(base_mol, patch_mol, target_res_id=res_id)
    stitched_mol = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_pfp_atoms, 
        target_reference=target_anchors[0], 
        target_anchors=target_anchors
    )
    
    # 6. Save outputs
    out_pdb = os.path.join(cdir, 'pdbs', f"patched_{res_id}.pdb")
    out_itp = os.path.join(cdir, 'topology_files', f"patched_{res_id}.itp")
    
    PdbBuilder(out_pdb, stitched_mol.records, headers, ters).write_pdb()
    TopologyBuilder(stitched_mol, out_itp, base_mol_name).write_itp()
    
    print(f"Successfully generated topology: {out_itp}")

    # 7. Box Sizing
    box_size, applied_buffer = utilities.get_optimal_box_size(
        stitched_mol.records, 
        buffer_percent=0.33,
        min_buffer_nm=3.0 
    )
    
    print(f"\nSuccessfully generated topology: {out_itp}")
    print(f"   --> Molecule Max Diagonal : {round(box_size - applied_buffer, 3)} nm")
    print(f"   --> Applied PBC Buffer    : {applied_buffer} nm")
    print(f"   --> Optimal Cubic Box     : {box_size} nm")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(description="Patch a ligand into a protein residue.")
        parser.add_argument("--pdb", required=True, help="Input PDB filename")
        parser.add_argument("--itp", required=True, help="Input ITP filename")
        parser.add_argument("--res", type=int, required=True, help="Target residue ID")
        parser.add_argument("--chain", default=" ", help="Chain ID")
        
        args = parser.parse_args()
        run_patch(args.pdb, args.res, args.chain, args.itp)
    else:
        run_patch(
            # For running outside of CLI
            pdb_file="step3_input.pdb", 
            res_id=136, 
            chain=" ",
            itp_file="PROE.itp"
        )