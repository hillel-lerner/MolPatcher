import argparse
import sys, os, shutil

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyBuilder
from mol_patcher.geometry import PatchAligner
from mol_patcher import utilities

def run_patch(pdb_file, res_id, chain, itp_file, outdir, ff_path=None):

    run_folder_name = f"patched_{chain}_{res_id}"
    final_outdir = os.path.join(outdir, run_folder_name)
    os.makedirs(final_outdir, exist_ok=True)

    cdir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(cdir, 'pdbs', pdb_file)
    itp_path = os.path.join(cdir, 'topology_files', itp_file)
    
    # Extract the molecule name for the ITP file (e.g., "PROE" from "PROE.itp")
    base_mol_name = itp_file.split('.')[0]
    
    # Load Patch Data 
    loader = Patchloader()
    pfp_atoms = loader.get_pfp_pdb()
    patch_mol = Mol("patch", pfp_atoms, [], [], [], [], [])
    loader.get_pfp_itp(patch_mol)

    # Load Base Protein Data
    headers, base_records, ters, t_name = PdbParser.read_file(pdb_path)
    base_mol = Mol(name=t_name, records=base_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
    base_mol.load_itp(itp_path)

    # Locate Anchors
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

    # Align the Patch
    aligner = PatchAligner(patch_mol.records, pfp_anchors, target_anchors)
    aligned_pfp_atoms = aligner.implement_align()

    # Execute the Stitcher
    # This automatically handles topology mapping, atom typing, and electrostatic balancing
    stitcher = Stitcher(base_mol, patch_mol, target_res_id=res_id)
    stitched_mol = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_pfp_atoms, 
        target_reference=target_anchors[0], 
        target_anchors=target_anchors
    )
    
    # Save outputs
    pdb_outfile = os.path.join(final_outdir, f"patched_{chain}_{res_id}.pdb")
    itp_outfile = os.path.join(final_outdir, f"patched_{chain}_{res_id}.itp")

    PdbBuilder(pdb_outfile, stitched_mol.records, headers, ters).write_pdb()
    TopologyBuilder(stitched_mol, itp_outfile, base_mol_name).write_itp()

    # Stage forcefield (optional)
    if ff_path is not None:
        if os.path.exists(ff_path):
            dest_ff_dir = os.path.join(final_outdir, "toppar")
            shutil.copytree(ff_path, dest_ff_dir, dirs_exist_ok=True)
            print(f"Master forcefield copied to {dest_ff_dir}")
        else:
            print(f"Warning: Forcefield path {ff_path} not found. Skipping copy.")


    
    box_size, applied_buffer = utilities.get_optimal_box_size(
        stitched_mol.records, 
        buffer_percent=0.33,
        min_buffer_nm=3.0 
    )
    
    print(f"\nSuccessfully generated topology: {itp_outfile}")
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
        parser.add_argument("-o", "--outdir", type=str, default=os.getcwd(), help="Directory to save the patched output. Defaults to current working directory.")
        parser.add_argument("--ff", "--forcefield", type=str, default=None, help="Optional: Path to the master forcefield directory to copy into the output folder.")
        
        args = parser.parse_args()
        # Passed the new outdir and ff arguments to the function
        run_patch(args.pdb, args.res, args.chain, args.itp, args.outdir, args.ff)
    else:
        run_patch(
            # For running outside of CLI
            pdb_file="step3_input.pdb", 
            res_id=136, 
            chain=" ",
            itp_file="PROE.itp",
            outdir=os.getcwd() # Added default output directory for non-CLI runs
        )