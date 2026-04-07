import argparse
import sys, os, shutil

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyBuilder
from mol_patcher.geometry import PatchAligner, MolGraph, StericChecker, RotamerSweeper
from mol_patcher import utilities

def run_patch(pdb_file, res_id, chain, itp_file, outdir, ff_path=None):

    run_folder_name = f"patched_lys_{res_id}"
    final_outdir = os.path.join(outdir, run_folder_name)
    os.makedirs(final_outdir, exist_ok=True)

    cdir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(cdir, 'pdbs', pdb_file)
    itp_path = os.path.join(cdir, 'topology_files', itp_file)
    
    # Extract the molecule name for the ITP file (e.g., "PROE" from "PROE.itp")
    base_mol_name = os.path.splitext(os.path.basename(itp_path))[0]
    
    # Load Patch Data 
    loader = Patchloader()
    pfp_atoms = loader.get_pfp_pdb()
    patch_mol = Mol("patch", pfp_atoms, [], [], [], [], [])
    loader.get_pfp_itp(patch_mol)

    # Load Base Protein Data
    parser = PdbParser()
    headers, base_records, t_name = parser.read_file(pdb_path) 
    base_mol = Mol(name=t_name, records=base_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
    base_mol.load_itp(itp_path)

    # Locate Anchors
    # Handle empty chains to ensure precise matching
    chain_check = chain.strip() if chain.strip() else ""
    target_residue_atoms = [r for r in base_mol.records 
                            if r.res_seq == res_id 
                            and (not chain_check or r.chain.strip() == chain_check)]

    # Check if the residue even exists in the PDB file
    if not target_residue_atoms:
        print(f"Error: Could not find residue {res_id} in chain '{chain_check}'.")
        return

    # Check if the residue is actually a Lysine
    actual_res_name = target_residue_atoms[0].res_name.strip()
    if actual_res_name != "LYS" and actual_res_name != "LYX":
        print(f"Error: Target residue {res_id} is '{actual_res_name}', not 'LYS' or 'LYX'. MolPatcher requires a Lysine residue.")
        return
    # -----------------------------

    target_anchors = []
    
    try:
        for name in ["CE", "CD", "NZ"]:
            anchor = next(r for r in base_mol.records 
                        if r.res_seq == res_id 
                        and r.name.strip() == name 
                        and (not chain_check or r.chain.strip() == chain_check))
            target_anchors.append(anchor)
    except StopIteration:
        print(f"Error: Lysine {res_id} is missing required anchor atoms (CE, CD, or NZ).")
        return

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
    
    # --- OPTIMIZATION (ROTAMER SWEEP) ---
    graph = MolGraph(stitched_mol, distXX=1.90)

    backbone_names = ['N', 'H', 'HN', 'CA', 'HA', 'C', 'O']
    moving_atoms = [a for a in stitched_mol.records 
                    if a.res_seq == res_id and a.name.strip() not in backbone_names]
    static_atoms = [a for a in stitched_mol.records if a.res_seq != res_id]

    checker = StericChecker(graph, moving_atoms, static_atoms, tolerance=0.75)
    sweeper = RotamerSweeper(stitched_mol, graph, res_id)

    print(f"Starting optimization for {len(moving_atoms)} atoms...")
    success = sweeper.run_sweep(checker)

    if not success:
        print(f"Error: Could not find a clash-free conformation for residue {res_id}. Exiting.")
        return


    # Save outputs
    pdb_outfile = os.path.join(final_outdir, f"patched_lys_{res_id}.pdb")
    itp_outfile = os.path.join(final_outdir, f"patched_lys_{res_id}.itp")

    PdbBuilder(pdb_outfile, stitched_mol.records, headers).write_pdb()
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

def main():
    parser = argparse.ArgumentParser(description="Patch a ligand into a protein residue.")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--itp", required=True, help="Input ITP file")
    parser.add_argument("--res", type=int, required=True, help="Target residue ID")
    parser.add_argument("--chain", default=" ", help="Chain ID")
    parser.add_argument("-o", "--outdir", default=os.getcwd(), help="Output directory")
    
    args = parser.parse_args()
    
    # Resolve absolute paths based on where the user is currently standing
    pdb_abs = os.path.abspath(args.pdb)
    itp_abs = os.path.abspath(args.itp)
    
    run_patch(pdb_abs, args.res, args.chain, itp_abs, args.outdir)

if __name__ == "__main__":
    main()