import argparse
import sys, os, shutil
import json

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyParser, TopologyBuilder
from mol_patcher.geometry import PatchAligner, MolGraph, StericChecker, RotamerSweeper
from mol_patcher import utilities

def run_patch(pdb_file, res_id, chain, itp_file, outdir, itp_chains, config_file, junction_config, copy_ff=False):

    """
    Orchestrates the entire MolPatcher pipeline: loading data, aligning the patch,
    stitching topologies, performing rotamer optimization, and saving final outputs.
    
    :param pdb_file: (str) Absolute path to the input PDB file.
    :param res_id: (int) The target residue ID for the patch.
    :param chain: (str) The target chain ID.
    :param itp_file: (str) Absolute path to the input ITP topology file.
    :param outdir: (str) Directory where the results folder will be created.
    :param copy_ff: (bool) If True, copies the master forcefield to the output directory.
    :return: None. Writes the patched PDB and ITP files to disk.
    """

    target_res = junction_config["target_res_name"]
    new_res = junction_config["new_res_name"]

    run_folder_name = f"patched_{target_res.lower()}_{res_id}"
    final_outdir = os.path.join(outdir, run_folder_name)
    os.makedirs(final_outdir, exist_ok=True)

    pdb_path = pdb_file
    itp_path = itp_file
    
    # Load Patch Data 
    loader = Patchloader(config_file)
    pfp_atoms = loader.get_pfp_pdb()
    patch_mol = Mol(name="patch", records=pfp_atoms)
    loader.get_pfp_itp(patch_mol)

    # Load Base Protein Data
    parser = PdbParser()
    headers, raw_records, t_name = parser.read_file(pdb_path) 
    fixed_chains = parser.fix_chain_id(raw_records)
    base_records = parser.fix_res_num(fixed_chains)

    base_mol = Mol(name=t_name, records=base_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
    
    safe_itp_path = TopologyParser.clean_itp_file(itp_path, final_outdir)
    base_mol.load_itp(safe_itp_path)
    
    # Remove temporary file
    if safe_itp_path != itp_path and os.path.exists(safe_itp_path):
        os.remove(safe_itp_path)

    # Locate Anchors
    chain_check = chain.strip() if chain.strip() else ""
    target_residue_atoms = [r for r in base_mol.records 
                            if r.res_seq == res_id 
                            and (not chain_check or r.chain.strip() == chain_check)]

    if not target_residue_atoms:
        print(f"Error: Could not find residue {res_id} in chain '{chain_check}'.")
        return

    actual_res_name = target_residue_atoms[0].res_name.strip()
    
    if actual_res_name == new_res:
        print(f"Warning: Target residue '{actual_res_name} {res_id}' is already modified by MolPatcher.")
    elif actual_res_name != target_res and actual_res_name != new_res:
        print(f"Error: Target residue {res_id} is '{actual_res_name}', not '{target_res}'. MolPatcher requires a {target_res} residue.")
        return

    target_anchors = []
    
    try:
        for name in junction_config["base_anchors"]:
            anchor = next(r for r in base_mol.records 
                        if r.res_seq == res_id 
                        and r.name.strip() == name 
                        and (not chain_check or r.chain.strip() == chain_check))
            target_anchors.append(anchor)
    except StopIteration:
        print(f"Error: {target_res} {res_id} is missing required anchor atoms ({', '.join(junction_config['base_anchors'])}).")
        return

    patch_anchors = [
        next(r for r in patch_mol.records if r.name.strip() == name)
        for name in junction_config["patch_anchors"]
    ]

    # Align the Patch
    aligner = PatchAligner(patch_mol.records, patch_anchors, target_anchors)
    aligned_pfp_atoms = aligner.implement_align()

    # Execute the Stitcher
    # CHANGED: Passed the config dictionary into Stitcher
    stitcher = Stitcher(base_mol, patch_mol, target_res_id=res_id, config=junction_config, itp_chains=itp_chains)
    stitched_mol = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_pfp_atoms, 
        target_reference=target_anchors[0], 
        target_anchors=target_anchors
    )
    
    # --- OPTIMIZATION (ROTAMER SWEEP) ---
    graph = MolGraph(stitched_mol, distXX=1.90)

    backbone_names = junction_config.get("rigid_backbone", ['N', 'CA', 'C', 'O'])
    moving_atoms = [a for a in stitched_mol.records 
                    if a.res_seq == res_id and a.name.strip() not in backbone_names]
    static_atoms = [a for a in stitched_mol.records if a.res_seq != res_id]

    checker = StericChecker(graph, moving_atoms, static_atoms, tolerance=0.75)
    sweeper = RotamerSweeper(stitched_mol, graph, res_id, junction_config, chain=chain_check)

    print(f"Starting optimization for {len(moving_atoms)} atoms...")
    success = sweeper.run_sweep(checker)

    if not success:
        print(f"Error: Could not find a clash-free conformation for residue {res_id}. Exiting.")
        return

    # Save outputs
    pdb_outfile = os.path.join(final_outdir, f"patched_{target_res.lower()}_{res_id}.pdb")
    itp_outfile = os.path.join(final_outdir, f"patched_{target_res.lower()}_{res_id}.itp")

    PdbBuilder(pdb_outfile, stitched_mol.records, headers).write_pdb()
    TopologyBuilder(stitched_mol, itp_outfile, stitched_mol.name).write_itp()
    
    box_size, applied_buffer = utilities.get_optimal_box_size(
        stitched_mol.records, 
        buffer_percent=0.33,
        min_buffer_nm=3.0 
    )
    
    print(f"\nSuccessfully generated topology: {itp_outfile}")
    print(f"   --> Molecule Max Diagonal : {round(box_size - applied_buffer, 3)} nm")
    print(f"   --> Applied PBC Buffer    : {applied_buffer} nm")
    print(f"   --> Optimal Cubic Box     : {box_size} nm")


    if copy_ff:
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        ff_src = os.path.join(project_root, 'forcefields', 'forcefield_master.itp')
        ff_dst = os.path.join(final_outdir, 'forcefield_master.itp')

        if os.path.exists(ff_src):
            shutil.copy2(ff_src, ff_dst)
            print(f"   --> Copied master forcefield to: {final_outdir}")
        else:
            print(f"Warning: Master forcefield not found at {ff_src}")

def main():

    """
    Command-line entry point. Parses CLI arguments and triggers the run_patch pipeline.
    """

    parser = argparse.ArgumentParser(description="Patch a ligand into a protein residue.")
    parser.add_argument("-pdb", "--pdb", required=True, help="Input PDB file")
    parser.add_argument("-itp", "--itp", required=True, help="Input ITP file")
    parser.add_argument("-res", "--res", type=int, required=True, help="Target residue ID")
    parser.add_argument("-c", "-chain", "--chain", default=" ", help="Chain ID")
    parser.add_argument("-itp_chains", nargs="+", default=[], help="Chains included in the ITP (e.g., -itp_chains B C)")
    
    parser.add_argument("-config", "--config", required=True, help="Path to JSON configuration file")
    
    parser.add_argument("-o", "--outdir", default=os.getcwd(), help="Output directory")
    parser.add_argument("-ff", "--ff", action="store_true", help="Copy master forcefield")
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        junction_config = json.load(f)
        
    config_basename = os.path.basename(args.config)
    
    pdb_abs = os.path.abspath(args.pdb)
    itp_abs = os.path.abspath(args.itp)

    itp_chains_list = args.itp_chains if args.itp_chains else [args.chain.strip()]
    
    run_patch(pdb_abs, args.res, args.chain, itp_abs, args.outdir, itp_chains_list, config_basename, junction_config, copy_ff=args.ff)

if __name__ == "__main__":
    main()