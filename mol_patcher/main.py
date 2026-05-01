import argparse
import sys, os, shutil
import json

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyParser, TopologyBuilder
from mol_patcher.geometry import PatchAligner, MolGraph, StericChecker, RotamerSweeper
from mol_patcher import utilities

def run_patch(patch_pdb, patch_itp, base_pdb, base_itp, patch_resid, base_resid, patch_chain, base_chain, outdir, itp_chains, config_file, junction_config, copy_ff=False):

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

    run_folder_name = f"patched_{target_res.lower()}_{base_resid}_PRO{base_chain.strip()}"
    final_outdir = os.path.join(outdir, run_folder_name)
    infiles_dir = os.path.join(final_outdir, "infiles")
    os.makedirs(infiles_dir, exist_ok=True)

    shutil.copy2(base_pdb, infiles_dir)
    shutil.copy2(base_itp, infiles_dir)
    shutil.copy2(patch_pdb, infiles_dir)
    shutil.copy2(patch_itp, infiles_dir)

    # Load Patch Data
    parser = PdbParser()
    _, patch_atoms, patch_name = parser.read_file(patch_pdb)
    patch_mol = Mol(name=patch_name, records=patch_atoms)
    patch_mol.load_itp(patch_itp)
    
    # Load Base Protein Data
    parser = PdbParser()
    headers, raw_records, t_name = parser.read_file(base_pdb) 
    fixed_chains = parser.fix_chain_id(raw_records)
    base_records = parser.fix_res_num(fixed_chains)

    base_mol = Mol(name=t_name, records=base_records, atoms=[], bonds=[], pairs=[], angles=[], dihs=[])
    
    safe_itp_path = TopologyParser.clean_itp_file(base_itp, final_outdir)
    base_mol.load_itp(safe_itp_path)
    
    # Remove temporary file
    if safe_itp_path != base_itp and os.path.exists(safe_itp_path):
        os.remove(safe_itp_path)

    # Locate Anchors
    chain_check = base_chain.strip() if base_chain.strip() else ""
    target_residue_atoms = [r for r in base_mol.records 
                            if r.res_seq == base_resid 
                            and (not chain_check or r.chain.strip() == chain_check)]

    if not target_residue_atoms:
        print(f"Error: Could not find residue {base_resid} in chain '{chain_check}'.")
        return

    actual_res_name = target_residue_atoms[0].res_name.strip()
    
    if actual_res_name == new_res:
        print(f"Warning: Target residue '{actual_res_name} {base_resid}' is already modified by MolPatcher.")
    elif actual_res_name != target_res and actual_res_name != new_res:
        print(f"Error: Target residue {base_resid} is '{actual_res_name}', not '{target_res}'. MolPatcher requires a {target_res} residue.")
        return

    target_anchors = []
    
    try:
        for name in junction_config["base_anchors"]:
            anchor = next(r for r in base_mol.records 
                        if r.res_seq == base_resid 
                        and r.name.strip() == name 
                        and (not chain_check or r.chain.strip() == chain_check))
            target_anchors.append(anchor)
    except StopIteration:
        print(f"Error: {target_res} {base_resid} is missing required anchor atoms ({', '.join(junction_config['base_anchors'])}).")
        return
    

    valid_patch_atoms = []
    for r in patch_mol.records:
        match_res = True if patch_resid is None else (r.res_seq == patch_resid)
        match_chain = True if patch_chain == " " else (r.chain.strip() == patch_chain.strip())
        
        if match_res and match_chain:
            valid_patch_atoms.append(r)
            
    if not valid_patch_atoms:
        print("Error: No patch atoms found matching the specified residue/chain.")
        return

    patch_anchors = [
        next(r for r in valid_patch_atoms if r.name.strip() == name) 
        for name in junction_config["patch_anchors"]
    ]

    # Align the Patch
    aligner = PatchAligner(valid_patch_atoms, patch_anchors, target_anchors)
    if "geometry" in junction_config:
        # Glycan attachment: preset torsions based on the template
        aligned_patch_atoms = aligner.align_from_template(base_records, valid_patch_atoms, junction_config, base_resid)
    else:
        # Small Molecule/PFP attachment: 3-point spatial translation + RotamerSweeper
        aligned_patch_atoms = aligner.implement_align()

    # Execute the Stitcher
    stitcher = Stitcher(base_mol, patch_mol, target_res_id=base_resid, config=junction_config, itp_chains=itp_chains)
    stitched_mol = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_patch_atoms, 
        target_reference=target_anchors[0], 
        target_anchors=target_anchors
    )
    
    # --- OPTIMIZATION (ROTAMER SWEEP) ---
    # Only execute if the JSON template provides specific dihedral definitions (e.g., PFP-Lysine)
    if "chi_definitions" in junction_config:
        graph = MolGraph(stitched_mol, distXX=1.90)
    
        backbone_names = junction_config.get("rigid_backbone", ['N', 'CA', 'C', 'O'])
        moving_atoms = [a for a in stitched_mol.records 
                        if a.res_seq == base_resid and a.name.strip() not in backbone_names]
        static_atoms = [a for a in stitched_mol.records if a.res_seq != base_resid]
    
        checker = StericChecker(graph, moving_atoms, static_atoms, tolerance=0.75)
        sweeper = RotamerSweeper(stitched_mol, graph, base_resid, junction_config, chain=chain_check)
    
        print(f"Starting optimization for {len(moving_atoms)} atoms...")
        success = sweeper.run_sweep(checker)
    
        if not success:
            print(f"Error: Could not find a clash-free conformation for residue {base_resid}. Exiting.")
            return
    else:
        # For glycans, bypass the sweep to preserve the specified linkage geometry
        print("Rigid template alignment detected. Bypassing rotamer optimization.")

    # Save outputs
    pdb_outfile = os.path.join(final_outdir, f"patched_{target_res.lower()}_{base_resid}.pdb")
    itp_outfile = os.path.join(final_outdir, f"patched_{target_res.lower()}_{base_resid}.itp")

    PdbBuilder(pdb_outfile, stitched_mol.records, headers).write_pdb()
    TopologyBuilder(stitched_mol, itp_outfile, stitched_mol.name).write_itp()
    
    box_size, applied_buffer = utilities.get_optimal_box_size(
        stitched_mol.records, 
        buffer_percent=0.33,
        min_buffer_nm=3.0 
    )

    real_mol_name = stitched_mol.name
    if stitched_mol.moltype_section:
        lines = stitched_mol.moltype_section.strip().split('\n')
        for line in lines:
            if line and not line.startswith('[') and not line.startswith(';'):
                real_mol_name = line.split()[0] # Grabs the first string (e.g., "PROB")
                break

    # ---------------------------------------------------------
    # TERMINAL OUTPUT & LOG FILE GENERATION
    # ---------------------------------------------------------
    print(f"\nSuccessfully generated topology: {itp_outfile}")
    print(f"   --> Molecule Max Diagonal : {round(box_size - applied_buffer, 3)} nm")
    print(f"   --> Applied PBC Buffer    : {applied_buffer} nm")
    print(f"   --> Optimal Cubic Box     : {box_size} nm")

    log_file_path = os.path.join(final_outdir, "patcher.log")
    with open(log_file_path, "w") as log:
        log.write(f"========== MOLPATCHER EXECUTION LOG ==========\n")
        log.write(f"Target : {target_res} {base_resid} (Chain {chain_check})\n")
        log.write(f"Patch  : {os.path.basename(patch_pdb)}\n")
        log.write(f"==============================================\n\n")
        
        log.write("--- GROMACS TOPOLOGY INSTRUCTIONS ---\n")
        log.write("1. Add the following to your master .top file (after the forcefield #includes):\n")
        log.write(f"   #include \"{os.path.basename(itp_outfile)}\"\n\n")
        log.write("2. Add the following to the [ molecules ] directive at the bottom of your .top file:\n")
        # GROMACS uses the specific moltype name from the [ moleculetype ] section
        log.write(f"   {stitched_mol.name:<15} 1\n\n")
        
        log.write("--- PERIODIC BOUNDARY CONDITIONS (PBC) ---\n")
        log.write(f"Molecule Max Diagonal : {round(box_size - applied_buffer, 3)} nm\n")
        log.write(f"Applied PBC Buffer    : {applied_buffer} nm\n")
        log.write(f"Optimal Cubic Box     : {box_size} nm\n\n")
        
        log.write("--- APPLIED CONFIGURATION ---\n")
        # Dumps the dictionary back into a formatted JSON string
        log.write(json.dumps(junction_config, indent=4))
        log.write("\n")

    print(f"   --> Wrote execution log to: {log_file_path}")


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
    Command-line entry. Parses CLI arguments and triggers the run_patch pipeline.
    """
    parser = argparse.ArgumentParser(description="Patch a ligand into a protein residue.")
    parser.add_argument("-config", "--config", required=True, help="Path to JSON configuration file")
    parser.add_argument("-base", nargs=4, metavar=('PDB', 'ITP', 'RES', 'CHAIN'), required=True, 
                        help="Base protein: [PDB] [ITP] [Residue ID] [Chain]")
    parser.add_argument("-patch", nargs='+', metavar=('PDB', 'ITP', '[RES]', '[CHAIN]'), required=True, 
                        help="Patch molecule: [PDB] [ITP] [Optional Res ID] [Optional Chain]")
    parser.add_argument("-o", "--outdir", default=os.getcwd(), help="Output directory")
    parser.add_argument("-ff", "--ff", action="store_true", help="Copy master forcefield")
    
    args = parser.parse_args()
    config_path = args.config

    if not os.path.exists(config_path):
        main_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(main_dir)
        internal_config_path = os.path.join(project_root, 'configs', args.config)

        if not internal_config_path.endswith('.json'):
            internal_config_path += '.json'

        if os.path.exists(internal_config_path):
            config_path = internal_config_path
        else:
            raise FileNotFoundError(f"Config file not found locally or in MolPatcher/configs: {args.config}")
    
    with open(config_path, 'r') as f:
        junction_config = json.load(f)
        
    config_basename = os.path.basename(config_path)
    
    base_pdb_abs = os.path.abspath(args.base[0])
    base_itp_abs = os.path.abspath(args.base[1])
    base_res = int(args.base[2])
    base_chain = args.base[3]

    patch_pdb_abs = os.path.abspath(args.patch[0])
    patch_itp_abs = os.path.abspath(args.patch[1])
    
    # Standard list length checks for the optional patch parameters
    patch_res = int(args.patch[2]) if len(args.patch) > 2 else None
    patch_chain = args.patch[3] if len(args.patch) > 3 else " "
    
    run_patch(
        patch_pdb=patch_pdb_abs, 
        patch_itp=patch_itp_abs, 
        base_pdb=base_pdb_abs, 
        base_itp=base_itp_abs, 
        patch_resid=patch_res, 
        base_resid=base_res, 
        patch_chain=patch_chain, 
        base_chain=base_chain, 
        outdir=args.outdir, 
        itp_chains=[],
        config_file=config_basename, 
        junction_config=junction_config, 
        copy_ff=args.ff
    )
if __name__ == "__main__":
    main()