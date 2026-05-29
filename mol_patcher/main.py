import argparse
import os
import shutil
import json

from networkx import config

from mol_patcher.stitcher import Stitcher, Patchloader
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyParser, TopologyBuilder
from mol_patcher.geometry import PatchAligner, MolGraph, StericChecker
from mol_patcher.sweeper import SweepConductor
from mol_patcher import utilities
from mol_patcher.combine_ff import ForceField


def run_patch(
    patch_pdb,
    patch_itp,
    base_pdb,
    base_itp,
    patch_resid,
    base_resid,
    patch_chain,
    base_chain,
    outdir,
    itp_chains,
    config_file,
    junction_config,
    copy_ff=None,
):
    """
    Orchestrates the entire MolPatcher pipeline: loading data, aligning the patch,
    stitching topologies, performing rotamer optimization, and saving final outputs.

    :param pdb_file: (str) Absolute path to the input PDB file.
    :param res_id: (int) The target residue ID for the patch.
    :param chain: (str) The target chain ID.
    :param itp_file: (str) Absolute path to the input ITP topology file.
    :param outdir: (str) Directory where the results folder will be created.
    :param copy_ff: (list) Optional list of forcefield files to copy to the output directory.
    :return: None. Writes the patched PDB and ITP files to disk.
    """

    if copy_ff is None:
        copy_ff

    target_res = junction_config["target_res_name"]
    new_res = junction_config["new_res_name"]

    run_folder_name = f"patched_{target_res.lower()}_{base_resid}_PRO{base_chain}"
    final_outdir = os.path.join(outdir, run_folder_name)
    infiles_dir = os.path.join(final_outdir, "infiles")
    toppar_dir = os.path.join(final_outdir, "toppar")
    os.makedirs(infiles_dir, exist_ok=True)
    os.makedirs(toppar_dir, exist_ok=True)

    if copy_ff:
        for ff_file in copy_ff:
            if os.path.exists(ff_file):
                shutil.copy2(ff_file, toppar_dir)
            else:
                print(f"Warning: Additional forcefield file not found: {ff_file}")

    shutil.copy2(base_pdb, infiles_dir)
    shutil.copy2(base_itp, infiles_dir)
    shutil.copy2(patch_pdb, infiles_dir)
    shutil.copy2(patch_itp, infiles_dir)

    # Load Patch Data
    parser = PdbParser()
    _, patch_atoms, patch_name = parser.read_file(patch_pdb)
    patch_mol = Mol(name=patch_name, records=patch_atoms)
    patch_mol.load_itp(patch_itp)

    parser = PdbParser()
    headers, raw_records, t_name = parser.read_file(base_pdb)
    base_records = parser.fix_chain_id(raw_records)

    base_mol = Mol(
        name=t_name,
        records=base_records,
        atoms=[],
        bonds=[],
        pairs=[],
        angles=[],
        dihs=[],
    )

    safe_itp_path = TopologyParser.clean_itp_file(base_itp, final_outdir)
    base_mol.load_itp(safe_itp_path)

    # Remove temporary file
    if safe_itp_path != base_itp and os.path.exists(safe_itp_path):
        os.remove(safe_itp_path)

    # Locate Anchors
    chain_check = base_chain.strip() if base_chain.strip() else ""
    target_residue_atoms = [
        r
        for r in base_mol.records
        if r.res_seq == base_resid
        and (not chain_check or r.chain.strip() == chain_check)
    ]

    if not target_residue_atoms:
        print(f"Error: Could not find residue {base_resid} in chain '{chain_check}'.")
        return

    actual_res_name = target_residue_atoms[0].res_name.strip()

    if actual_res_name == new_res:
        print(
            f"Warning: Target residue '{actual_res_name} {base_resid}' is already modified by MolPatcher."
        )
    elif actual_res_name != target_res and actual_res_name != new_res:
        print(
            f"Error: Target residue {base_resid} is '{actual_res_name}', not '{target_res}'. MolPatcher requires a {target_res} residue."
        )
        return

    target_anchors = []

    try:
        for name in junction_config["base_anchors"]:
            anchor = next(
                r
                for r in base_mol.records
                if r.res_seq == base_resid
                and r.name.strip() == name
                and (not chain_check or r.chain.strip() == chain_check)
            )
            target_anchors.append(anchor)
    except StopIteration:
        print(
            f"Error: {target_res} {base_resid} is missing required anchor atoms ({', '.join(junction_config['base_anchors'])})."
        )
        return

    valid_patch_atoms = []

    is_polymer = "patch_merged_atom" not in junction_config

    for r in patch_mol.records:
        if is_polymer:
            match_res = True
        else:
            match_res = True if patch_resid is None else (r.res_seq == patch_resid)

        match_chain = (
            True if patch_chain == " " else (r.chain.strip() == patch_chain.strip())
        )

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
        aligned_patch_atoms = aligner.align_from_template(
            base_records, valid_patch_atoms, junction_config, base_resid
        )
    else:
        # Small Molecule/PFP attachment: 3-point spatial translation + RotamerSweeper
        aligned_patch_atoms = aligner.implement_align()

    # Execute the Stitcher
    stitcher = Stitcher(
        base_mol,
        patch_mol,
        target_res_id=base_resid,
        config=junction_config,
        itp_chains=itp_chains,
    )
    stitched_mol, graph, junction_log = stitcher.stitch_molecules(
        aligned_patch_atoms=aligned_patch_atoms,
        target_reference=target_anchors[0],
        target_anchors=target_anchors,
    )

    # --- OPTIMIZATION (ROTAMER SWEEP) ---
    if "chi_definitions" in junction_config:
        print(
            f"Flexible template detected. Starting optimization for {target_res} {base_resid}..."
        )

        moving_atoms, static_atoms = utilities.identify_optimization_clusters(
            stitched_mol, base_resid, chain_check, junction_config, base_records
        )

        # Initialize the unified scoring conductor
        conductor = SweepConductor(
            stitched_mol, graph, base_resid, junction_config, chain_check
        )

        if not conductor.optimize(moving_atoms, static_atoms):
            print(f"Warning: Could not calculate any poses for residue {base_resid}.")
    else:
        print("Rigid template alignment detected. Bypassing rotamer optimization.")

    # Save outputs
    pdb_outfile = os.path.join(
        final_outdir, f"patched_{target_res.lower()}_{base_resid}.pdb"
    )
    itp_outfile = os.path.join(
        toppar_dir, f"patched_{target_res.lower()}_{base_resid}.itp"
    )

    PdbBuilder(pdb_outfile, stitched_mol.records, headers).write_pdb()
    TopologyBuilder(stitched_mol, itp_outfile, stitched_mol.name).write_itp()

    box_size, applied_buffer = utilities.get_optimal_box_size(
        stitched_mol.records, buffer_percent=0.33, min_buffer_nm=3.0
    )

    ff_refs = junction_config["charmm_forcefield_files"]
    print("Scraping junction parameters from CHARMM database...")
    ff_scraper = ForceField(ff_refs)
    ff_scraper.read_database()

    ff_requests = stitcher.generate_ff_requests(stitched_mol, junction_log)
    ff_scraper.extract_junction_params(ff_requests)
    junction_ff_path = os.path.join(
        toppar_dir, f"junction_ff_{target_res.lower()}_{base_resid}.itp"
    )
    ff_scraper.write_ff(junction_ff_path)

    print(f"   --> Wrote junction parameters to: {junction_ff_path}")

    real_mol_name = stitched_mol.name
    if stitched_mol.moltype_section:
        lines = stitched_mol.moltype_section.strip().split("\n")
        for line in lines:
            if line and not line.startswith("[") and not line.startswith(";"):
                real_mol_name = line.split()[0]  # Grabs the first string (e.g., "PROB")
                break

    # ---------------------------------------------------------
    # TERMINAL OUTPUT & LOG FILE GENERATION
    # ---------------------------------------------------------
    print(f"\nSuccessfully generated topology: {itp_outfile}")
    print(f"   --> Molecule Max Diagonal : {round(box_size - applied_buffer, 3)} nm")
    print(f"   --> Applied PBC Buffer    : {applied_buffer} nm")
    print(f"   --> Optimal Cubic Box     : {box_size} nm")

    log_file_path = os.path.join(
        final_outdir, f"patched_{target_res.lower()}_{base_resid}.log"
    )
    with open(log_file_path, "w") as log:
        log.write("========== MOLPATCHER EXECUTION LOG ==========\n")
        log.write(f"Target : {target_res} {base_resid} (Chain {chain_check})\n")
        log.write(f"Patch  : {os.path.basename(patch_pdb)}\n")
        log.write("==============================================\n\n")

        log.write("--- FILE I/O ---\n")
        log.write(
            f"INFILES:\n  Base:  {base_pdb} & {base_itp}\n  Patch: {patch_pdb} & {patch_itp}\n"
        )
        log.write(f"OUTFILES:\n  PDB:   {pdb_outfile}\n  ITP:   {itp_outfile}\n\n")

        log.write("--- JUNCTION TOPOLOGY ---\n")
        log.write(
            f"Bond Created: {junction_log['bond'][0]} - {junction_log['bond'][1]}\n"
        )

        log.write(f"Angles Created ({len(junction_log['angles'])}):\n")
        for a1, a2, a3 in junction_log["angles"]:
            log.write(f"  {a1} - {a2} - {a3}\n")

        log.write(f"Dihedrals Created ({len(junction_log['dihs'])}):\n")
        for d1, d2, d3, d4 in junction_log["dihs"]:
            log.write(f"  {d1} - {d2} - {d3} - {d4}\n")
        log.write("\n")

        log.write("--- GROMACS TOPOLOGY INSTRUCTIONS ---\n")
        log.write(
            "1. Add the following to your master .top file (after the forcefield #includes):\n"
        )
        log.write(f'   #include "{os.path.basename(itp_outfile)}"\n\n')
        log.write(
            "2. Add the following to the [ molecules ] directive at the bottom of your .top file:\n"
        )
        log.write(f"   {stitched_mol.name:<15} 1\n\n")

        log.write("--- PERIODIC BOUNDARY CONDITIONS (PBC) ---\n")
        log.write(f"Molecule Max Diagonal : {round(box_size - applied_buffer, 3)} nm\n")
        log.write(f"Applied PBC Buffer    : {applied_buffer} nm\n")
        log.write(f"Optimal Cubic Box     : {box_size} nm\n\n")

        log.write("--- APPLIED CONFIGURATION ---\n")
        log.write(
            f"Residue Mutation : {junction_config.get('target_res_name')} -> {junction_config.get('new_res_name')}\n"
        )
        log.write(
            f"Anchors Used     : Base ({junction_config.get('base_bridge')}), Patch ({junction_config.get('patch_bridge')})\n"
        )

        if "geometry" in junction_config:
            geom = junction_config["geometry"]
            log.write(
                f"Preset Geometry  : Bond={geom.get('bond_length')}Å, Angle={geom.get('angle_target')}°, Omega={geom.get('omega_target')}°, Phi={geom.get('phi_target')}°\n"
            )
        else:
            log.write("Preset Geometry  : None (Flexible alignment used)\n")
        log.write("\n")

        log.write("--- JUNCTION FORCEFIELD PARAMETERS ADDED ---\n")
        if hasattr(ff_scraper, "junction_output"):
            params_found = False
            for section, lines in ff_scraper.junction_output.items():
                if lines:
                    params_found = True
                    log.write(f"[ {section} ]\n")
                    for line in sorted(list(lines)):
                        log.write(line)
                    log.write("\n")
            if not params_found:
                log.write("No junction parameters were extracted or added.\n\n")

    print(f"   --> Wrote execution log to: {log_file_path}")


def main():
    """
    Command-line entry. Parses CLI arguments and triggers the run_patch pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Patch a ligand into a protein residue."
    )
    parser.add_argument(
        "-config", "--config", required=True, help="Path to JSON configuration file"
    )
    parser.add_argument(
        "-base",
        nargs=4,
        metavar=("PDB", "ITP", "RES", "CHAIN"),
        required=True,
        help="Base protein: [PDB] [ITP] [Residue ID] [Chain]",
    )
    parser.add_argument(
        "-patch",
        nargs="+",
        metavar="ARG",
        required=True,
        help="Patch molecule: [PDB] [ITP] [Optional Res ID] [Optional Chain]",
    )
    parser.add_argument("-o", "--outdir", default=os.getcwd(), help="Output directory")
    parser.add_argument(
        "-ff",
        "--forcefield",
        dest="copy_ff",
        nargs="*",
        default=[],
        help="Optional list of additional forcefield files to copy to toppar directory",
    )

    args = parser.parse_args()
    config_path = args.config

    if not os.path.exists(config_path):
        main_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(main_dir)
        internal_config_path = os.path.join(project_root, "configs", args.config)

        if not internal_config_path.endswith(".json"):
            internal_config_path += ".json"

        if os.path.exists(internal_config_path):
            config_path = internal_config_path
        else:
            raise FileNotFoundError(
                f"Config file not found locally or in MolPatcher/configs: {args.config}"
            )

    with open(config_path, "r") as f:
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
        itp_chains=[base_chain.strip()] if base_chain.strip() else [],
        config_file=config_basename,
        junction_config=junction_config,
        copy_ff=args.copy_ff,
    )


if __name__ == "__main__":
    main()
