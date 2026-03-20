import argparse
import sys, os
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, PdbBuilder
from mol_patcher.topology_io import TopologyParser, TopologyBuilder
from mol_patcher.align_geom import PatchAligner
from mol_patcher import mol_stitcher, utilities

def run_patch(pdb_file, res_id, chain, itp_file):

    # Paths (adjust to your structure)
    cdir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(cdir, 'pdbs', pdb_file)
    itp_path = os.path.join(cdir, 'topology_files', itp_file)
    
    # Extract the molecule name for the ITP file (e.g., "PROE" from "PROE.itp")
    base_mol_name = itp_file.split('.')[0]
    
    # Load Data
    headers, target_atoms, ters, t_name = PdbParser.read_file(pdb_path)
    target_mol = Mol(name=t_name, records=target_atoms)
    target_mol.load_itp(itp_path)

    pfp_atoms = mol_stitcher.get_pfp_pdb()
    pfp_mol = Mol(name="PFP", records=pfp_atoms)
    mol_stitcher.get_pfp_itp(pfp_mol)

    # Anchors (Reverted to the last working state for stable geometry)
    pfp_anchors = [
        pfp_mol.get_atom(1, " ", "PFP", "C10"),
        pfp_mol.get_atom(1, " ", "PFP", "C11"),
        pfp_mol.get_atom(1, " ", "PFP", "N")
    ]
    target_anchors = [
        target_mol.get_atom(res_id, chain, "LYS", "CE"),
        target_mol.get_atom(res_id, chain, "LYS", "CD"),
        target_mol.get_atom(res_id, chain, "LYS", "NZ")
    ]

    # Align & Stitch
    aligner = PatchAligner(pfp_mol.records, pfp_anchors, target_anchors)
    aligned_pfp_atoms = aligner.implement_align()

    stitched_mol = mol_stitcher.stitch_molecules(
        target_mol, aligned_pfp_atoms, target_anchors[0], target_anchors, pfp_mol
    )

    # Balance charges
    pfp_junction_charges = {"CD": -0.245, "CE": -0.006, "NZ": -0.092}
    charge_sinks = ["HD1", "HD2", "HE1", "HE2", "CG"] 
    stitched_mol = mol_stitcher.balance_charges(stitched_mol, res_id, pfp_junction_charges, charge_sinks)
    stitched_mol = mol_stitcher.force_integer_charge(stitched_mol, res_id)


    # Update topology atom types
    stitched_mol = mol_stitcher.update_atom_types(
        stitched_mol, 
        res_id, 
        type_updates={
            "NZ": "NG311",
            "HZ1": "HGPAM1", 
            "HZ2": "HGPAM1", 
            "HZ3": "HGPAM1"
        }
    )  
    
    # Save outputs
    out_pdb = os.path.join(cdir, 'pdbs', f"patched_{res_id}.pdb")
    out_itp = os.path.join(cdir, 'topology_files', f"patched_{res_id}.itp")
    
    PdbBuilder(out_pdb, stitched_mol.records, headers, ters).write_pdb()
    
    # Pass the required base_mol_name to TopologyBuilder
    TopologyBuilder(stitched_mol, out_itp, base_mol_name).write_itp()
    print(f"Successfully generated topology: {out_itp}")

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