import argparse
import sys, os
from mol_patcher.mol_record import Mol
from mol_patcher.pdb_io import PdbParser, BuildPdb
from mol_patcher.itp_io import ItpParser, BuildItp
from mol_patcher.align_geom import PatchAligner
from mol_patcher import mol_stitcher

def run_patch(pdb_file, res_id, chain, itp_file):

    # Paths (adjust to your structure)
    cdir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(cdir, 'pdbs', pdb_file)
    itp_path = os.path.join(cdir, 'itps', itp_file)
    
    # Load Data
    headers, target_atoms, ters, t_name = PdbParser.read_file(pdb_path)
    target_mol = Mol(name=t_name, records=target_atoms)
    target_mol.load_itp(itp_path)

    pfp_atoms = mol_stitcher.get_pfp_pdb()
    pfp_mol = Mol(name="PFP", records=pfp_atoms)
    mol_stitcher.get_pfp_itp(pfp_mol)

    # Anchors
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

    # Save
    out_pdb = os.path.join(cdir, 'pdbs', f"patched_{res_id}.pdb")
    out_itp = os.path.join(cdir, 'itps', f"patched_{res_id}.itp")
    BuildPdb(out_pdb, stitched_mol.records, headers, ters).write_pdb()
    BuildItp(stitched_mol, out_itp).write_itp()

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