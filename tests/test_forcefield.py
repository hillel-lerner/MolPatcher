import os
import sys

if os.path.dirname(os.path.dirname(os.path.abspath(__file__))) not in sys.path:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mol_patcher import stitcher
from mol_patcher.combine_ff import ForceField
from mol_patcher.stitcher import Stitcher

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run_test():
    """
    Tests the ForceField class by mocking a request from the Stitcher
    and verifying the output file generation.
    """

    #    pdb_dir = os.path.join(project_root, "pdbs")
    #   itp_dir = os.path.join(project_root, "topology_files")

    #    base_pdb = os.path.join(pdb_dir, "nipah_fixed.pdb")
    #    base_itp = os.path.join(itp_dir, "nipah.itp")

    #   patch_pdb = os.path.join(pdb_dir, "m9.pdb")
    #    patch_itp = os.path.join(itp_dir, "CARB-m9.itp")

    print("--- Starting ForceField Extraction Test ---\n")

    charmm_files = [
        "/home/hillel/Data_Files/hillel/charmm36_files/toppar/par_all36_carb.prm",
        "/home/hillel/Data_Files/hillel/charmm36_files/toppar/par_all36m_prot.prm",
        "/home/hillel/Data_Files/hillel/charmm36_files/toppar/stream/carb/toppar_all36_carb_glycopeptide.str",
    ]

    output_itp = os.path.join(project_root, "forcefields", "test_junction.itp")

    mock_requests = {
        "atomtypes": [
            ("CC2O1",),
            ("NC2D1",),
            ("HCP1",),
            ("CC3162",),
            ("HCA1",),
            ("OC3C61",),
            ("CC3161",),
        ],
        "bonds": [("NC2D1", "CC3162")],
        "angles": [("CC2O1", "NC2D1", "CC3162")],
        "dihedrals_proper": [("CC2O1", "NC2D1", "CC3162", "HCA1")],
        "dihedrals_improper": [],
    }

    print("1. Initializing and reading database...")
    ff_scraper = ForceField(charmm_files)
    ff_scraper.read_database()

    print("2. Extracting requested parameters...")
    ff_scraper.extract_junction_params(mock_requests)

    print("3. Writing output file...")
    ff_scraper.write_ff(output_itp)
    if os.path.exists(output_itp):
        print(f"\nSUCCESS: {output_itp} was generated!")
        print("--- File Contents ---")

        with open(output_itp, "r") as f:
            print(f.read())
    else:
        print(f"\nFAILURE: {output_itp} was not found.")


if __name__ == "__main__":
    run_test()
