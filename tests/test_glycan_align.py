import os
import sys
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser
from mol_patcher.mol_record import Mol
from mol_patcher.geometry import MolGraph
from mol_patcher.utilities import get_distance


if __name__ == "__main__":
    print("Testing Glycan Bond Alignment...\n")

    pdb_dir = os.path.join(project_root, 'pdbs')
    patch_pdb = os.path.join(pdb_dir, "glycan_only.pdb")
    base_pdb = os.path.join(pdb_dir, 'no_glycan.pdb')
    out_pdb = "aligned_glycan_test.pdb"

    pdb_parser = PdbParser()
    _, records, _ = pdb_parser.read_file(base_pdb)
    _, records, _ = pdb_parser.read_file(patch_pdb)
