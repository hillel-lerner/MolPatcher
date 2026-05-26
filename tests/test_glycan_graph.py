import os
import sys

# Ensure MolPatcher modules can be imported
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

from mol_patcher.pdb_io import PdbParser
from mol_patcher.topology_io import TopologyParser
from mol_patcher.utilities import get_distance

def run_diagnostics():
    pdb_path = os.path.join(project_root, "pdbs", "m9.pdb")
    itp_path = os.path.join(project_root, "topology_files", "CARB-m9.itp")

    # 1. Load the raw files independently
    parser = PdbParser()
    _, pdb_records, _ = parser.read_file(pdb_path)
    moltype, itp_atoms, itp_bonds, _, _, _ = TopologyParser.read_file(itp_path)

    print(f"Loaded {len(pdb_records)} PDB atoms and {len(itp_atoms)} ITP atoms.")

    # 2. Map ITP definitions to PDB serial numbers
    # Assumes a 1:1 sequential alignment for the isolated patch
    itp_to_pdb = {atom.number: pdb_records[i] for i, atom in enumerate(itp_atoms)}

    itp_bond_pairs = set()
    for bond in itp_bonds:
        a1 = itp_to_pdb[bond.a1].serial
        a2 = itp_to_pdb[bond.a2].serial
        # Sort tuples so (1, 2) and (2, 1) evaluate as identical
        itp_bond_pairs.add(tuple(sorted((a1, a2))))

    print(f"ITP explicitly defines {len(itp_bond_pairs)} bonds.")

    # 3. Calculate distance-based geometric bonds
    dist_bond_pairs = set()
    for i in range(len(pdb_records)):
        for j in range(i + 1, len(pdb_records)):
            r1 = pdb_records[i]
            r2 = pdb_records[j]
            
            dist = get_distance((r1.x, r1.y, r1.z), (r2.x, r2.y, r2.z))
            is_h = r1.name.strip().startswith('H') or r2.name.strip().startswith('H')
            max_dist = 1.25 if is_h else 1.65
            
            if dist <= max_dist:
                dist_bond_pairs.add(tuple(sorted((r1.serial, r2.serial))))

    print(f"Distance calculation found {len(dist_bond_pairs)} physical bonds.")

    # 4. Compare the two sets
    missing_in_itp = dist_bond_pairs - itp_bond_pairs
    missing_in_dist = itp_bond_pairs - dist_bond_pairs

    print("\n--- Discrepancy Report ---")
    print(f"Bonds found by Distance but MISSING in ITP: {len(missing_in_itp)}")
    for a1, a2 in missing_in_itp:
        r1 = next(r for r in pdb_records if r.serial == a1)
        r2 = next(r for r in pdb_records if r.serial == a2)
        print(f"  -> {r1.res_name} {r1.res_seq} | {r1.name} -- {r2.name}")

    print(f"\nBonds defined in ITP but MISSING in Distance: {len(missing_in_dist)}")
    for a1, a2 in missing_in_dist:
        r1 = next(r for r in pdb_records if r.serial == a1)
        r2 = next(r for r in pdb_records if r.serial == a2)
        dist = get_distance((r1.x, r1.y, r1.z), (r2.x, r2.y, r2.z))
        print(f"  -> {r1.res_name} {r1.res_seq} | {r1.name} -- {r2.name} (Dist: {dist:.2f} A)")

if __name__ == "__main__":
    run_diagnostics()