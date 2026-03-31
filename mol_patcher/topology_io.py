from .mol_record import ItpAtom, ItpBond, ItpAngle, ItpDih, ItpPair
from .topology_tools import reindex_topology

class TopologyParser:
    @staticmethod
    def read_file(filepath):
        atoms, bonds, pairs, angles, dihs = [], [], [], [], []
        current_section = None

        with open(filepath, 'r') as f:
            for line in f:
                stripped = line.strip()
                
                # Ignore blank lines and comments 
                if not stripped or stripped.startswith(';'):
                    continue

                # Resets state exactly when a new section starts
                if stripped.startswith('['):
                    current_section = stripped.strip('[] ').lower()
                    continue

                parts = line.split()
                
                # Parse data 
                if current_section == "atoms" and len(parts) >= 8 and parts[0].isdigit():
                    atoms.append(ItpAtom(
                        int(parts[0]), 
                        parts[1], 
                        int(parts[2]), 
                        parts[3], 
                        parts[4], 
                        int(parts[5]),
                        float(parts[6]), 
                        float(parts[7])
                    ))

                elif current_section == "bonds" and len(parts) >= 3:
                    bonds.append(ItpBond(int(parts[0]), int(parts[1]), int(parts[2])))

                elif current_section == "pairs" and len(parts) >= 3:
                    pairs.append(ItpPair(int(parts[0]), int(parts[1]), int(parts[2])))

                elif current_section == "angles" and len(parts) >= 4:
                    angles.append(ItpAngle(int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])))

                elif current_section == "dihedrals" and len(parts) >= 5:
                    dihs.append(ItpDih(int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])))

        return atoms, bonds, pairs, angles, dihs

class TopologyBuilder:
    def __init__(self, mol, filename, mol_name):
        """Initializes the builder with the molecule, destination path, and strict mol_name."""
        self.mol = mol
        self.filename = filename
        self.mol_name = mol_name
    
    def write_itp(self):
        """Writes a GROMACS ITP file after ensuring internal indices are synced."""
        # Trigger the Mol object's internal reindexing to sync a1, a2, etc.
        reindex_topology(self.mol)

        with open(self.filename, 'w') as f:
            # Molecule Header
            f.write("[ moleculetype ]\n; name       nrexcl\n")
            f.write(f"{self.mol_name:<10s} 3\n\n")

            # Atoms Section
            if self.mol.atoms:
                f.write("[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n")
                for atom in self.mol.atoms:
                    # Check if the atom object has a comment, default to empty string if not
                    comment = getattr(atom, 'comment', '') 
                    f.write(f"{atom.number:>6d} {atom.type.strip():>10s} {atom.res_n:>6d} {atom.res:>6s} "
                            f"{atom.atom:>6s} {atom.cgnr:>6d} {atom.charge:>10.4f} {atom.mass:>10.4f} {comment}\n")
                f.write("\n")

            # Bonds Section 
            if self.mol.bonds:
                f.write("[ bonds ]\n;  ai    aj  funct\n")
                for b in self.mol.bonds:
                    f.write(f"{b.a1:>7d} {b.a2:>7d} {b.type:>7d}\n")
                f.write("\n")

            # Pairs Section 
            if self.mol.pairs:
                f.write("[ pairs ]\n;  ai    aj  funct\n")
                for p in self.mol.pairs:
                    f.write(f"{p.a1:>7d} {p.a2:>7d} {p.type:>7d}\n")
                f.write("\n")

            # Angles Section
            if self.mol.angles:
                f.write("[ angles ]\n;  ai    aj    ak  funct\n")
                for a in self.mol.angles:
                    f.write(f"{a.a1:>7d} {a.a2:>7d} {a.a3:>7d} {a.type:>7d}\n")
                f.write("\n")

            # Dihedrals Section 
            if self.mol.dihs:
                f.write("[ dihedrals ]\n;  ai    aj    ak    al  funct\n")
                for d in self.mol.dihs:
                    f.write(f"{d.a1:>7d} {d.a2:>7d} {d.a3:>7d} {d.a4:>7d} {d.type:>7d}\n")
                f.write("\n")