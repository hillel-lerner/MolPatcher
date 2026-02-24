from .mol_record import ItpAtom, ItpBond, ItpAngle, ItpDih, ItpPair

class ItpParser:
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

class BuildItp:
    def __init__(self, mol, filename):
        """Initializes the builder with the molecule and destination path."""
        self.mol = mol
        self.filename = filename

    def write_itp(self):
        with open(self.filename, 'w') as f:
            f.write("[ moleculetype ]\n")
            f.write(f"{self.mol.name:<10s} 3\n\n")

            if self.mol.atoms:
                f.write("[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n")
                for atom in self.mol.atoms:
                    f.write(f"{atom.number:>6d} {atom.type:>10s} {atom.res_n:>6d} {atom.res:>6s} "
                            f"{atom.atom:>6s} {atom.cgnr:>6d} {atom.charge:>10.4f} {atom.mass:>10.4f}\n")
                f.write("\n")

            # Bonds
            if self.mol.bonds:
                f.write("[ bonds ]\n")
                for bond in self.mol.bonds:
                    f.write(f"{bond.a1:5d} {bond.a2:5d} {bond.type:5d}\n")

            # Pairs
            if self.mol.pairs:
                f.write("[ pairs ]\n")
                for pair in self.mol.pairs:
                    f.write(f"{pair.a1:5d} {pair.a2:5d} {pair.type:5d}\n")
            
            # Angles
            if self.mol.angles:
                f.write("[ angles ]\n")
                for angle in self.mol.angles:
                    f.write(f"{angle.a1:5d} {angle.a2:5d} {angle.a3:5d} {angle.type:5d}\n")

            # Dihedrals
            if self.mol.dihs:
                f.write("[ dihedrals ]\n")
                for dih in self.mol.dihs:
                    f.write(f"{dih.a1:5d} {dih.a2:5d} {dih.a3:5d} {dih.a4:5d} {dih.type:5d}\n")