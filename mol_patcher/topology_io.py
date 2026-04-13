import os
from .mol_record import ItpAtom, ItpBond, ItpAngle, ItpDih, ItpPair
from .topology_tools import reindex_topology

class TopologyParser:
    """
    Handles parsing and cleanup of GROMACS .itp topology files.
    """

    @staticmethod
    def clean_itp_file(itp_path, temp):
        """
        Scans the ITP file for non-integer residue designations (e.g., '52A'). 
        If found, generates a sanitized temporary copy with sequential integers to prevent parser crashes.
        
        :param itp_path: (str) Path to the input .itp file.
        :param temp: (str) Path to the output directory for the temporary sanitized file.
        :return: (str) Path to the clean .itp file (returns original path if no fix was needed).
        """

        needs_fix = False
        with open(itp_path, 'r') as f:
            in_atoms = False
            for line in f:
                stripped = line.strip()
                if stripped.startswith("[ atoms ]"):
                    in_atoms = True
                    continue
                elif stripped.startswith("[") and in_atoms:
                    break
                
                if in_atoms and stripped and not stripped.startswith(';'):
                    parts = stripped.split()
                    if len(parts) > 2 and not parts[2].isdigit():
                        needs_fix = True
                        break

        # If all residues are integers, proceed with the original file
        if not needs_fix:
            return itp_path

        print(f"Warning: Non-integer residue detected in {os.path.basename(itp_path)}. Cleaning topology before loading...")
        clean_itp_path = os.path.join(temp, "temp_sanitized.itp")
        
        with open(itp_path, 'r') as infile, open(clean_itp_path, 'w') as outfile:
            in_atoms = False
            res_counter = 1
            last_seen_raw = None
            
            for line in infile:
                if line.strip().startswith("[ atoms ]"):
                    in_atoms = True
                    outfile.write(line)
                    continue
                elif line.strip().startswith("[") and in_atoms:
                    in_atoms = False

                if in_atoms and line.strip() and not line.strip().startswith(';'):
                    parts = line.split()
                    raw_res = parts[2]
                    
                    base_int_str = "".join(filter(str.isdigit, raw_res))
                    base_int = int(base_int_str) if base_int_str else 0

                    if last_seen_raw is None:
                        last_seen_raw = raw_res
                        res_counter = 1
                    elif base_int == 1 and raw_res != last_seen_raw:
                        res_counter = 1
                        last_seen_raw = raw_res
                    elif raw_res != last_seen_raw:
                        res_counter += 1
                        last_seen_raw = raw_res

                    parts[2] = str(res_counter)
                    
                    formatted_line = f"{parts[0]:>6} {parts[1]:>10} {parts[2]:>6} {parts[3]:>6} {parts[4]:>6} {parts[5]:>6} {parts[6]:>12} {parts[7]:>12}"
                    if len(parts) > 8:
                        formatted_line += "   " + " ".join(parts[8:])
                    outfile.write(formatted_line + "\n")
                else:
                    outfile.write(line)
                    
        return clean_itp_path


    @staticmethod
    def read_file(filepath):
        """
        Reads a clean GROMACS .itp file and extracts topological parameters.
        
        :param filepath: (str) Path to the .itp file.
        :return: (tuple) Lists containing (ItpAtom, ItpBond, ItpPair, ItpAngle, ItpDih) objects.
        """

        atoms, bonds, pairs, angles, dihs = [], [], [], [], []
        current_section = None

        moltype_lines = []

        with open(filepath, 'r') as f:
            for line in f:
                stripped = line.strip()

                if stripped.startswith('['):
                    current_section = stripped.strip('[] ').lower()
                    if current_section == 'moleculetype':
                        moltype_lines.append(line)
                    continue

                if current_section == 'moleculetype':
                    moltype_lines.append(line)
                    continue
                
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

        moltype_section = "".join(moltype_lines)
        return moltype_section, atoms, bonds, pairs, angles, dihs

class TopologyBuilder:
    """
    Handles formatting and writing topology objects back into a GROMACS-compatible .itp file.
    """

    def __init__(self, mol, filename, mol_name):
        """
        Constructs the TopologyBuilder.
        
        :param mol: (Mol) The molecule object containing the topology to write.
        :param filename: (str) The output file path.
        :param mol_name: (str) The strict molecule name to be written in the [ moleculetype ] header.
        """

        self.mol = mol
        self.filename = filename
        self.mol_name = mol_name
    
    def write_itp(self):
        # Writes a GROMACS ITP file after ensuring internal indices are synced.
        # Trigger the Mol object's internal reindexing to sync a1, a2, etc.
        reindex_topology(self.mol)

        with open(self.filename, 'w') as f:
            # Molecule Header
            if self.mol.moltype_section:
                f.write(self.mol.moltype_section)
                if not self.mol.moltype_section.endswith('\n\n'):
                    f.write('\n')
            else:
                f.write("[ moleculetype ]\n; Name            nrexcl\n")
                f.write(f"{self.mol.name:<15} 3\n\n")

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