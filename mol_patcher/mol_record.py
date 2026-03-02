from dataclasses import dataclass, field, replace
from typing import List, Optional

"""
For more information on PDB/ITP file formatting see the following:
ITP: https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html
PDB: https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf (p. 187)
"""

@dataclass
class PdbRecord:
    """Represents a single line in a PDB file."""
    pdb_name: str
    record_type: str
    serial: int
    name: str
    res_name: str
    chain: str
    res_seq: int
    x: float
    y: float
    z: float
    occupancy: float = 1.00
    temp_factor: float = 0.00
    seg_id: str = ""


@dataclass
class ItpAtom:
    number: int
    type: str
    res_n: int
    res: str
    atom: str
    cgnr: int
    charge: float
    mass: float

@dataclass
class ItpBond:
    a1: int
    a2: int
    type: int

@dataclass
class ItpPair:
    a1: int
    a2: int
    type: int

@dataclass
class ItpAngle:
    a1: int
    a2: int
    a3: int
    type: int

@dataclass
class ItpDih:
    a1: int
    a2: int
    a3: int
    a4: int
    type: int

@dataclass
class Mol:
    name: str
    records: List[PdbRecord] = field(default_factory=list)
    atoms: List[ItpAtom] = field(default_factory=list)
    bonds: List[ItpBond] = field(default_factory=list)
    pairs: List[ItpPair] = field(default_factory=list)
    angles: List[ItpAngle] = field(default_factory=list)
    dihs: List[ItpDih] = field(default_factory=list)
    
    def load_itp(self, itp_path: str):
        from .itp_io import ItpParser
        a, b, p, ang, d = ItpParser.read_file(itp_path)
        self.atoms, self.bonds, self.pairs, self.angles, self.dihs = a, b, p, ang, d

    def get_atom(self, res_seq: int, chain: str, res_name: str, atom_name: str) -> PdbRecord:
        for atom in self.records:
            if atom.res_seq == res_seq and atom.chain == chain and atom.res_name == res_name and atom.name == atom_name:
                return atom
        raise ValueError(f"Atom {atom_name} not found in {res_name} {res_seq}")
    
    def reindex(self):
        """
        Synchronizes atom numbers and all bonded interactions
        using the index_map lookup dictionary.
        """
        index_map = {}
        
        for i, atom in enumerate(self.atoms):
            new_nr = i + 1
            # Record: {Old ID: New ID}
            index_map[atom.number] = new_nr
            # Update the atom objects directly
            atom.number = new_nr
            atom.cgnr = new_nr  
        # Sync PDB serials to match ITP atom numbers
        for i, record in enumerate(self.records):
            record.serial = i + 1


        # Update Bonds (a1, a2)
        new_bonds = []
        for bond in self.bonds:
            try:
                # We explicitly look up the new names for a1 and a2
                new_a1 = index_map[bond.a1]
                new_a2 = index_map[bond.a2]
                new_bonds.append(replace(bond, a1=new_a1, a2=new_a2))
            except KeyError:
                # If an atom was deleted during surgery, skip this bond
                continue
        self.bonds = new_bonds

        # Update Pairs (a1, a2)
        new_pairs = []
        for pair in self.pairs:
            try:
                new_a1 = index_map[pair.a1]
                new_a2 = index_map[pair.a2]
                new_pairs.append(replace(pair, a1=new_a1, a2=new_a2))
            except KeyError:
                continue
        self.pairs = new_pairs

        # Update Angles (a1, a2, a3)
        new_angles = []
        for angle in self.angles:
            try:
                new_a1 = index_map[angle.a1]
                new_a2 = index_map[angle.a2]
                new_a3 = index_map[angle.a3]
                new_angles.append(replace(angle, a1=new_a1, a2=new_a2, a3=new_a3))
            except KeyError:
                continue
        self.angles = new_angles

        # Update Dihedrals (a1, a2, a3, a4)
        new_dihs = []
        for dih in self.dihs:
            try:
                new_a1 = index_map[dih.a1]
                new_a2 = index_map[dih.a2]
                new_a3 = index_map[dih.a3]
                new_a4 = index_map[dih.a4]
                new_dihs.append(replace(dih, a1=new_a1, a2=new_a2, a3=new_a3, a4=new_a4))
            except KeyError:
                continue
        self.dihs = new_dihs