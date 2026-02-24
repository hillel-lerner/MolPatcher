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
        Synchronizes everything: atom numbers, charge groups, PDB serials, 
        and all bonded interaction pointers.
        """
        index_map = {}
        
        # 1. Update Atoms, CGNR, and create the Map
        for i, atom in enumerate(self.atoms):
            new_nr = i + 1
            # Map the OLD number to the NEW number
            index_map[atom.number] = new_nr
            
            atom.number = new_nr
            # Setting cgnr to match the new atom number for consistency
            atom.cgnr = new_nr  
            
        # 2. Sync PDB Records serials
        for i, record in enumerate(self.records):
            record.serial = i + 1

        # 3. Update all bonded interactions using the index_map
        # This helper replaces the old indices with the new ones from our map
        def remap_list(interaction_list, attr_count):
            new_list = []
            for obj in interaction_list:
                try:
                    # Target only the relevant attributes (a1, a2, etc.)
                    attrs = ['a1', 'a2', 'a3', 'a4'][:attr_count]
                    updated_params = {attr: index_map[getattr(obj, attr)] for attr in attrs}
                    # This is where 'replace' is used
                    new_list.append(replace(obj, **updated_params))
                except KeyError:
                    # If an atom was deleted during surgery, the interaction is skipped
                    continue
            return new_list

        self.bonds = remap_list(self.bonds, 2)
        self.pairs = remap_list(self.pairs, 2)
        self.angles = remap_list(self.angles, 3)
        self.dihs = remap_list(self.dihs, 4)