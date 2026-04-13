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
    ins_code: str = " "


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
    """
    A standardized container for molecular topology and coordinate data.
    Aggregates PDB records and GROMACS ITP topology lists into a single accessible object.
    """
    name: str
    moltype_section: str = ""
    records: List[PdbRecord] = field(default_factory=list)
    atoms: List[ItpAtom] = field(default_factory=list)
    bonds: List[ItpBond] = field(default_factory=list)
    pairs: List[ItpPair] = field(default_factory=list)
    angles: List[ItpAngle] = field(default_factory=list)
    dihs: List[ItpDih] = field(default_factory=list)
    
    def load_itp(self, itp_path: str):
        """
        Parses a GROMACS ITP file and populates the molecule's topology lists.
        
        :param itp_path: (str) Absolute or relative path to the .itp file.
        :return: None. Updates internal lists.
        """
        from .topology_io import TopologyParser
        moltype, a, b, p, ang, d = TopologyParser.read_file(itp_path)
        self.moltype_section = moltype
        self.atoms, self.bonds, self.pairs, self.angles, self.dihs = a, b, p, ang, d

    def get_pdb_record(self, res_seq: int, chain: str, res_name: str, atom_name: str) -> PdbRecord:
        """
        Retrieves a specific PdbRecord from the molecule based on its identifiers.
        
        :param res_seq: (int) The residue sequence number.
        :param chain: (str) The chain identifier.
        :param res_name: (str) The 3-letter residue name.
        :param atom_name: (str) The specific atom name (e.g., 'CA', 'NZ').
        :return: (PdbRecord) The matching atom record.
        """
        for atom in self.records:
            if (atom.res_seq == res_seq and 
                atom.chain.strip() == chain.strip() and 
                atom.res_name.strip() == res_name.strip() and 
                atom.name.strip() == atom_name.strip()):
                return atom
        raise ValueError(f"Atom {atom_name} not found in {res_name} {res_seq}")