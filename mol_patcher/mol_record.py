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
    name: str
    records: List[PdbRecord] = field(default_factory=list)
    atoms: List[ItpAtom] = field(default_factory=list)
    bonds: List[ItpBond] = field(default_factory=list)
    pairs: List[ItpPair] = field(default_factory=list)
    angles: List[ItpAngle] = field(default_factory=list)
    dihs: List[ItpDih] = field(default_factory=list)
    
    def load_itp(self, itp_path: str):
        from .topology_io import TopologyParser
        a, b, p, ang, d = TopologyParser.read_file(itp_path)
        self.atoms, self.bonds, self.pairs, self.angles, self.dihs = a, b, p, ang, d

    def load_pdb(self, res_seq: int, chain: str, res_name: str, atom_name: str) -> PdbRecord:
        for atom in self.records:
            if atom.res_seq == res_seq and atom.chain == chain and atom.res_name == res_name and atom.name == atom_name:
                return atom
        raise ValueError(f"Atom {atom_name} not found in {res_name} {res_seq}")