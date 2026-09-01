"""
Standardized containers for molecular topology and coordinate data.

For more information on PDB/ITP file formatting see:
- ITP: https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html
- PDB: https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf (p. 187)
"""

from dataclasses import dataclass, field
from typing import List


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
    """Represents a single atom entry in a GROMACS [ atoms ] directive."""

    number: int
    type: str
    res_n: int
    res: str
    atom: str
    cgnr: int
    charge: float
    mass: float
    comment: str = ""


@dataclass
class ItpBond:
    """Represents a connection in a GROMACS [ bonds ] directive."""

    a1: int
    a2: int
    type: int


@dataclass
class ItpPair:
    """Represents a 1-4 interaction in a GROMACS [ pairs ] directive."""

    a1: int
    a2: int
    type: int


@dataclass
class ItpAngle:
    """Represents an angle in a GROMACS [ angles ] directive."""

    a1: int
    a2: int
    a3: int
    type: int


@dataclass
class ItpDih:
    """Represents a proper or improper dihedral in a GROMACS [ dihedrals ] directive."""

    a1: int
    a2: int
    a3: int
    a4: int
    type: int


@dataclass
class Notes:
    note: str


@dataclass
class Mol:
    """
    A standardized container for molecular topology and coordinate data.

    Aggregates PdbRecord objects and GROMACS ITP topology directives into a single
    accessible structure for modification by the MolPatcher pipeline.
    """

    name: str
    moltype_section: str = ""
    records: List[PdbRecord] = field(default_factory=list)
    atoms: List[ItpAtom] = field(default_factory=list)
    bonds: List[ItpBond] = field(default_factory=list)
    pairs: List[ItpPair] = field(default_factory=list)
    angles: List[ItpAngle] = field(default_factory=list)
    dihs: List[ItpDih] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def load_itp(self, itp_path: str):
        """
        Parses a GROMACS ITP file and populates the molecule's topology lists.

        :param itp_path: Absolute or relative path to the .itp file.
        :type itp_path: str
        """
        from .topology_io import TopologyParser

        moltype, a, b, p, ang, d = TopologyParser.read_file(itp_path)
        self.moltype_section = moltype
        self.atoms, self.bonds, self.pairs, self.angles, self.dihs = a, b, p, ang, d

    def get_pdb_record(
        self, res_seq: int, chain: str, res_name: str, atom_name: str
    ) -> PdbRecord:
        """
        Retrieves a specific PdbRecord from the molecule based on its topological identifiers.

        :param res_seq: The residue sequence number.
        :type res_seq: int
        :param chain: The chain identifier.
        :type chain: str
        :param res_name: The 3-letter (or 4-letter) residue name.
        :type res_name: str
        :param atom_name: The specific atom name (e.g., 'CA', 'NZ').
        :type atom_name: str
        :raises ValueError: If the specified atom cannot be found in the record list.
        :return: The matching atom coordinate record.
        :rtype: PdbRecord
        """
        for atom in self.records:
            if (
                atom.res_seq == res_seq
                and atom.chain.strip() == chain.strip()
                and atom.res_name.strip() == res_name.strip()
                and atom.name.strip() == atom_name.strip()
            ):
                return atom
        raise ValueError(f"Atom {atom_name} not found in {res_name} {res_seq}")
