"""Renumber atom serial numbers in PDB files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

from Bio import PDB

from casptoolkit.PDBOps._utils import sort_chains

LOGGER = logging.getLogger(__name__)


def renumber_atoms(
    structure: PDB.Structure.Structure,
    output_path: str,
    chain_order: Optional[List[str]] = None,
) -> None:
    """Renumber atom serial numbers, resetting the counter for each chain.

    Handles structures where atom serial numbers would exceed the PDB format
    limit of 99999 by writing each chain separately with independent numbering.

    Args:
        structure: BioPython Structure object.
        output_path: Path to the output PDB file.
        chain_order: Explicit chain output order. If None, letter chains
            precede digit chains, sorted lexicographically within each group.
    """
    out = Path(output_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    chains = sort_chains(structure[0], chain_order)

    pdb_io = PDB.PDBIO()
    with open(out, "w") as fh:
        for chain in chains:
            tmp = PDB.Structure.Structure("s")
            m = PDB.Model.Model(0)
            m.add(chain.copy())
            tmp.add(m)
            pdb_io.set_structure(tmp)
            pdb_io.save(fh, write_end=False, preserve_atom_numbering=False)
        fh.write("END\n")
