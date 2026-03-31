"""Merge multiple PDB structures into one."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

from Bio import PDB

from casptoolkit.PDBOps._utils import print_cli_settings
from casptoolkit.PDBOps.renumber_atoms import renumber_atoms

LOGGER = logging.getLogger(__name__)

_CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


def merge_pdb_files(
    input_files: List[str],
    output_file: str,
    renumber: bool = False,
) -> None:
    """Merge PDB files into a single structure.

    Chains are reassigned sequential IDs from the pool A-Z, a-z, 0-9.

    Args:
        input_files: Ordered list of input PDB file paths.
        output_file: Path to the output PDB file.
        renumber: If True, renumber atom serial numbers after merging.

    Raises:
        ValueError: If the total number of chains exceeds 62.
    """
    out = Path(output_file).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    parser = PDB.PDBParser(QUIET=True)
    structure = PDB.Structure.Structure("structure")
    model = PDB.Model.Model(0)
    structure.add(model)

    chain_id_pool = list(_CHAIN_IDS)
    for input_file in input_files:
        sub_model = parser.get_structure("s", input_file)[0]
        for chain in sub_model:
            if not chain_id_pool:
                raise ValueError("Too many chains: exceeded the 62-chain PDB limit.")
            new_chain = chain.copy()
            new_chain.id = chain_id_pool.pop(0)
            model.add(new_chain)

    if renumber:
        renumber_atoms(structure, str(out))
    else:
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(str(out))

    LOGGER.info("Merged %d files into %s", len(input_files), out)


def main(args) -> None:
    input_files = sorted(Path(args.input_dir).iterdir())
    input_files = [str(f) for f in input_files if f.suffix == ".pdb"]
    merge_pdb_files(input_files, args.output_file, args.renumber_atoms)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Merge multiple PDB files into one. This will reassign chain IDs automatically.")
    parser.add_argument("input_dir", help="Input directory containing PDB files.")
    parser.add_argument("output_file", help="Output merged PDB file path.")
    parser.add_argument("--renumber-atoms", action="store_true", help="Renumber atom serial numbers.")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
