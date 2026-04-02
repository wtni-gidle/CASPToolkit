"""Radially expand multi-chain structures to reduce steric clashes."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
from Bio import PDB

from casptoolkit.PDBOps._utils import print_cli_settings
from casptoolkit.PDBOps.renumber_atoms import renumber_atoms

LOGGER = logging.getLogger(__name__)


def _iter_standard_atoms(chain):
    for residue in chain:
        if residue.id[0] != " ":
            continue
        yield from residue


def calculate_structure_center(model) -> np.ndarray:
    center = np.zeros(3)
    num_atoms = 0

    for chain in model:
        for atom in _iter_standard_atoms(chain):
            center += atom.coord
            num_atoms += 1

    if num_atoms == 0:
        raise ValueError("No atoms found in structure when calculating center.")

    return center / num_atoms


def calculate_chain_center(chain) -> np.ndarray | None:
    center = np.zeros(3)
    num_atoms = 0

    for atom in _iter_standard_atoms(chain):
        center += atom.coord
        num_atoms += 1

    if num_atoms == 0:
        return None

    return center / num_atoms


def relax_prism(
    input_path: str,
    output_path: str,
    translation_distance: float,
    renumber: bool = True,
) -> None:
    input_file = Path(input_path).resolve()
    output_file = Path(output_path).resolve()

    if not input_file.is_file():
        raise FileNotFoundError(f"Input PDB file does not exist: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_file.as_posix())
    model = structure[0]

    structure_center = calculate_structure_center(model)

    for chain in model:
        chain_center = calculate_chain_center(chain)

        if chain_center is None:
            LOGGER.warning("Chain %s has no standard residues; skipped.", chain.id)
            continue

        direction = chain_center - structure_center
        norm = np.linalg.norm(direction)

        if norm < 1e-6:
            LOGGER.warning("Chain %s too close to structure center; skipped.", chain.id)
            continue

        translation_vector = (direction / norm) * translation_distance

        for atom in chain.get_atoms():
            atom.coord += translation_vector

    if renumber:
        renumber_atoms(structure, output_file.as_posix())
    else:
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_file.as_posix())


def main(args) -> None:
    relax_prism(
        args.input_path,
        args.output_path,
        args.distance,
        renumber=args.renumber_atoms,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description="Radially expand chains in a multi-chain PDB to reduce steric clashes."
    )
    parser.add_argument("input_path", type=str, help="Input PDB file path.")
    parser.add_argument("output_path", type=str, help="Output PDB file path.")
    parser.add_argument("distance", type=float, help="Radial translation distance for each chain.")
    parser.add_argument("--renumber-atoms", action="store_true", help="Renumber atoms in the merged structure.")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)