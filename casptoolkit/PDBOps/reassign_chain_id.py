"""Reassign chain IDs in PDB files."""

from __future__ import annotations

import argparse
import logging
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional

from Bio import PDB

from casptoolkit.PDBOps._utils import print_cli_settings, sort_chains
from casptoolkit.PDBOps.renumber_atoms import renumber_atoms

LOGGER = logging.getLogger(__name__)


def reassign_chain_id(
    input_path: str,
    output_path: str,
    chain_map: Dict[str, str],
    chain_order: Optional[List[str]] = None,
    renumber: bool = True,
) -> None:
    """Reassign chain IDs of a PDB file according to *chain_map*.

    Args:
        input_path: Path to the input PDB file.
        output_path: Path to the output PDB file.
        chain_map: Mapping from original chain ID to new chain ID.
            Every chain ID present in the input structure must appear as a key.
            Values must be unique single-character IDs.
        chain_order: Explicit output chain order. If None, letter chains
            precede digit chains, sorted lexicographically within each group.
        renumber: If True, renumber atom serial numbers after reassignment.

    Raises:
        ValueError: If *chain_map* keys do not match the chains in the input,
            or if new chain IDs are not unique single characters.
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    parser = PDB.PDBParser(QUIET=True)
    model = parser.get_structure("structure", input_path)[0]

    source_chain_ids = {c.id for c in model}
    map_keys = set(chain_map)
    if source_chain_ids != map_keys:
        missing = source_chain_ids - map_keys
        extra = map_keys - source_chain_ids
        mismatch_details = []
        if missing:
            mismatch_details.append(f"chains not in chain_map: {sorted(missing)}")
        if extra:
            mismatch_details.append(f"chain_map keys not in model: {sorted(extra)}")
        raise ValueError("chain_map mismatch — " + "; ".join(mismatch_details))

    target_chain_ids = list(chain_map.values())
    if len(target_chain_ids) != len(set(target_chain_ids)):
        raise ValueError("chain_map contains duplicate target chain IDs.")
    if any(len(v) != 1 for v in target_chain_ids):
        raise ValueError("All target chain IDs in chain_map must be single characters.")

    new_structure = PDB.Structure.Structure("new_structure")
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)

    for chain in model:
        new_chain = chain.copy()
        new_chain.id = chain_map[chain.id]
        new_model.add(new_chain)

    if renumber:
        renumber_atoms(new_structure, str(out), chain_order=chain_order)
    else:
        ordered_chains = sort_chains(new_model, chain_order)
        for chain_id in [c.id for c in new_model.get_list()]:
            new_model.detach_child(chain_id)
        for chain in ordered_chains:
            new_model.add(chain)
        io = PDB.PDBIO()
        io.set_structure(new_structure)
        io.save(str(out))


def reassign_chain_id_in_parallel(
    input_dir: str,
    output_dir: str,
    chain_map: Dict[str, str],
    chain_order: Optional[List[str]] = None,
    renumber: bool = True,
    n_cpu: int = 1,
) -> None:
    """Reassign chain IDs of all PDB files in *input_dir* in parallel.

    Args:
        input_dir: Directory containing input PDB files.
        output_dir: Directory to write output PDB files.
        chain_map: Mapping from original chain ID to new chain ID.
        chain_order: Explicit output chain order.
        renumber: If True, renumber atom serial numbers after reassignment.
        n_cpu: Number of worker processes.
    """
    in_dir = Path(input_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    task_args = [
        (str(in_dir / f), str(out_dir / f), chain_map, chain_order, renumber)
        for f in sorted(in_dir.iterdir())
        if f.suffix == ".pdb"
    ]

    with Pool(n_cpu) as pool:
        pool.starmap(reassign_chain_id, task_args)


def _parse_chain_mapping(chain_mapping: str) -> Dict[str, str]:
    """Parse 'ABC:DEF' into {'A': 'D', 'B': 'E', 'C': 'F'}.

    Raises:
        ValueError: If the format is invalid or both sides differ in length.
    """
    if ":" not in chain_mapping:
        raise ValueError(
            f"Invalid chain_mapping '{chain_mapping}': expected format 'ABC:DEF'."
        )
    src_ids, dst_ids = chain_mapping.split(":", 1)
    if len(src_ids) != len(dst_ids):
        raise ValueError(
            f"chain_mapping '{chain_mapping}': both sides must have the same length "
            f"({len(src_ids)} vs {len(dst_ids)})."
        )

    if len(set(src_ids)) != len(src_ids) or len(set(dst_ids)) != len(dst_ids):
        raise ValueError(
            f"chain_mapping '{chain_mapping}': chain IDs must be unique on both sides."
        )

    return dict(zip(src_ids, dst_ids))


def main(args: argparse.Namespace) -> None:
    chain_map = _parse_chain_mapping(args.chain_mapping)
    chain_order = list(args.chain_order) if args.chain_order else None

    input_path = Path(args.input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    if input_path.is_dir():
        reassign_chain_id_in_parallel(
            str(input_path),
            args.output_path,
            chain_map,
            chain_order,
            args.renumber_atoms,
            args.num_workers,
        )
    else:
        reassign_chain_id(
            str(input_path), args.output_path, chain_map, chain_order, args.renumber_atoms
        )
    LOGGER.info("Done.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description="Modify chain names in a PDB file using a mapping of original "
        "chain IDs to new chain IDs."
    )
    parser.add_argument("input_path", type=str, help="Path to the input PDB file or directory.")
    parser.add_argument("output_path", type=str, help="Path to the output PDB file or directory.")
    parser.add_argument("chain_mapping", type=str, help="Chain ID mapping, e.g. 'ABC:DEF' maps A→D, B→E, C→F.")
    parser.add_argument("--renumber_atoms", action="store_true", help="Renumber atom serial numbers.")
    parser.add_argument("--chain_order", type=str, default=None, help="Output chain order.")
    parser.add_argument("--num_workers", type=int, default=1, help="Number of worker processes (default: 1).")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
