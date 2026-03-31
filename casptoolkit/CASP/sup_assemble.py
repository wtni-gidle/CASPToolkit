"""Superpose and assemble structures using USalign."""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

from Bio import PDB

from casptoolkit.config import USALIGN_ENV_VAR, USALIGN_PATH
from casptoolkit.PDBOps._utils import print_cli_settings
from casptoolkit.PDBOps.merge_structure import merge_pdb_files

LOGGER = logging.getLogger(__name__)


def run_usalign(
    model: str,
    reference: str,
    output_prefix: str,
    usalign_command: str,
    extra_args: Optional[List[str]] = None,
) -> None:
    """Run USalign for one model-reference pair and write superposed output."""
    command = [usalign_command, model, reference, "-o", output_prefix]
    if extra_args:
        command.extend(extra_args)

    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise subprocess.SubprocessError(
            f"Error occurred while running USalign: {result.stderr.strip()}"
        )


def split_chains(input_file: str, output_dir: str) -> None:
    """Split a PDB/mmCIF structure into per-chain PDB files."""
    input_path = Path(input_file)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input structure file does not exist: {input_path}")

    suffix = input_path.suffix.lower()
    if suffix == ".pdb":
        parser = PDB.PDBParser(QUIET=True)
    elif suffix in (".cif", ".mmcif"):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Please provide a PDB or mmCIF file.")

    io = PDB.PDBIO()
    model = parser.get_structure("structure", input_path.as_posix())[0]
    for chain in model.get_chains():
        io.set_structure(chain)
        io.save((out_dir / f"chain_{chain.id}.pdb").as_posix())


def sup_assemble(
    source_file: str,
    output_path: str,
    usalign_command: str,
    target_dir: str,
    renumber_atoms: bool = False,
    extra_args: Optional[List[str]] = None,
) -> None:
    """Superpose source onto each target chain/file, then merge outputs.

    Args:
        source_file: Path to the source PDB file.
        output_path: Path to save the merged PDB file.
        usalign_command: Path to the USalign executable.
        target_dir: Directory containing target PDB files.
        renumber_atoms: If True, renumber atoms in the merged structure.
        extra_args: Additional arguments forwarded to USalign.
    """
    source_path = Path(source_file).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"source_file does not exist: {source_path}")

    target_dir_path = Path(target_dir).resolve()
    if not target_dir_path.is_dir():
        raise NotADirectoryError(f"target_dir does not exist: {target_dir_path}")

    out = Path(output_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        sup_structures: List[str] = []
        for target_path in sorted(target_dir_path.iterdir()):
            if not target_path.is_file() or target_path.suffix.lower() != ".pdb":
                continue
            output_prefix = temp_path / f"sup_{target_path.stem}"
            run_usalign(
                source_path.as_posix(),
                target_path.as_posix(),
                output_prefix.as_posix(),
                usalign_command,
                extra_args,
            )
            sup_structures.append(f"{output_prefix.as_posix()}.pdb")

        if not sup_structures:
            raise ValueError(f"No PDB files found in target directory: {target_dir_path}")

        merge_pdb_files(sup_structures, out.as_posix(), renumber=renumber_atoms)

    LOGGER.info("Done.")


def main(args) -> None:
    usalign_command = shutil.which(USALIGN_PATH)
    if not usalign_command:
        raise FileNotFoundError(
            f"Cannot find USalign executable. "
            f"Set {USALIGN_ENV_VAR} or add '{USALIGN_PATH}' to PATH."
        )

    extra_args = args.extra_args.split() if args.extra_args else None

    if args.target_file is not None:
        with tempfile.TemporaryDirectory() as temp_dir:
            split_dir = Path(temp_dir) / "target_chains"
            split_chains(args.target_file, split_dir.as_posix())
            sup_assemble(
                args.source_file,
                args.output_path,
                usalign_command,
                target_dir=split_dir.as_posix(),
                renumber_atoms=args.renumber_atoms,
                extra_args=extra_args,
            )
    else:
        sup_assemble(
            args.source_file,
            args.output_path,
            usalign_command,
            target_dir=args.target_dir,
            renumber_atoms=args.renumber_atoms,
            extra_args=extra_args,
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description=(
            "Superimpose a source structure onto target structures and assemble the results. "
            "Targets can be supplied as a directory of PDB files (--target-dir) or as a single "
            "structure file that will be split into chains automatically (--target-file)."
        )
    )
    parser.add_argument("source_file", type=str, help="Path to the source PDB file.")
    parser.add_argument("output_path", type=str, help="Path to save the merged PDB file.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--target-dir", type=str, help="Directory containing target PDB files.")
    group.add_argument("--target-file", type=str, help="Target structure file to split into chains first.")
    parser.add_argument("--renumber-atoms", action="store_true", help="Renumber atoms in the merged structure.")
    parser.add_argument("--extra-args", type=str, default=None, help="Additional arguments for USalign, e.g. '-mm 1 -ter 0'.")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
