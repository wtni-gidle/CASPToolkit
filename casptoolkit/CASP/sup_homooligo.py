"""Superpose homooligomers (An to Am)."""

from __future__ import annotations

import argparse
import logging
import shutil
import tempfile
from pathlib import Path
from typing import List, Optional

from casptoolkit.config import USALIGN_ENV_VAR, USALIGN_PATH
from casptoolkit.PDBOps._utils import print_cli_settings
from casptoolkit.CASP.sup_assemble import split_chains, sup_assemble

LOGGER = logging.getLogger(__name__)


def sup_homooligomers(
    source_file: str,
    target_file: str,
    output_dir: str,
    usalign_command: str,
    renumber_atoms: bool = False,
    extra_args: Optional[List[str]] = None,
    output_prefix: str = "sup_",
) -> None:
    """Superpose each chain of *source_file* onto all chains of *target_file*.

    Args:
        source_file: Source homooligomer PDB/CIF file.
        target_file: Target homooligomer PDB/CIF file.
        output_dir: Directory to write superposed structures.
        usalign_command: Path to the USalign executable.
        renumber_atoms: If True, renumber atoms in merged output.
        extra_args: Additional arguments forwarded to USalign.
        output_prefix: Prefix used for output file names.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        source_chain_dir = temp_path / "source_chains"
        target_chain_dir = temp_path / "target_chains"
        split_chains(source_file, source_chain_dir.as_posix())
        split_chains(target_file, target_chain_dir.as_posix())

        for i, source_chain in enumerate(source_chain_dir.iterdir()):
            if not source_chain.is_file() or source_chain.suffix.lower() != ".pdb":
                continue
            output_path = out_dir / f"{output_prefix}{i}.pdb"
            sup_assemble(
                source_chain.as_posix(),
                output_path.as_posix(),
                usalign_command,
                target_chain_dir.as_posix(),
                renumber_atoms=renumber_atoms,
                extra_args=extra_args,
            )

    LOGGER.info("Done.")


def main(args) -> None:
    usalign_command = shutil.which(USALIGN_PATH)
    if not usalign_command:
        raise FileNotFoundError(
            f"Cannot find USalign executable. "
            f"Set {USALIGN_ENV_VAR} or add '{USALIGN_PATH}' to PATH."
        )

    source_file = Path(args.source_file).resolve()
    target_file = Path(args.target_file).resolve()
    if not source_file.is_file():
        raise FileNotFoundError(f"Source file does not exist: {source_file}")
    if not target_file.is_file():
        raise FileNotFoundError(f"Target file does not exist: {target_file}")

    extra_args = args.extra_args.split() if args.extra_args else None

    sup_homooligomers(
        source_file.as_posix(),
        target_file.as_posix(),
        args.output_dir,
        usalign_command,
        renumber_atoms=args.renumber_atoms,
        extra_args=extra_args,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description=(
            "Superpose each chain of a source homooligomer onto all chains of a target "
            "homooligomer. Produces one merged structure per source chain. Used for RNA multimers."
        )
    )
    parser.add_argument("source_file", type=str, help="Source homooligomer PDB/CIF file.")
    parser.add_argument("target_file", type=str, help="Target homooligomer PDB/CIF file.")
    parser.add_argument("output_dir", type=str, help="Output directory for superposed structures.")
    parser.add_argument("--renumber_atoms", action="store_true", help="Renumber atoms in the merged structure.")
    parser.add_argument("--extra_args", type=str, default=None, help="Additional arguments for USalign, e.g. '-mm 1 -ter 0'.")
    parser.add_argument("--output_prefix", type=str, default="sup_", help="Prefix for output file names (default: sup_).")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
