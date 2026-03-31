"""Convert mmCIF files to PDB format."""

import argparse
import logging
from multiprocessing import Pool
from pathlib import Path

from Bio import PDB

from casptoolkit.PDBOps._utils import print_cli_settings
from casptoolkit.PDBOps.renumber_atoms import renumber_atoms

LOGGER = logging.getLogger(__name__)


def cif_to_pdb(input_path: str, output_path: str, renumber: bool = False) -> None:
    """Convert a single CIF file to PDB."""
    input_file = Path(input_path)
    output_file = Path(output_path)

    if not input_file.is_file():
        raise FileNotFoundError(f"Input CIF file does not exist: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    parser = PDB.MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("structure", input_file.as_posix())
    except Exception as exc:
        raise ValueError(f"Failed to parse CIF file: {input_file}") from exc

    if renumber:
        renumber_atoms(structure, output_file.as_posix())
    else:
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_file.as_posix())


def cif_to_pdb_in_parallel(
    input_dir: str,
    output_dir: str,
    renumber: bool = False,
    n_cpu: int = 1,
) -> None:
    """Convert all CIF files in a directory to PDB files."""
    input_directory = Path(input_dir)
    output_directory = Path(output_dir)

    if not input_directory.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_directory}")

    output_directory.mkdir(parents=True, exist_ok=True)

    total_args = []
    for cif_file in input_directory.iterdir():
        if cif_file.suffix.lower() == ".cif":
            output_path = output_directory / f"{cif_file.stem}.pdb"
            total_args.append((cif_file.as_posix(), output_path.as_posix(), renumber))

    if not total_args:
        raise ValueError(f"No CIF files found in input directory: {input_directory}")

    with Pool(n_cpu) as pool:
        pool.starmap(cif_to_pdb, total_args)


def _resolve_single_output_path(input_path: str, output_path: str) -> str:
    """Resolve output path for single-file conversion.

    If output_path does not end with '.pdb', treat it as a directory and
    write <input_stem>.pdb under it.
    """
    input_file = Path(input_path)
    if output_path.lower().endswith(".pdb"):
        return Path(output_path).resolve().as_posix()

    output_dir = Path(output_path).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    return (output_dir / f"{input_file.stem}.pdb").as_posix()


def main(args) -> None:
    input_path = Path(args.input_path).resolve()
    output_path = Path(args.output_path).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    if input_path.is_file():
        resolved_output_path = _resolve_single_output_path(str(input_path), str(output_path))
        cif_to_pdb(str(input_path), resolved_output_path, args.renumber_atoms)
    else:
        if output_path.is_file():
            raise ValueError(
                "When input_path is a directory, output_path must be a directory path."
            )
        cif_to_pdb_in_parallel(
            str(input_path), str(output_path), args.renumber_atoms, args.num_workers
        )
    LOGGER.info("Done.")


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description='Convert CIF files to PDB files.')
    parser.add_argument('input_path', type=str, help='Path to the input CIF file or directory.')
    parser.add_argument('output_path', type=str, help='Path to the output PDB file or directory.')
    parser.add_argument('--renumber-atoms', action='store_true', help='Renumber atoms in the structure.')
    parser.add_argument('--num-workers', type=int, default=1, help='Number of worker processes.')
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
