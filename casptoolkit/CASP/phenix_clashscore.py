"""Calculate phenix clashscore for PDB files."""

from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional

from casptoolkit.config import PHENIX_CLASHSCORE_ENV_VAR, PHENIX_CLASHSCORE_PATH
from casptoolkit.PDBOps._utils import print_cli_settings

LOGGER = logging.getLogger(__name__)


def calc_clashscore(file_path: str, phenix_command: str) -> Optional[float]:
    """Calculate clashscore for a single PDB file."""
    result = subprocess.run(
        [phenix_command, file_path, "nuclear=True", "keep_hydrogens=True"],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        LOGGER.error("Error processing %s: %s", file_path, result.stderr.strip())
        return None

    match = re.search(r"clashscore\s*=\s*([\d.]+)", result.stdout)
    if not match:
        LOGGER.error(
            "Unable to parse clashscore from output for %s. Raw output: %s",
            file_path,
            result.stdout.strip(),
        )
        return None

    clashscore = float(match.group(1))
    LOGGER.info("Clashscore for %s: %.2f", file_path, clashscore)
    return clashscore


def _calc_clashscore_worker(file_path: str, phenix_command: str) -> tuple[str, Optional[float]]:
    return Path(file_path).stem, calc_clashscore(file_path, phenix_command)


def calc_clashscore_in_parallel(
    file_list: List[str],
    phenix_command: str,
    output_path: Optional[str] = None,
    num_workers: int = 1,
) -> None:
    """Calculate clashscores for multiple PDB files.

    Args:
        file_list: List of PDB file paths.
        phenix_command: Path to the phenix.clashscore executable.
        output_path: If provided, write results to this JSON file.
        num_workers: Number of parallel worker processes.

    """
    tasks = [(fp, phenix_command) for fp in file_list]
    with Pool(num_workers) as pool:
        results: Dict[str, Optional[float]] = dict(pool.starmap(_calc_clashscore_worker, tasks))

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=4)


def main(args) -> None:
    phenix_command = shutil.which(PHENIX_CLASHSCORE_PATH)
    if not phenix_command:
        raise FileNotFoundError(
            "Cannot find phenix.clashscore executable. "
            f"Set {PHENIX_CLASHSCORE_ENV_VAR} or add '{PHENIX_CLASHSCORE_PATH}' to PATH."
        )

    input_path = Path(args.input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    output_path = str(Path(args.output_path).resolve()) if args.output_path else None

    if input_path.is_file():
        if output_path is not None:
            LOGGER.warning("--output-path is ignored when input_path is a file.")
        calc_clashscore(input_path.as_posix(), phenix_command)
    else:
        files = [p.as_posix() for p in input_path.glob("*.pdb") if p.is_file()]
        if not files:
            raise ValueError(f"No .pdb files found in directory: {input_path}")
        calc_clashscore_in_parallel(files, phenix_command, output_path, args.num_workers)

    LOGGER.info("Done.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Calculate phenix clashscore for PDB files.")
    parser.add_argument("input_path", type=str, help="Path to a PDB file or directory of PDB files.")
    parser.add_argument("--output-path", type=str, default=None, help="Path to the output JSON file. Effective only when input_path is a directory.")
    parser.add_argument("--num-workers", type=int, default=1, help="Number of worker processes (default: 1).")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
