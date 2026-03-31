"""Superpose models to a template and calculate TM-score."""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from casptoolkit.config import USALIGN_ENV_VAR, USALIGN_PATH
from casptoolkit.PDBOps._utils import print_cli_settings

LOGGER = logging.getLogger(__name__)


def run_usalign(
    model: str,
    reference: str,
    usalign_command: str,
    output_prefix: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
) -> Optional[float]:
    """Run USalign on a single model and return TM-score (normalised by reference)."""
    command = [usalign_command, model, reference]
    if output_prefix:
        command.extend(["-o", output_prefix])
    if extra_args:
        command.extend(extra_args)

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        LOGGER.error("Error processing %s: %s", model, result.stderr.strip())
        return None

    matches = re.findall(r"TM-score\s*=\s*([0-9.]+)", result.stdout)
    if len(matches) < 2:
        LOGGER.error("Error parsing TM-score for %s: %s", model, result.stdout.strip())
        return None

    tmscore = float(matches[1])
    LOGGER.info("TM-score for %s: %.4f", model, tmscore)
    return tmscore


def _run_usalign_worker(
    model: str,
    reference: str,
    usalign_command: str,
    output_prefix: Optional[str],
    extra_args: Optional[List[str]],
) -> tuple[str, Optional[float]]:
    return Path(model).stem, run_usalign(model, reference, usalign_command, output_prefix, extra_args)


def superpose_in_parallel(
    model_dir: str,
    reference: str,
    usalign_command: str,
    sup_dir: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    num_workers: int = 1,
) -> Dict[str, Optional[float]]:
    """Superpose all PDB/CIF files in *model_dir* to *reference* in parallel.

    Args:
        model_dir: Directory containing model PDB/CIF files.
        reference: Path to the reference PDB file.
        usalign_command: Path to the USalign executable.
        sup_dir: If provided, write superposed structures here.
        extra_args: Additional arguments forwarded to USalign.
        num_workers: Number of parallel worker processes.

    Returns:
        Mapping of model stem to TM-score (None if alignment failed).
    """
    models = [p.as_posix() for p in Path(model_dir).iterdir() if p.suffix.lower() in (".pdb", ".cif")]
    if not models:
        raise ValueError(f"No .pdb or .cif files found in model directory: {model_dir}")

    tasks = [
        (
            model,
            reference,
            usalign_command,
            str(Path(sup_dir) / f"{Path(model).stem}_sup") if sup_dir else None,
            extra_args,
        )
        for model in models
    ]

    with Pool(num_workers) as pool:
        results: Dict[str, Optional[float]] = dict(pool.starmap(_run_usalign_worker, tasks))

    LOGGER.info("Processed %d models.", len(results))
    return results


def main(args) -> None:
    usalign_command = shutil.which(USALIGN_PATH)
    if not usalign_command:
        raise FileNotFoundError(
            f"Cannot find USalign executable. "
            f"Set {USALIGN_ENV_VAR} or add '{USALIGN_PATH}' to PATH."
        )

    model_dir = Path(args.model_dir).resolve()
    reference = Path(args.reference).resolve()
    if not model_dir.is_dir():
        raise NotADirectoryError(f"Model directory does not exist: {model_dir}")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference file does not exist: {reference}")

    sup_dir: Optional[Path] = None
    if args.sup_dir:
        sup_dir = Path(args.sup_dir).resolve()
        sup_dir.mkdir(parents=True, exist_ok=True)

    extra_args = args.extra_args.split() if args.extra_args else None

    tmscore_dict = superpose_in_parallel(
        model_dir.as_posix(),
        reference.as_posix(),
        usalign_command,
        sup_dir.as_posix() if sup_dir else None,
        extra_args,
        args.num_workers,
    )

    if sup_dir:
        for pml in sup_dir.glob("*.pml"):
            pml.unlink()

    if args.output_file is None:
        return

    output_path = Path(args.output_file).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {"model": tmscore_dict.keys(), "tmscore": tmscore_dict.values()}
    ).sort_values(by="tmscore", ascending=False, na_position="last").to_csv(output_path, index=False)
    LOGGER.info("Results saved to %s", output_path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description="Superpose models to a template and calculate TM-score."
    )
    parser.add_argument("model_dir", type=str, help="Directory containing model PDB/CIF files.")
    parser.add_argument("reference", type=str, help="Reference PDB file.")
    parser.add_argument("--sup-dir", type=str, default=None, help="Directory to save superposed structures.")
    parser.add_argument("--output-file", type=str, default=None, help="Output CSV file for TM-score results.")
    parser.add_argument("--extra-args", type=str, default=None, help="Additional arguments for USalign.")
    parser.add_argument("--num-workers", type=int, default=1, help="Number of worker processes (default: 1).")
    args = parser.parse_args()

    if args.output_file is None and args.sup_dir is None:
        parser.error("Specify at least one of --output-file or --sup-dir.")

    print_cli_settings(args)
    main(args)
