import argparse
import json
import logging
import os
import re
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path
from typing import List

from casptoolkit.config import PHENIX_CLASHSCORE_ENV_VAR, PHENIX_CLASHSCORE_PATH

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def _validate_n_cpu(n_cpu: int) -> None:
    if n_cpu < 1:
        raise ValueError(f"n_cpu must be >= 1, got {n_cpu}.")


def _resolve_phenix_command() -> str:
    """Resolve phenix.clashscore command from configured path or PATH."""
    configured = PHENIX_CLASHSCORE_PATH

    if os.path.sep in configured:
        if os.path.isfile(configured) and os.access(configured, os.X_OK):
            return configured
        raise FileNotFoundError(
            "Configured phenix.clashscore path is not executable: "
            f"{configured}. Set {PHENIX_CLASHSCORE_ENV_VAR} to a valid executable."
        )

    resolved = shutil.which(configured)
    if resolved:
        return resolved

    raise FileNotFoundError(
        "Cannot find phenix.clashscore executable. "
        f"Set {PHENIX_CLASHSCORE_ENV_VAR} or add '{configured}' to PATH."
    )


def calc_clashscore(file_path: str, phenix_command: str):
    """Calculate clashscore for a single PDB file."""
    command = [
        phenix_command,
        file_path,
        "nuclear=True",
        "keep_hydrogens=True",
    ]
    result = subprocess.run(command, capture_output=True, text=True)

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
    LOGGER.info("Clashscore for %s: %.3f", file_path, clashscore)

    return clashscore


def wrapper(file_path: str, phenix_command: str):
    clashscore = calc_clashscore(file_path, phenix_command)
    return file_path, clashscore


def _collect_input_files(args: argparse.Namespace) -> List[str]:
    if args.file:
        file_path = Path(args.file).resolve()
        if not file_path.is_file():
            raise FileNotFoundError(f"Input file does not exist: {file_path}")
        return [file_path.as_posix()]

    if args.directory:
        input_dir = Path(args.directory).resolve()
        if not input_dir.is_dir():
            raise NotADirectoryError(f"Input directory does not exist: {input_dir}")
        files = sorted(str(p.resolve()) for p in input_dir.glob("*.pdb") if p.is_file())
        if not files:
            raise ValueError(f"No .pdb files found in directory: {input_dir}")
        return files

    list_path = Path(args.list).resolve()
    if not list_path.is_file():
        raise FileNotFoundError(f"Input list file does not exist: {list_path}")

    with list_path.open("r", encoding="utf-8") as file_list:
        files = [line.strip() for line in file_list if line.strip()]

    if not files:
        raise ValueError(f"Input list file is empty: {list_path}")

    missing_files = [file for file in files if not Path(file).is_file()]
    if missing_files:
        raise FileNotFoundError(
            "Some files listed in input list do not exist. "
            f"First missing file: {missing_files[0]}"
        )

    return [str(Path(file).resolve()) for file in files]


def process_in_parallel(file_list, output_path, n_cpu, phenix_command):
    _validate_n_cpu(n_cpu)
    tasks = [(file_path, phenix_command) for file_path in file_list]
    with Pool(n_cpu) as pool:
        results = pool.starmap(wrapper, tasks)
    results = dict(results)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=4)


def main(args):
    phenix_command = _resolve_phenix_command()
    files = _collect_input_files(args)

    output_path = os.path.abspath(args.output_path)
    dirname = os.path.dirname(output_path)
    os.makedirs(dirname, exist_ok=True)

    process_in_parallel(files, output_path, args.n_cpu, phenix_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate phenix clashscore.')
    parser.add_argument('-f', '--file', type=str, help='Single PDB file to process.')
    parser.add_argument('-d', '--directory', type=str, help='Directory containing PDB files.')
    parser.add_argument('-l', '--list', type=str, help='File containing list of PDB files.')
    parser.add_argument('output_path', type=str, help='Path to the output file.')
    parser.add_argument('--n_cpu', type=int, default=1, help='Number of CPUs to use for parallel processing.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)

    options = [args.file, args.directory, args.list]
    if options.count(None) != 2:
        raise ValueError("You must specify exactly one of --file, --directory, or --list.")

    main(args)
