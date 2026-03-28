"""Superpose models to a template and calculate TM-score."""

import argparse
import logging
import os
import re
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from casptoolkit.config import USALIGN_ENV_VAR, USALIGN_PATH

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def _validate_n_cpu(n_cpu: int) -> None:
    if n_cpu < 1:
        raise ValueError(f"n_cpu must be >= 1, got {n_cpu}.")


def _resolve_usalign_command() -> str:
    """Resolve USalign command from configured path or PATH."""
    configured = USALIGN_PATH

    if os.path.sep in configured:
        if os.path.isfile(configured) and os.access(configured, os.X_OK):
            return configured
        raise FileNotFoundError(
            "Configured USalign path is not executable: "
            f"{configured}. Set {USALIGN_ENV_VAR} to a valid executable."
        )

    resolved = shutil.which(configured)
    if resolved:
        return resolved

    raise FileNotFoundError(
        f"Cannot find USalign executable. Set {USALIGN_ENV_VAR} or add '{configured}' to PATH."
    )


def run_usalign(
    usalign_command,
    model,
    reference,
    output_prefix=None,
    extra_args=None,
):
    """Run USalign and return TM-score parsed from command output."""
    command = [
        usalign_command,
        model,
        reference,
    ]
    if output_prefix:
        command.extend(["-o", output_prefix])
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        LOGGER.error("Error processing %s: %s", model, result.stderr.strip())
        return None

    tmscore_pattern = re.compile(r"TM-score\s*=\s*([0-9.]+)")
    matches = tmscore_pattern.findall(result.stdout)
    if len(matches) < 2:
        LOGGER.error("Error parsing TM-score for %s: %s", model, result.stdout.strip())
        return None

    tmscore = float(matches[1])
    LOGGER.info("TM-score for %s: %.4f", model, tmscore)
    return tmscore


def wrapper(usalign_command, model, reference, output_prefix, extra_args):
    tmscore = run_usalign(usalign_command, model, reference, output_prefix, extra_args)
    return model, tmscore


def process_in_parallel(
    usalign_command,
    model_dir,
    reference_file,
    sup_dir=None,
    extra_args=None,
    n_cpu=1,
):
    _validate_n_cpu(n_cpu)

    model_list = [
        os.path.join(model_dir, model)
        for model in sorted(os.listdir(model_dir))
        if model.lower().endswith((".pdb", ".cif"))
    ]
    if not model_list:
        raise ValueError(f"No .pdb or .cif files found in model directory: {model_dir}")

    reference_list = [reference_file for _ in model_list]
    if sup_dir:
        output_prefix_list = [
            os.path.join(sup_dir, f"{Path(model).stem}_sup")
            for model in model_list
        ]
    else:
        output_prefix_list = [None for _ in model_list]
    extra_args_list = [extra_args for _ in model_list]
    usalign_command_list = [usalign_command for _ in model_list]
    total_args = zip(
        usalign_command_list,
        model_list,
        reference_list,
        output_prefix_list,
        extra_args_list,
    )

    with Pool(n_cpu) as pool:
        results = pool.starmap(wrapper, total_args)

    LOGGER.info("Processed %d models.", len(results))

    return dict(results)


def main(args):
    model_dir = Path(args.model_dir)
    reference_file = Path(args.reference)
    if not model_dir.is_dir():
        raise NotADirectoryError(f"Model directory does not exist: {model_dir}")
    if not reference_file.is_file():
        raise FileNotFoundError(f"Reference file does not exist: {reference_file}")

    _validate_n_cpu(args.n_cpu)
    usalign_command = _resolve_usalign_command()

    if args.sup_dir:
        sup_dir = Path(args.sup_dir)
        sup_dir.mkdir(parents=True, exist_ok=True)
    else:
        sup_dir = None

    extra_args = args.extra_args.split() if args.extra_args else None

    tmscore_dict = process_in_parallel(
        usalign_command,
        model_dir.resolve().as_posix(), 
        reference_file.resolve().as_posix(), 
        sup_dir.resolve().as_posix() if sup_dir else None,
        extra_args,
        args.n_cpu
    )
    tmscore_dict = {
        os.path.basename(k): v
        for k, v in tmscore_dict.items()
        if v is not None
    }
    if not tmscore_dict:
        raise RuntimeError("No valid TM-score results were generated.")

    tmscore_df = pd.DataFrame({"model": tmscore_dict.keys(), "tmscore": tmscore_dict.values()})
    tmscore_df.sort_values(by="tmscore", ascending=False, inplace=True)

    if args.output_file is not None:
        output_path = Path(args.output_file)
    else:
        output_path = sup_dir / "tmscore.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmscore_df.to_csv(output_path, index=False)
    LOGGER.info("Results saved to %s", output_path)

    if sup_dir:
        for file in sup_dir.glob("*.pml"):
            os.remove(file)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Superpose models to template and calculate TM-score.")
    parser.add_argument("model_dir", help="Directory containing model files to process.")
    parser.add_argument("reference", help="Reference PDB file.")
    parser.add_argument("--sup_dir", help="Directory to save the output PDB files.")
    parser.add_argument("--output_file", help="Output file to save the TM-score results.")
    parser.add_argument('--extra_args', type=str, default=None, 
                        help='Additional arguments for USalign.')
    parser.add_argument("--n_cpu", type=int, default=1, help="Number of CPUs to use.")
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)
    
    options = [args.output_file, args.sup_dir]
    if options.count(None) == 2:
        raise ValueError("You must specify at least one of --output_file, --sup_dir.")

    main(args)
