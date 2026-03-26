"""
Superpose models to template and calculate tmscore.
"""
import os
import re
import subprocess
import pandas as pd
import logging
from multiprocessing import Pool
import argparse
from pathlib import Path

from PDBToolkit.config import USALIGN_PATH

logging.basicConfig(level=logging.INFO)


def run_usalign(model, reference, output_prefix = None, extra_args = None):
    command = [
        USALIGN_PATH,
        model, 
        reference
    ]
    if output_prefix:
        command.extend(["-o", output_prefix])
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f"Error processing {model}: {result.stderr.strip()}")
        return None
    
    tmscore_pattern = re.compile(r'TM-score\s*=\s*([0-9.]+)')
    try:
        tmscore = float(tmscore_pattern.findall(result.stdout)[1])
        logging.info(f"TM-score for {model}: {tmscore}")
        return tmscore
    except:
        logging.error(f"Error parsing TM-score for {model}: {result.stdout.strip()}")
        return None


def wrapper(model, reference, output_prefix, extra_args):
    tmscore = run_usalign(model, reference, output_prefix, extra_args)
    return model, tmscore


def process_in_parallel(model_dir, reference_file, sup_dir = None, extra_args = None, n_cpu = 1):
    model_list = [os.path.join(model_dir, model) for model in os.listdir(model_dir) if model.endswith('.pdb') or model.endswith('.cif')]
    reference_list = [reference_file for model in model_list]
    if sup_dir:
        output_prefix_list = [os.path.join(sup_dir, os.path.basename(model).replace(".pdb", "_sup")) for model in model_list]
    else:
        output_prefix_list = [None for model in model_list]
    extra_args_list = [extra_args for model in model_list]
    total_args = zip(model_list, reference_list, output_prefix_list, extra_args_list)
    
    with Pool(n_cpu) as pool:
        results = pool.starmap(wrapper, total_args)
    
    logging.info(f"Processed {len(results)} models.")
    
    return dict(results)


def main(args):
    model_dir = Path(args.model_dir)
    reference_file = Path(args.reference)
    if args.sup_dir:
        sup_dir = Path(args.sup_dir)
        sup_dir.mkdir(parents=True, exist_ok=True)
    else:
        sup_dir = None

    tmscore_dict = process_in_parallel(
        model_dir.resolve().as_posix(), 
        reference_file.resolve().as_posix(), 
        sup_dir.resolve().as_posix(), 
        args.extra_args.split(" ") if args.extra_args else None, 
        args.n_cpu
    )
    tmscore_dict = {os.path.basename(k): v for k, v in tmscore_dict.items()}
    tmscore_df = pd.DataFrame({"model": tmscore_dict.keys(), "tmscore": tmscore_dict.values()})

    if args.output_file is not None:
        output_path = Path(args.output_file)
    else:
        output_path = sup_dir / "tmscore.csv"
    tmscore_df.to_csv(output_path, index=False)
    logging.info(f"Results saved to {output_path}")
    
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
