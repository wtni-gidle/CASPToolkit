"""
Calculate QA scores and rank for CASP models.
"""
import argparse
import json
import logging
import os
import shutil
import zipfile
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio import PDB

from casptoolkit.PDBOps.cif2pdb import cif_to_pdb_in_parallel

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def _validate_n_cpu(n_cpu: int) -> None:
    if n_cpu < 1:
        raise ValueError(f"n_cpu must be >= 1, got {n_cpu}.")


def extract_files(zip_file, target_dir):
    """Extract CIF and summary JSON files from a zip archive."""
    try:
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            for file in zip_ref.namelist():
                if file.endswith(".cif") or (
                    "summary" in file and file.endswith(".json")
                ):
                    zip_ref.extract(file, target_dir)
    except zipfile.BadZipFile as exc:
        raise ValueError(f"Invalid zip file: {zip_file}") from exc


def calc_qa(directory, only_ptm=False):
    """Build a QA table from AF3 summary JSON files."""
    records = []
    for file in sorted(os.listdir(directory)):
        name, ext = os.path.splitext(file)
        if ext == ".json" and "summary_confidences" in name:
            pdb_file = name.replace("summary_confidences", "model") + ".pdb"
            with open(os.path.join(directory, file), "r") as f:
                metrics = json.load(f)

            missing_keys = [k for k in ["iptm", "ptm", "has_clash"] if k not in metrics]
            if missing_keys:
                raise KeyError(
                    f"Missing keys {missing_keys} in metrics file: {os.path.join(directory, file)}"
                )

            iptm = float(metrics["iptm"])
            ptm = float(metrics["ptm"])
            has_clash = bool(metrics["has_clash"])
            qa = ptm if only_ptm else (iptm * 0.8 + ptm * 0.2)

            records.append(
                {
                    "file": pdb_file,
                    "iptm": iptm,
                    "ptm": ptm,
                    "has_clash": has_clash,
                    "qa": qa,
                }
            )

    if not records:
        raise ValueError(
            f"No summary_confidences JSON files were found in directory: {directory}"
        )

    data = pd.DataFrame(records)
    data.sort_values(by="qa", ascending=False, inplace=True)
    data["rank"] = [f"rank_{i}.pdb" for i in range(1, len(data) + 1)]

    return data


def calc_plddt(pdb_file):
    """Calculate mean pLDDT from PDB B-factors."""
    parser = PDB.PDBParser(QUIET=True)
    model = parser.get_structure('structure', pdb_file)[0]

    b_factors = [
        atom.bfactor 
        for chain in model
        for residue in chain
        for atom in residue
    ]

    return np.mean(b_factors) / 100

def calc_plddt_wrapper(file):
    plddt = calc_plddt(file)
    return file, plddt


def _format_output_table(data, only_ptm=False):
    """Format numeric columns for stable output presentation."""
    output = data.copy()
    output["qa"] = output["qa"].map("{:.3f}".format)
    output["ptm"] = output["ptm"].map("{:.2f}".format)
    if only_ptm:
        output.drop(columns=["iptm"], inplace=True)
    else:
        output["iptm"] = output["iptm"].map("{:.2f}".format)
    output["plddt"] = output["plddt"].map("{:.4f}".format)
    return output


def qa_pipeline(input_dir, output_dir, renumber=True, no_clash=False, only_ptm=False, n_cpu=1):
    _validate_n_cpu(n_cpu)

    # unzip
    for file in sorted(os.listdir(input_dir)):
        if file.endswith(".zip"):
            extract_files(os.path.join(input_dir, file), output_dir)

    # cif to pdb
    has_cif = any(file.lower().endswith(".cif") for file in os.listdir(output_dir))
    if has_cif:
        LOGGER.info("Converting CIF files to PDB")
        cif_to_pdb_in_parallel(output_dir, output_dir, renumber=renumber, n_cpu=n_cpu)
    else:
        LOGGER.info("No CIF files found. Skipping CIF to PDB conversion.")

    # qa
    data = calc_qa(output_dir, only_ptm=only_ptm)

    # plddt
    pdb_paths = [os.path.join(output_dir, f) for f in data["file"].to_list()]
    missing_pdb = [path for path in pdb_paths if not os.path.isfile(path)]
    if missing_pdb:
        raise FileNotFoundError(
            "PDB files referenced by summary JSON are missing after conversion. "
            f"First missing file: {missing_pdb[0]}"
        )

    with Pool(n_cpu) as pool:
        plddts = pool.map(calc_plddt_wrapper, pdb_paths)
    plddt_dict = {os.path.basename(file): plddt for file, plddt in plddts}
    data["plddt"] = data["file"].map(plddt_dict)

    # rank
    for pdb_file, rank in zip(data["file"], data["rank"]):
        shutil.copy(os.path.join(output_dir, pdb_file), os.path.join(output_dir, rank))

    # clash
    if no_clash:
        clash_files = data[data["has_clash"]]["rank"]
        for file in clash_files:
            os.remove(os.path.join(output_dir, file))
        data = data[~data["has_clash"]].copy()

    output_table = _format_output_table(data, only_ptm=only_ptm)
    output_table.to_csv(os.path.join(output_dir, "qa.csv"), index=False, sep="\t")


def main(args):
    _validate_n_cpu(args.n_cpu)
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.isdir(input_dir):
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")

    os.makedirs(output_dir, exist_ok=True)
    qa_pipeline(
        input_dir,
        output_dir,
        not args.no_renumber,
        args.no_clash,
        args.only_ptm,
        args.n_cpu,
    )
    LOGGER.info("QA calculation completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="The directory containing the zip files to be processed")
    parser.add_argument("output_dir", help="The directory where the processed files will be saved")
    parser.add_argument('--no-renumber', action='store_true', help='Do not renumber atoms in the structure.')
    parser.add_argument('--no-clash', action='store_true', help=
                        'Do not include structures with clashes in the final ranking.')
    parser.add_argument('--only-ptm', action='store_true', help='Only calculate the ptm score.')
    parser.add_argument('--n-cpu', type=int, default=1, help='Number of CPUs to use for processing.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)
    
    main(args)
