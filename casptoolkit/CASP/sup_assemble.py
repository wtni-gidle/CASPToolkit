"""
Superpose and assemble.
"""
import os
import subprocess
import tempfile
import argparse
import logging

from casptoolkit.PDBOps.merge_structure import merge_pdb_files
from casptoolkit.config import USALIGN_PATH

logging.basicConfig(level=logging.INFO)


def run_usalign(model, reference, output_prefix, extra_args = None):
    command = [
        USALIGN_PATH,
        model, 
        reference, 
        "-o", output_prefix
    ]
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise subprocess.SubprocessError(f"Error occurred while running USalign: {result.stderr.strip()}")


def sup_assemble(source_file, target_dir, output_path, renumber = True, extra_args = None):
    output_path = os.path.abspath(output_path)
    dirname = os.path.dirname(output_path)
    os.makedirs(dirname, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp_dir:
        sup_structures = []
        for filename in os.listdir(target_dir):
            if filename.endswith('.pdb'):
                target_file = os.path.join(target_dir, filename)
                output_prefix = os.path.join(temp_dir, f"sup_{filename[:-4]}")
                run_usalign(source_file, target_file, output_prefix, extra_args)
                sup_structures.append(output_prefix + ".pdb")

        merge_pdb_files(sup_structures, output_path, renumber=renumber)


def main(args):
    sup_assemble(
        args.source_file, 
        args.target_dir, 
        args.output_path, 
        not args.no_renumber,
        args.extra_args
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Superimpose source structure onto target structures listed in a directory "
        "and then assemble the resulted structures. When predicting a large complex composed of "
        "multiple identical subunits, you might first predict a small complex composed of truncated "
        "subunits, each of which is a truncated version of the original subunit. "
        "Then, each complete subunit can be superimposed onto the truncated subunits one by one "
        "and then assemble resulted structures into the final structure. "
        "You can use this program for this process. "
        "Before using the program, ensure that your small complex is split into individual subunits."
    )
    parser.add_argument('source_file', type=str, help="Path to the structure A PDB file.")
    parser.add_argument('target_dir', type=str, help="Directory containing structure B PDB files.")
    parser.add_argument('output_path', type=str, help="Path to save the merged PDB file.")
    parser.add_argument('--no_renumber', action='store_true', help='Do not renumber atoms in the structure.')
    parser.add_argument('--extra_args', nargs='*', default=None, help='Additional arguments for USalign.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)

    main(args)
