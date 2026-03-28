"""
Merge structures into one.
"""
from typing import List
from Bio import PDB
import os
import argparse
import logging

from casptoolkit.PDBOps.renumber_atom import renumber_atom


logging.basicConfig(level=logging.INFO)


def merge_structures(input_files: List, output_file, renumber = True):
    output_file = os.path.abspath(output_file)
    directory = os.path.dirname(output_file)
    os.makedirs(directory, exist_ok=True)
    
    parser = PDB.PDBParser(QUIET=True)
    structure = PDB.Structure.Structure('structure')
    model = PDB.Model.Model(0)
    structure.add(model)
    chain_ids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')

    for input_file in input_files:
        sub_model = parser.get_structure('sub_model', input_file)[0]
        for chain in sub_model:
            new_chain = chain.copy()
            new_chain.id = chain_ids.pop(0)
            model.add(new_chain)

    if renumber:
        renumber_atom(structure, output_file)
    else:
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_file)
    
    logging.info(f"Successfully merged structures into {output_file}")

def main(args):
    input_files = [os.path.join(args.input_dir, file) for file in os.listdir(args.input_dir)]
    merge_structures(input_files, args.output_file, not args.no_renumber)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge multiple PDB files into one')
    parser.add_argument('input_dir', help='Input directory')
    parser.add_argument('output_file', help='Output merged PDB file')
    parser.add_argument('--no_renumber', action='store_true', help='Do not renumber atoms in the structure.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)
    
    main(args)
