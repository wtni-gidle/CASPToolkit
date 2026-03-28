"""
Superpose An to Am.
"""
import os
from Bio import PDB
import argparse
import tempfile
import logging

from casptoolkit.CASP.sup_assemble import sup_assemble

logging.basicConfig(level=logging.INFO)


def split_chains(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    file_extension = os.path.splitext(input_file)[1].lower()
    if file_extension == '.pdb':
        parser = PDB.PDBParser()
    elif file_extension == '.cif' or file_extension == '.mmcif':
        parser = PDB.MMCIFParser()
    else:
        raise ValueError("Unsupported file format. Please provide a PDB or mmCIF file.")
    io = PDB.PDBIO()
    
    model = parser.get_structure('structure', input_file)[0]
    
    for chain in model.get_chains():
        chain_id = chain.id
        chain_file = os.path.join(output_dir, f"chain_{chain_id}.pdb")
        io.set_structure(chain)
        io.save(chain_file)


def sup_homooligomers(source_file, target_file, output_dir, renumber = True, extra_args = None):
    os.makedirs(output_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        source_dir = os.path.join(temp_dir, "source_chains")
        target_dir = os.path.join(temp_dir, "target_chains")
        split_chains(source_file, source_dir)
        split_chains(target_file, target_dir)
        for i, source_chain in enumerate(os.listdir(source_dir)):
            source_chain = os.path.join(source_dir, source_chain)
            output_path = os.path.join(output_dir, f"sup_{i}.pdb")
            sup_assemble(source_chain, target_dir, output_path, renumber, extra_args)
    
    logging.info("Superimposition completed.")


def main(args):
    sup_homooligomers(
        args.source_file, 
        args.target_file, 
        args.output_dir,  
        not args.no_renumber,
        args.extra_args
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Superimposition between Homooligomers (An to Am). " 
        "Superimpose each chain from structure An onto all chains from Am to obtain a new structure. "
        "A total of n new structures can be obtained. pdb and cif formats are supported. "
    )
    parser.add_argument("source_file", help="Input PDB file containing Homooligomer")
    parser.add_argument("target_file", help="Input PDB file containing Homooligomer")
    parser.add_argument("output_dir", help="Output directory for combined structures")
    parser.add_argument('--no_renumber', action='store_true', help='Do not renumber atoms in the structure.')
    parser.add_argument('--extra_args', nargs='*', default=None, help='Additional arguments for USalign.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)
    
    main(args)
