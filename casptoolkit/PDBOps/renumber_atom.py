"""
Renumber atom serial numbers in PDB files.
"""
from Bio import PDB
import os
import tempfile


def save_chain_as_structure(chain, output_file):
    """
    Save a given chain as a new PDB structure file.
    """
    new_structure = PDB.Structure.Structure("structure")
    new_model = PDB.Model.Model(0)

    new_chain = chain.copy()
    new_model.add(new_chain)
    new_structure.add(new_model)

    pdb_io = PDB.PDBIO()
    pdb_io.set_structure(new_structure)
    # This will automatically renumber the atom serial numbers.
    pdb_io.save(output_file, preserve_atom_numbering=False)


def merge_files(input_dir, output_file):
    """
    Directly merge individual PDB files into a single PDB file.
    """
    pdb_files = sorted(os.listdir(input_dir))
    with open(output_file, 'w') as outfile:
        for file in pdb_files:
            path = os.path.join(input_dir, file)
            with open(path, 'r') as infile:
                for line in infile:
                    if not line.startswith("END"):
                        outfile.write(line)
        outfile.write("END\n")


def renumber_atom(structure, output_path, chain_order = None):
    """
    Renumber the atom serial numbers in the structure. Supports custom chain ordering.

    This function handles large structures where atom serial numbers might exceed 
    the PDB format limit of 100000. It renumbers the atoms sequentially within each chain,
    ensuring that the numbering does not exceed this limit.
    """
    output_path = os.path.abspath(output_path)
    directory = os.path.dirname(output_path)
    os.makedirs(directory, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        model = structure[0]
        if chain_order:
            sorted_chains = sorted(model, key=lambda chain: chain_order.index(chain.id))
        else:
            sorted_chains = sorted(model, key=lambda chain: (chain.id.isdigit(), chain.id))

        for i, chain in enumerate(sorted_chains, start=1):
            chain_file = os.path.join(temp_dir, f"{i:02}_{chain.id}.pdb")
            save_chain_as_structure(chain, chain_file)

        merge_files(temp_dir, output_path)
