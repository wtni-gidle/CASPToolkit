"""
Modify chain ids in PDB files.
"""
from Bio import PDB
import os
import argparse
import logging
from multiprocessing import Pool
from casptoolkit.PDBOps.renumber_atom import renumber_atom

logging.basicConfig(level=logging.INFO)


def sort_chains(model, chain_order = None):
    """
    Sort chains by given order.
    """
    if chain_order:
        sorted_chains = sorted(model, key=lambda chain: chain_order.index(chain.id))
    else:
        sorted_chains = sorted(model, key=lambda chain: (chain.id.isdigit(), chain.id))

    for chain in model:
        model.detach_child(chain.id)

    for chain in sorted_chains:
        model.add(chain)


def reassign_chain_id(input_path, output_path, chain_map, chain_order = None, renumber = True):
    """
    Reassign chain ids of a PDB file according to the given chain map.
    """
    directory = os.path.dirname(output_path)
    os.makedirs(directory, exist_ok=True)

    parser = PDB.PDBParser(QUIET=True)
    model = parser.get_structure('structure', input_path)[0]

    new_structure = PDB.Structure.Structure('new_structure')
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)

    for chain in model:
        new_chain_id = chain_map[chain.id]
        new_chain = PDB.Chain.Chain(new_chain_id)
        for residue in chain:
            new_chain.add(residue.copy())
        new_model.add(new_chain)

    if renumber:
        renumber_atom(new_structure, output_path, chain_order=chain_order)
    else:
        sort_chains(new_model, chain_order=chain_order)
        io = PDB.PDBIO()
        io.set_structure(new_structure)
        io.save(output_path)


def reassign_chain_id_in_parallel(input_dir, output_dir, chain_map, chain_order = None, renumber = True, n_cpu = 1):
    """
    Reassign chain ids of PDB files in a directory in parallel.
    """
    os.makedirs(output_dir, exist_ok=True)
    total_args = []
    for file in os.listdir(input_dir):
        if file.endswith('.pdb'):
            input_path = os.path.join(input_dir, file)
            output_path = os.path.join(output_dir, file)
            total_args.append((input_path, output_path, chain_map, chain_order, renumber))

    with Pool(n_cpu) as pool:
        pool.starmap(reassign_chain_id, total_args)


def main(args):
    chain_map = dict(zip(args.orig_chain_ids, args.new_chain_ids))
    chain_order = args.chain_order
    renumber = not args.no_renumber
    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_path)
    n_cpu = args.n_cpu

    if os.path.isdir(input_path):
        reassign_chain_id_in_parallel(input_path, output_path, chain_map, chain_order, renumber, n_cpu)
    else:
        reassign_chain_id(input_path, output_path, chain_map, chain_order, renumber)
    logging.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=
        "Modify chain names in a PDB file using a mapping of original chain IDs "
        "to new chain IDs."
    )
    parser.add_argument('input_path', type=str, help='Path to the input PDB file.')
    parser.add_argument('output_path', type=str, help='Path to save the modified PDB file.')
    parser.add_argument('orig_chain_ids', type=str, help='Original chain IDs.')
    parser.add_argument('new_chain_ids', type=str, help='New chain IDs.')
    parser.add_argument('--no_renumber', action='store_true', help='Do not renumber atoms in the structure.')
    parser.add_argument('--chain_order', type=str, default=None, help='Order of chains in the structure.')
    parser.add_argument('--n_cpu', type=int, default=1, help='Number of CPUs to use for parallel processing.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)

    main(args)
