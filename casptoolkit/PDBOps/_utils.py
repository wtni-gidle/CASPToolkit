"""Internal utilities shared across PDBOps modules."""

from __future__ import annotations

import argparse
from typing import List, Optional

from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model


def print_cli_settings(args) -> None:
    """Print user settings to stdout in a consistent format."""
    sep = "-" * 59
    print(sep, flush=True)
    print("User settings:", flush=True)
    width = max(len(k) for k in vars(args))
    for key, value in vars(args).items():
        print(f"  {key:<{width}} : {value}", flush=True)
    print(sep, flush=True)


def sort_chains(model: Model, chain_order: Optional[List[str]] = None) -> List[Chain]:
    """Return chains of *model* in the requested order.

    Args:
        model: BioPython Model object.
        chain_order: Explicit list of chain IDs defining output order.
            If None, letter chains precede digit chains; within each group
            chains are sorted lexicographically.
            All chain IDs present in the model must appear in *chain_order*,
            and all IDs in *chain_order* must exist in the model.

    Returns:
        Ordered list of Chain objects.

    Raises:
        ValueError: If *chain_order* does not match the chains in *model*.
    """
    if chain_order is not None:
        model_ids = {c.id for c in model}
        order_ids = set(chain_order)
        if model_ids != order_ids:
            missing = model_ids - order_ids
            extra = order_ids - model_ids
            mismatch_details = []
            if missing:
                mismatch_details.append(f"missing from chain_order: {sorted(missing)}")
            if extra:
                mismatch_details.append(f"not in model: {sorted(extra)}")
            raise ValueError("chain_order mismatch — " + "; ".join(mismatch_details))
        return sorted(model, key=lambda c: chain_order.index(c.id))
    return sorted(model, key=lambda c: (c.id.isdigit(), c.id))
