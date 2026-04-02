"""Convert AlphaFold3 Server zip files to AF3 Pro output directory format."""
import argparse
import json
import logging
import re
import zipfile
from pathlib import Path
import functools
from typing import Any, Sequence

import zstandard as zstd
from casptoolkit.PDBOps._utils import print_cli_settings

LOGGER = logging.getLogger(__name__)

# AF3 Server sequence key -> AF3 local type key
_SERVER_TO_LOCAL_KEY = {
    "proteinChain": "protein",
    "rnaSequence": "rna",
    "dnaSequence": "dna",
    "ion": "ligand",
    "ligand": "ligand",
}


@functools.lru_cache(maxsize=256)
def _int_to_chain_id(num: int) -> str:
    """Encodes a number as a string, using reverse spreadsheet style naming.

    Args:
        num: A positive integer.

    Returns:
        A string that encodes the positive integer using reverse spreadsheet style,
        naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
        usual way to encode chain IDs in mmCIF files.
    """
    num = num - 1  # 1-based indexing.
    output = []
    while num >= 0:
        output.append(chr(num % 26 + ord('A')))
        num = num // 26 - 1
    return ''.join(output)


def _get_chain_type(entry: dict[str, Any]):
    """Return (local_key, sub_dict) for a Server sequences entry."""
    for server_key, local_key in _SERVER_TO_LOCAL_KEY.items():
        if server_key in entry:
            return local_key, entry[server_key]
    raise ValueError(f"Unknown sequence type in job_request: {list(entry.keys())}")


def _expand_server_sequences(job_request: list[dict[str, Any]]):
    """Expand sequences from Server job_request, honouring the count field.

    Returns list of (chain_ids, local_key, sub_dict) in chain-assignment order.
    """
    chains = []
    next_chain_num = 1
    for entry in job_request[0]["sequences"]:
        local_key, sub_dict = _get_chain_type(entry)
        count = sub_dict.get("count", 1)
        chain_ids = []
        for _ in range(count):
            chain_ids.append(_int_to_chain_id(next_chain_num))
            next_chain_num += 1
        chains.append((chain_ids, local_key, sub_dict))
    return chains


def collect_zips(dirs: Sequence[str]):
    """Collect all .zip files from the given directories, sorted for determinism."""
    zips = []
    for d in dirs:
        p = Path(d)
        found = sorted(p.glob("*.zip"))
        LOGGER.info("Found %d zip(s) in %s", len(found), d)
        zips.extend(found)
    if not zips:
        raise ValueError("No zip files found in the provided directories.")
    return zips


def read_msa_first_sequence(a3m_bytes: bytes):
    """Extract the first sequence from A3M format bytes (strips insertion chars)."""
    lines = a3m_bytes.decode("utf-8").splitlines()
    seq_lines = []
    in_first = False
    for line in lines:
        if line.startswith(">"):
            if in_first:
                break
            in_first = True
        elif in_first:
            seq_lines.append(re.sub(r"[a-z]", "", line).strip())
    return "".join(seq_lines)


def compress_to_zst(data: bytes):
    """Compress bytes using zstandard."""
    cctx = zstd.ZstdCompressor()
    return cctx.compress(data)


def build_sequences_section(job_request: list[dict[str, Any]]):
    """Build initial Pro `sequences` from job_request and assign chain IDs by count."""
    expanded = _expand_server_sequences(job_request)
    sequences_out = []
    seq_to_chain_ids = {}          # sequence_str -> (group_name, local_key)
    group_to_entry = {}            # group_name -> nested sequence dict (protein/rna/dna/ligand payload)
    group_to_local_key = {}        # group_name -> local_key

    for chain_ids, local_key, sub_dict in expanded:
        group_name = "_".join(chain_ids)
        pro_id = chain_ids[0] if len(chain_ids) == 1 else chain_ids
        group_to_local_key[group_name] = local_key

        if local_key == "protein":
            seq = sub_dict.get("sequence", "")
            if seq in seq_to_chain_ids:
                raise ValueError(
                    "Duplicate sequences are not supported in this converter; "
                    f"found duplicate protein sequence for groups "
                    f"'{seq_to_chain_ids[seq][0]}' and '{group_name}'."
                )
            modifications = [
                {"ptmType": m["ptmType"].removeprefix("CCD_"), "ptmPosition": m["ptmPosition"]}
                for m in sub_dict.get("modifications", [])
            ]
            protein_entry = {
                "id": pro_id,
                "sequence": seq,
                "modifications": modifications,
                "templates": [],
            }
            sequences_out.append({"protein": protein_entry})
            group_to_entry[group_name] = protein_entry
            seq_to_chain_ids[seq] = (group_name, local_key)

        elif local_key == "rna":
            seq = sub_dict.get("sequence", "")
            if seq in seq_to_chain_ids:
                raise ValueError(
                    "Duplicate sequences are not supported in this converter; "
                    f"found duplicate RNA sequence for groups "
                    f"'{seq_to_chain_ids[seq][0]}' and '{group_name}'."
                )
            modifications = [
                {"modificationType": m["modificationType"].removeprefix("CCD_"), "basePosition": m["basePosition"]}
                for m in sub_dict.get("modifications", [])
            ]
            rna_entry = {
                "id": pro_id,
                "sequence": seq,
                "modifications": modifications,
                "unpairedMsa": "",
            }
            sequences_out.append({"rna": rna_entry})
            group_to_entry[group_name] = rna_entry
            seq_to_chain_ids[seq] = (group_name, local_key)

        elif local_key == "dna":
            modifications = [
                {"modificationType": m["modificationType"].removeprefix("CCD_"), "basePosition": m["basePosition"]}
                for m in sub_dict.get("modifications", [])
            ]
            dna_entry = {
                "id": pro_id,
                "sequence": sub_dict.get("sequence", ""),
                "modifications": modifications,
            }
            sequences_out.append({"dna": dna_entry})
            group_to_entry[group_name] = dna_entry

        elif local_key == "ligand":
            if "ion" in sub_dict:
                ccd_codes = [sub_dict["ion"]]
            elif "ligand" in sub_dict:
                ccd_codes = [sub_dict["ligand"].removeprefix("CCD_")]
            else:
                raise ValueError(f"Cannot extract ligand CCD code from: {sub_dict}")
            ligand_entry = {
                "id": pro_id,
                "ccdCodes": ccd_codes,
            }
            sequences_out.append({"ligand": ligand_entry})
            group_to_entry[group_name] = ligand_entry

    return sequences_out, seq_to_chain_ids, group_to_entry, group_to_local_key


def enrich_sequences_with_paths(
    group_to_entry: dict[str, dict[str, Any]],
    group_to_local_key: dict[str, str],
    query_to_hit_map: dict[str, list[dict[str, Any]]],
    target_name: str,
    output_dir: str | Path,
):
    """Attach resolved MSA/template path fields to the prebuilt sequence entries."""
    output_dir = Path(output_dir)

    for group_name, local_key in group_to_local_key.items():
        entry = group_to_entry[group_name]

        if local_key == "protein":
            unpaired_name = f"{target_name}__{group_name}_unpairedmsa.a3m.zst"
            if (output_dir / unpaired_name).exists():
                entry["unpairedMsaPath"] = unpaired_name

            paired_name = f"{target_name}__{group_name}_pairedmsa.a3m.zst"
            if (output_dir / paired_name).exists():
                entry["pairedMsaPath"] = paired_name

            templates = []
            for n, hit in enumerate(query_to_hit_map.get(group_name, [])):
                mmcif_name = f"{target_name}__{group_name}_template_{n}.cif.zst"
                if (output_dir / mmcif_name).exists():
                    templates.append({
                        "queryIndices": hit["queryIndices"],
                        "templateIndices": hit["templateIndices"],
                        "mmcifPath": mmcif_name,
                    })
            entry["templates"] = templates

        elif local_key == "rna":
            entry.setdefault("unpairedMsa", "")
            unpaired_name = f"{target_name}__{group_name}_unpairedmsa.a3m.zst"
            if (output_dir / unpaired_name).exists():
                entry["unpairedMsaPath"] = unpaired_name


def process_zips(zips: Sequence[Path], target_name: str, output_dir: str | Path):
    """Extract and reorganise all zips into AF3 Pro directory format.

    Flow:
        1) Use only the first zip to build Pro JSON skeleton and MSA/template mapping.
        2) Extract model/summary/full_data from every zip.
        3) Emit final _data.json with resolved paths and template mappings.
    """
    target_dir = Path(output_dir)
    models_dir = target_dir / "models"
    summary_dir = target_dir / "summary_confidences"
    full_data_dir = target_dir / "full_data"
    for d in (models_dir, summary_dir, full_data_dir):
        d.mkdir(parents=True, exist_ok=True)

    msa_written = set()
    job_request = None
    seq_to_chain_ids = {}           # sequence_str -> (group_name, local_key)
    zip_token_to_group = {}         # zip token (e.g. a_b) -> group_name (from unpaired MSA)
    query_to_hit_map = {}           # pro group_name -> list of hit dicts
    data_json = None
    group_to_entry = {}
    group_to_local_key = {}

    model_re = re.compile(r".*_model_(\d+)\.cif$")
    summary_re = re.compile(r".*_summary_confidences_(\d+)\.json$")
    full_re = re.compile(r".*_full_data_(\d+)\.json$")
    msa_paired_re = re.compile(r"msas/.*_paired_msa_chains_([a-z_]+)\.a3m$")
    msa_unpaired_re = re.compile(r"msas/.*_unpaired_msa_chains_([a-z_]+)\.a3m$")
    tmpl_re = re.compile(r"templates/.*_template_hit_(\d+)_chains_([a-z_]+)\.cif$")
    qth_re = re.compile(r"templates/.*_template_hits_chains_([a-z_]+)_query_to_hit\.json$")

    first_zip_path = zips[0]

    for zip_path in zips:
        LOGGER.info("Processing %s", Path(zip_path).name)
        try:
            zf = zipfile.ZipFile(zip_path, "r")
        except zipfile.BadZipFile as exc:
            raise ValueError(f"Invalid zip file: {zip_path}") from exc

        with zf:
            names = set(zf.namelist())

            req_candidates = sorted(n for n in names if n.endswith("_job_request.json"))
            if not req_candidates:
                if zip_path == first_zip_path:
                    raise ValueError(
                        f"The first zip file '{Path(zip_path).name}' has no job_request.json; "
                        "cannot build output JSON."
                    )
                LOGGER.warning(
                    "No job_request.json found in %s, skipping model/summary/full_data extraction for this zip.",
                    Path(zip_path).name,
                )
                continue
            this_request = json.loads(zf.read(req_candidates[0]))
            if zip_path == first_zip_path:
                job_request = this_request
                sequences_out, seq_to_chain_ids, group_to_entry, group_to_local_key = build_sequences_section(job_request)
                data_json = {
                    "dialect": "alphafold3",
                    "version": 3,
                    "name": target_name,
                    "sequences": sequences_out,
                    "modelSeeds": job_request[0].get("modelSeeds", []),
                    "bondedAtomPairs": None,
                    "userCCD": None,
                }
                LOGGER.info("Initial Pro skeleton sequences: %d", len(sequences_out))

            seed = str(this_request[0]["modelSeeds"][0])

            # ---- model / summary / full_data ----
            for name in names:
                m = model_re.match(name)
                if m:
                    (models_dir / f"seed-{seed}_sample-{m.group(1)}_model.cif").write_bytes(zf.read(name))
                    continue
                m = summary_re.match(name)
                if m:
                    (summary_dir / f"seed-{seed}_sample-{m.group(1)}_summary_confidences.json").write_bytes(zf.read(name))
                    continue
                m = full_re.match(name)
                if m:
                    (full_data_dir / f"seed-{seed}_sample-{m.group(1)}_full_data.json").write_bytes(zf.read(name))

            # MSA/template resolution is defined only by the first zip.
            if zip_path != first_zip_path:
                continue

            # ---- Pass 1: use unpaired MSA to build token -> group mapping ----
            for name in sorted(names):
                m = msa_unpaired_re.match(name)
                if m:
                    ch = m.group(1)
                    if f"unpaired_{ch}" not in msa_written:
                        raw = zf.read(name)
                        msa_seq = read_msa_first_sequence(raw)
                        match = seq_to_chain_ids.get(msa_seq)
                        if not match:
                            raise ValueError(
                                f"Unpaired MSA for zip chain '{ch}' could not be matched to any "
                                f"sequence in job_request. First 50 chars: {msa_seq[:50]}"
                            )
                        group_name, local_key = match
                        existing = zip_token_to_group.get(ch)
                        if existing is not None and existing != group_name:
                            raise ValueError(
                                f"Zip chain token '{ch}' mapped inconsistently: {existing} vs {group_name}"
                            )
                        zip_token_to_group[ch] = group_name

                        compressed = compress_to_zst(raw)
                        if local_key in ("protein", "rna"):
                            (target_dir / f"{target_name}__{group_name}_unpairedmsa.a3m.zst").write_bytes(compressed)
                        msa_written.add(f"unpaired_{ch}")

            # ---- Pass 2: paired MSA files (reuse token -> group mapping) ----
            for name in sorted(names):
                m = msa_paired_re.match(name)
                if m:
                    ch = m.group(1)
                    if f"paired_{ch}" not in msa_written:
                        group_name = zip_token_to_group.get(ch)
                        if group_name is None:
                            raise ValueError(
                                f"Paired MSA for zip chain '{ch}' has no mapping from unpaired MSA. "
                                "Please ensure unpaired_msa_chains_* exists for this token."
                            )
                        local_key = group_to_local_key.get(group_name)
                        if local_key != "protein":
                            raise ValueError(
                                f"Paired MSA for zip chain '{ch}' mapped to non-protein group '{group_name}'."
                            )

                        raw = zf.read(name)
                        compressed = compress_to_zst(raw)
                        (target_dir / f"{target_name}__{group_name}_pairedmsa.a3m.zst").write_bytes(compressed)
                        msa_written.add(f"paired_{ch}")

            # ---- Pass 3: template files (reuse token -> group mapping) ----
            for name in sorted(names):
                m = tmpl_re.match(name)
                if m:
                    n_str, ch = m.group(1), m.group(2)
                    group_name = zip_token_to_group.get(ch)
                    if group_name is None:
                        raise ValueError(
                            f"Template for zip chain '{ch}' has no mapping from unpaired MSA. "
                            "Please ensure unpaired_msa_chains_* exists for this token."
                        )
                    local_key = group_to_local_key.get(group_name)
                    if local_key != "protein":
                        raise ValueError(
                            f"Template for zip chain '{ch}' mapped to non-protein group '{group_name}'."
                        )

                    key = f"tmpl_{group_name}_{n_str}"
                    if key not in msa_written:
                        (target_dir / f"{target_name}__{group_name}_template_{n_str}.cif.zst").write_bytes(
                            compress_to_zst(zf.read(name))
                        )
                        msa_written.add(key)
                    continue

                m = qth_re.match(name)
                if m:
                    ch = m.group(1)
                    group_name = zip_token_to_group.get(ch)
                    if group_name is None:
                        raise ValueError(
                            f"Template hit map for zip chain '{ch}' has no mapping from unpaired MSA. "
                            "Please ensure unpaired_msa_chains_* exists for this token."
                        )
                    local_key = group_to_local_key.get(group_name)
                    if local_key != "protein":
                        raise ValueError(
                            f"Template hit map for zip chain '{ch}' mapped to non-protein group '{group_name}'."
                        )

                    hits = json.loads(zf.read(name))
                    if group_name not in query_to_hit_map:
                        query_to_hit_map[group_name] = hits

    if job_request is None:
        raise ValueError(
            f"No valid zip files were processed for target_name '{target_name}'. "
            "Ensure the directory contains at least one valid AF3 Server zip with a job_request file."
        )

    enrich_sequences_with_paths(
        group_to_entry=group_to_entry,
        group_to_local_key=group_to_local_key,
        query_to_hit_map=query_to_hit_map,
        target_name=target_name,
        output_dir=target_dir,
    )
    data_path = target_dir / f"{target_name}_data.json"
    data_path.write_text(json.dumps(data_json, indent=2))
    LOGGER.info("Wrote %s", data_path)
    LOGGER.info("Done. Output directory: %s", target_dir)


def main(args):
    zips = collect_zips(args.dirs)
    LOGGER.info("Total zip files to process: %d", len(zips))
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    process_zips(zips, args.target_name, args.output_dir)
    LOGGER.info("Conversion completed.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Convert AF3 Server zip files to AF3 Pro output directory format.")
    parser.add_argument("target_name", help="Name for the output target and file prefixes")
    parser.add_argument("dirs", nargs="+", help="Input directories containing AF3 Server zip files")
    parser.add_argument("-o", "--output-dir", required=True, help="Output parent directory")
    args = parser.parse_args()

    print_cli_settings(args)
    main(args)
