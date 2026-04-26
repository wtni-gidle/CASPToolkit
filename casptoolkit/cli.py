import argparse
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent

COMMANDS = {
    "qa-af3": BASE / "CASP" / "qa_af3.py",
    "clashscore": BASE / "CASP" / "phenix_clashscore.py",
    "af3server-to-pro": BASE / "CASP" / "af3server_to_pro.py",
    "relax-prism": BASE / "CASP" / "relax_prism.py",
    "sup-assemble": BASE / "CASP" / "sup_assemble.py",
    "sup-oligo": BASE / "CASP" / "sup_homooligo.py",
    "sup-template": BASE / "CASP" / "sup_template.py",

    "cif2pdb": BASE / "PDBOps" / "cif2pdb.py",
    "merge-structure": BASE / "PDBOps" / "merge_structure.py",
    "reassign-chainid": BASE / "PDBOps" / "reassign_chain_id.py",
}


def main():
    parser = argparse.ArgumentParser(prog="casp", add_help=True)
    parser.add_argument("command", nargs="?", help="subcommand")

    args, rest = parser.parse_known_args()

    if args.command is None:
        parser.print_help()
        print("\nAvailable commands:")
        for k in COMMANDS:
            print(f"  {k}")
        sys.exit(0)

    if args.command not in COMMANDS:
        print(f"Unknown command: {args.command}")
        print("Available:", ", ".join(COMMANDS))
        sys.exit(1)

    script = COMMANDS[args.command]

    cmd = [sys.executable, str(script)] + rest

    result = subprocess.run(cmd)
    sys.exit(result.returncode)