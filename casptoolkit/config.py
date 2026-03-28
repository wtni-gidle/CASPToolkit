"""Runtime configuration for external binaries.

This module keeps backward-compatible constant names while allowing users to
override tool paths via environment variables.
"""

from __future__ import annotations

import os


USALIGN_ENV_VAR = "CASPTOOLKIT_USALIGN_PATH"
PHENIX_CLASHSCORE_ENV_VAR = "CASPTOOLKIT_PHENIX_CLASHSCORE_PATH"

LEGACY_USALIGN_ENV_VAR = "PDBTOOLKIT_USALIGN_PATH"
LEGACY_PHENIX_CLASHSCORE_ENV_VAR = "PDBTOOLKIT_PHENIX_CLASHSCORE_PATH"

_DEFAULT_USALIGN_PATH = "/projects/bbgs/nwentao/tools/casp_tools/casp-rna/bins/us-align/USalign"
_DEFAULT_PHENIX_CLASHSCORE_PATH = "/projects/bbgs/nwentao/tools/casp_tools/casp-rna/bins/phenix/phenix.clashscore"


def _resolve_path(
    primary_env_var_name: str,
    legacy_env_var_name: str,
    default_path: str,
) -> str:
    """Return binary path from new env var, legacy env var, or default path."""
    return os.getenv(primary_env_var_name) or os.getenv(legacy_env_var_name, default_path)


USALIGN_PATH = _resolve_path(
    USALIGN_ENV_VAR,
    LEGACY_USALIGN_ENV_VAR,
    _DEFAULT_USALIGN_PATH,
)
PHENIX_CLASHSCORE_PATH = _resolve_path(
    PHENIX_CLASHSCORE_ENV_VAR,
    LEGACY_PHENIX_CLASHSCORE_ENV_VAR,
    _DEFAULT_PHENIX_CLASHSCORE_PATH,
)
