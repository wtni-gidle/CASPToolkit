"""Runtime configuration for external binaries.

Override tool paths via environment variables:
  CASPTOOLKIT_USALIGN_PATH
  CASPTOOLKIT_PHENIX_CLASHSCORE_PATH
"""

from __future__ import annotations

import os


USALIGN_ENV_VAR = "CASPTOOLKIT_USALIGN_PATH"
PHENIX_CLASHSCORE_ENV_VAR = "CASPTOOLKIT_PHENIX_CLASHSCORE_PATH"

_DEFAULT_USALIGN_PATH = "xxx"
_DEFAULT_PHENIX_CLASHSCORE_PATH = "xxx"


def _resolve_path(env_var_name: str, default_path: str) -> str:
    """Return binary path from env var or default path."""
    return os.getenv(env_var_name, default_path)


USALIGN_PATH = _resolve_path(USALIGN_ENV_VAR, _DEFAULT_USALIGN_PATH)
PHENIX_CLASHSCORE_PATH = _resolve_path(PHENIX_CLASHSCORE_ENV_VAR, _DEFAULT_PHENIX_CLASHSCORE_PATH)
