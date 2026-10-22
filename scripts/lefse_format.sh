#!/usr/bin/env bash
set -euo pipefail
unset R_HOME R_LIBS_USER R_LIBS_SITE
if command -v conda >/dev/null 2>&1; then
  exec conda run --no-capture-output -n "${PAGLEBA_LEFSE_ENV:-pagleba-lefse}" lefse_format_input.py "$@"
else
  exec micromamba run -n "${PAGLEBA_LEFSE_ENV:-pagleba-lefse}" lefse_format_input.py "$@"
fi
