#!/bin/bash
# Invoke from any directory; participant files are written only under --work/GCS.
set -euo pipefail
PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PACKAGE_DIR}/.."
exec "${PHASING_PYTHON:-python}" -u -m phasing_and_gnomix.run "$@"
