#!/usr/bin/env bash
# Prepare a Colab session for one register entry (colab/README.md). Run from the repository root:
#   !bash colab/setup.sh            # the analysis stack + this package
#   !bash colab/setup.sh acquisition   # + ASF / HyP3 / STAC clients (needs Earthdata credentials)
# Heavy SAR tools (Dolphin, COMPASS, ISCE2) are installed by the launcher notebook that needs them, not here.
set -euo pipefail
git config user.email "${GIT_USER_EMAIL:-colab@localhost}"
git config user.name "${GIT_USER_NAME:-colab}"
pip install -q -r colab/requirements.txt
if [ "${1:-}" = "acquisition" ]; then
  pip install -q -e ".[acquisition]"
else
  pip install -q -e .
fi
python -c "from insar_wetlands import paths; print('drive', paths.resolve_drive()); print('cube', paths.cube_root()); print('s1', paths.s1_local_root())"
echo ">>> ready"
