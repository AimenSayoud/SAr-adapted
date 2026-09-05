#!/usr/bin/env bash
# =========================================================
# Installation de l'environnement dans une session Colab.
# A executer au debut de CHAQUE session (le disque Colab
# est efface entre les sessions) :
#   !bash environment/colab_setup.sh
# =========================================================
set -e

# Reproducibility: use the pinned lock file when it exists, so a rerun installs
# the stack that produced the published numbers. Fall back to the loose list.
if [ -f environment/requirements-lock.txt ]; then
    echo ">>> Installation depuis requirements-lock.txt (versions epinglees)..."
    pip install -q -r environment/requirements-lock.txt
    pip install -q -e .
    return 0 2>/dev/null || exit 0
fi

echo ">>> Pas de lock file — installation des dependances pip (non epinglees)..."
pip install -q --upgrade \
    asf_search \
    hyp3_sdk \
    pystac-client \
    cdsapi \
    rioxarray \
    rasterio \
    geopandas \
    shapely \
    pyproj \
    xarray \
    h5py \
    netcdf4 \
    pyyaml \
    planetary-computer

# MintPy n'est necessaire qu'a partir de la Phase 8 — installation separee
# (plus lourde) pour ne pas ralentir les phases 1-7 :
#   pip install -q mintpy

echo ">>> Installation du package local insar_wetlands..."
pip install -q -e .

echo ">>> Environnement pret."
