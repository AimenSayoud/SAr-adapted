#!/usr/bin/env bash
# =========================================================
# Installation de l'environnement dans une session Colab.
# A executer au debut de CHAQUE session (le disque Colab
# est efface entre les sessions) :
#   !bash environment/colab_setup.sh
# =========================================================
set -e

# Configure git identity for commits from Colab (C-013)
git config user.email "aimen.sayoud.polska@gmail.com"
git config user.name "Aymen Sayoud"

echo ">>> Installation des dependances InSAR..."
pip install -q --upgrade -r environment/requirements.txt

echo ">>> Installation du package local insar_wetlands..."
pip install -q -e .

echo ">>> Environnement pret."
