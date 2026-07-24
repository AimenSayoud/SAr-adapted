"""Lecture de la configuration et des secrets (Colab userdata ou variables d'env)."""

from __future__ import annotations

import os
from pathlib import Path

import yaml


def repo_root() -> Path:
    """Racine du repo = dossier contenant config/config.yaml, en remontant depuis cwd."""
    p = Path.cwd().resolve()
    for cand in [p, *p.parents]:
        if (cand / "config" / "config.yaml").exists():
            return cand
    raise FileNotFoundError(
        "config/config.yaml introuvable — lancer le notebook depuis le repo clone."
    )


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path else repo_root() / "config" / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def get_secret(name: str, required: bool = True) -> str | None:
    """Recupere un secret : Colab userdata d'abord, puis variable d'environnement."""
    value = None
    try:
        from google.colab import userdata  # type: ignore

        try:
            value = userdata.get(name)
        except Exception:
            value = None
    except ImportError:
        pass
    if not value:
        value = os.environ.get(name)
    if required and not value:
        raise RuntimeError(
            f"Secret '{name}' introuvable. L'ajouter dans Colab (icone cle) "
            f"et autoriser l'acces pour ce notebook."
        )
    return value


def setup_earthdata() -> None:
    """Ecrit ~/.netrc pour asf_search / hyp3_sdk (NASA Earthdata)."""
    user = get_secret("EARTHDATA_USERNAME")
    pwd = get_secret("EARTHDATA_PASSWORD")
    netrc = Path.home() / ".netrc"
    entry = f"machine urs.earthdata.nasa.gov login {user} password {pwd}\n"
    existing = netrc.read_text() if netrc.exists() else ""
    if "urs.earthdata.nasa.gov" not in existing:
        netrc.write_text(existing + entry)
    netrc.chmod(0o600)


def setup_cdsapi() -> None:
    """Ecrit ~/.cdsapirc pour l'API Climate Data Store (ERA5)."""
    key = get_secret("CDSAPI_KEY")
    rc = Path.home() / ".cdsapirc"
    rc.write_text(f"url: https://cds.climate.copernicus.eu/api\nkey: {key}\n")
    rc.chmod(0o600)


def setup_git(repo_dir: str | Path | None = None) -> None:
    """Configure git (identite + remote avec token) pour pousser les outputs depuis Colab."""
    import subprocess

    email = get_secret("GIT_USER_EMAIL")
    token = get_secret("GITHUB_TOKEN")
    cwd = str(repo_dir or repo_root())

    def run(*args: str) -> None:
        subprocess.run(["git", *args], cwd=cwd, check=True)

    run("config", "user.email", email)
    run("config", "user.name", email.split("@")[0])
    url = subprocess.run(
        ["git", "remote", "get-url", "origin"], cwd=cwd, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    if "@" not in url and url.startswith("https://"):
        run("remote", "set-url", "origin",
            url.replace("https://", f"https://x-access-token:{token}@"))
