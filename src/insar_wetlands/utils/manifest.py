"""Manifeste standard ecrit par chaque phase dans outputs/phaseXX/."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


def write_manifest(phase: str, outdir: str | Path, params: dict,
                   products: dict) -> Path:
    """Ecrit outputs/<phase>/manifest.json (parametres + produits + horodatage)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "phase": phase,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "parameters": params,
        "products": products,
    }
    path = outdir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, default=str))
    return path
