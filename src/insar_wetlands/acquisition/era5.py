"""Telechargement ERA5 (precipitation, temperature, vapeur d'eau) via l'API CDS."""

from __future__ import annotations

from pathlib import Path


def download_era5(cfg: dict, out_path: str | Path) -> Path:
    """Telecharge ERA5 single-levels horaire sur la petite bbox du config.

    Fichier resultant : quelques dizaines de Mo pour 5 ans sur 1x1 deg.
    La requete CDS peut rester en file d'attente plusieurs minutes.
    """
    import cdsapi

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        return out_path

    e5 = cfg["era5"]
    y0 = int(cfg["time"]["start"][:4])
    y1 = int(cfg["time"]["end"][:4])

    c = cdsapi.Client()
    c.retrieve(
        e5["dataset"],
        {
            "product_type": "reanalysis",
            "variable": e5["variables"],
            "year": [str(y) for y in range(y0, y1 + 1)],
            "month": [f"{m:02d}" for m in range(1, 13)],
            "day": [f"{d:02d}" for d in range(1, 32)],
            "time": [f"{h:02d}:00" for h in range(0, 24, 6)],
            "area": e5["area"],  # N, W, S, E
            "format": "netcdf",
        },
        str(out_path),
    )
    return out_path
