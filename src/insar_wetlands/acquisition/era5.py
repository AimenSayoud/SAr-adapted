"""Telechargement ERA5 (precipitation, temperature, vapeur d'eau) via l'API CDS.

Le nouveau backend CDS impose des limites de cout par requete ("cost limits
exceeded" si on demande plusieurs annees d'un coup). Strategie :
  1 requete par ANNEE (fichiers annuels caches sur Drive, idempotent) ;
  si une annee est encore trop grosse -> repli : 1 requete par (annee, variable).
Les fichiers partiels sont ensuite fusionnes en un seul era5_rzecin.nc.
"""

from __future__ import annotations

from pathlib import Path


def _retrieve(cfg: dict, out: Path, year: int,
              variables: list[str]) -> None:
    import cdsapi

    e5 = cfg["era5"]
    c = cdsapi.Client()
    c.retrieve(
        e5["dataset"],
        {
            "product_type": "reanalysis",
            "variable": variables,
            "year": [str(year)],
            "month": [f"{m:02d}" for m in range(1, 13)],
            "day": [f"{d:02d}" for d in range(1, 32)],
            "time": [f"{h:02d}:00" for h in range(0, 24, 6)],
            "area": e5["area"],  # N, W, S, E
            "format": "netcdf",
        },
        str(out),
    )


def _download_year(cfg: dict, cache_dir: Path, year: int) -> Path:
    """Fichier annuel, avec repli par variable si limite de cout depassee."""
    import requests
    import xarray as xr

    out = cache_dir / f"era5_{year}.nc"
    if out.exists():
        return out
    variables = cfg["era5"]["variables"]
    try:
        _retrieve(cfg, out, year, variables)
        return out
    except requests.HTTPError as e:
        if "cost" not in str(e).lower():
            raise
        print(f"  {year}: requete annuelle trop grosse -> repli par variable")
    parts = []
    for var in variables:
        pv = cache_dir / f"era5_{year}_{var}.nc"
        if not pv.exists():
            _retrieve(cfg, pv, year, [var])
        parts.append(pv)
    merged = xr.merge([xr.open_dataset(p) for p in parts])
    merged.to_netcdf(out)
    for p in parts:
        p.unlink(missing_ok=True)
    return out


def download_era5(cfg: dict, out_path: str | Path) -> Path:
    """Telecharge ERA5 6-horaire par annee et fusionne. Idempotent.

    Chaque requete CDS peut rester en file d'attente plusieurs minutes.
    """
    import xarray as xr

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        return out_path

    cache_dir = out_path.parent / "era5_yearly"
    cache_dir.mkdir(parents=True, exist_ok=True)
    y0 = int(cfg["time"]["start"][:4])
    y1 = int(cfg["time"]["end"][:4])

    yearly = []
    for year in range(y0, y1 + 1):
        print(f"ERA5 {year}...")
        yearly.append(_download_year(cfg, cache_dir, year))

    datasets = [xr.open_dataset(p) for p in yearly]
    time_dim = "valid_time" if "valid_time" in datasets[0].dims else "time"
    merged = xr.concat(datasets, dim=time_dim).sortby(time_dim)
    tmp = out_path.with_suffix(".tmp.nc")
    merged.to_netcdf(tmp)
    for ds in datasets:
        ds.close()
    tmp.replace(out_path)
    return out_path
