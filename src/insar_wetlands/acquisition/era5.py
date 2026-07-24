"""Telechargement ERA5 (precipitation, temperature, vapeur d'eau) via l'API CDS.

Deux pieges du backend CDS geres ici :
  1. Limite de cout par requete ("cost limits exceeded" si plusieurs annees
     d'un coup) -> 1 requete par ANNEE, avec repli par variable si besoin.
  2. Conversion GRIB->netCDF cote serveur (cfgrib) qui separe les variables
     "instantanees" (t2m, tcwv) des variables "accumulees" (tp) : quand elles
     sont demandees ensemble, le fichier resultant ne garde parfois que le
     groupe instantane et perd silencieusement 'total_precipitation', SANS
     erreur HTTP. On verifie donc apres coup les variables presentes et on
     rattrape celles manquantes par une requete separee.
Tout est idempotent (fichiers annuels caches sur Drive) et relancable.
"""

from __future__ import annotations

from pathlib import Path

# Mapping nom CDS long -> nom court dans le netCDF resultant (cfgrib).
VAR_SHORTNAME = {
    "total_precipitation": "tp",
    "2m_temperature": "t2m",
    "total_column_water_vapour": "tcwv",
}


def _ensure_cdsapirc() -> None:
    """Reecrit ~/.cdsapirc si absent (ex: kernel Colab redemarre entre-temps)."""
    from pathlib import Path as _Path

    if not (_Path.home() / ".cdsapirc").exists():
        from ..config import setup_cdsapi

        setup_cdsapi()


def _retrieve(cfg: dict, out: Path, year: int, variables: list[str]) -> None:
    import cdsapi

    _ensure_cdsapirc()
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
            # Le nouveau backend CDS zippe par defaut meme avec un nom de
            # sortie en .nc -> forcer un fichier netCDF brut, non archive.
            "data_format": "netcdf",
            "download_format": "unarchived",
        },
        str(out),
    )
    _unwrap_if_zipped(out)


def _unwrap_if_zipped(path: Path) -> None:
    """Le CDS renvoie parfois un zip malgre 'unarchived' (bug connu) ->
    extrait le .nc qu'il contient et remplace le fichier."""
    import zipfile

    if not zipfile.is_zipfile(path):
        return
    with zipfile.ZipFile(path) as zf:
        members = [n for n in zf.namelist() if n.endswith(".nc")]
        if not members:
            raise ValueError(f"zip ERA5 sans fichier .nc: {path} ({zf.namelist()})")
        data = zf.read(members[0])
    path.write_bytes(data)


def _missing_variables(path: Path, variables: list[str]) -> list[str]:
    """Variables demandees (noms longs CDS) absentes du netCDF (nom court)."""
    import xarray as xr

    with xr.open_dataset(path) as ds:
        present = set(ds.data_vars)
    return [v for v in variables
            if VAR_SHORTNAME.get(v, v) not in present]


def _backfill_missing_vars(cfg: dict, out: Path, year: int,
                           missing: list[str]) -> None:
    """Retelecharge separement chaque variable manquante (stepType different,
    ex: precipitation accumulee vs temperature instantanee) et fusionne."""
    import xarray as xr

    print(f"  {year}: variable(s) manquante(s) apres conversion CDS -> "
          f"{missing} (repli par variable)")
    with xr.open_dataset(out) as base:
        base = base.load()
    for var in missing:
        part = out.parent / f"era5_{year}_{var}.nc"
        if not part.exists() or _missing_variables(part, [var]):
            _retrieve(cfg, part, year, [var])
        with xr.open_dataset(part) as ds_var:
            base = xr.merge([base, ds_var], join="outer", compat="override")
    tmp = out.with_suffix(".backfill.nc")
    base.to_netcdf(tmp)
    base.close()
    tmp.replace(out)


def _download_year(cfg: dict, cache_dir: Path, year: int) -> Path:
    """Fichier annuel complet (toutes variables), cache et verifie. Idempotent."""
    import requests

    out = cache_dir / f"era5_{year}.nc"
    variables = cfg["era5"]["variables"]

    if not out.exists():
        try:
            _retrieve(cfg, out, year, variables)
        except requests.HTTPError as e:
            if "cost" not in str(e).lower():
                raise
            print(f"  {year}: requete annuelle trop grosse -> repli par variable")
            _retrieve(cfg, out, year, [variables[0]])
            _backfill_missing_vars(cfg, out, year, variables[1:])
            return out
    else:
        _unwrap_if_zipped(out)

    missing = _missing_variables(out, variables)
    if missing:
        _backfill_missing_vars(cfg, out, year, missing)
    return out


def download_era5(cfg: dict, out_path: str | Path) -> Path:
    """Telecharge ERA5 6-horaire par annee et fusionne. Idempotent et
    auto-reparant (variables manquantes rattrapees a chaque appel).

    Chaque requete CDS peut rester en file d'attente plusieurs minutes.
    """
    import xarray as xr

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    variables = cfg["era5"]["variables"]

    if out_path.exists() and not _missing_variables(out_path, variables):
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
