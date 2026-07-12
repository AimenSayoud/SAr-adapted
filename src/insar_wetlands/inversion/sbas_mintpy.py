"""Phase 8 — Inversion SBAS classique via MintPy (smallbaselineApp).

MintPy lit directement les produits HyP3 croppes (processor = hyp3).
La contrainte 'indice de qualite' est appliquee en amont : on ne donne a
MintPy que le masque W >= seuil (mintpy.network + mask).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

MINTPY_CFG_TEMPLATE = """\
mintpy.load.processor        = hyp3
mintpy.load.unwFile          = {data_dir}/*/*_unw_phase.tif
mintpy.load.corFile          = {data_dir}/*/*_corr.tif
mintpy.load.connCompFile     = {data_dir}/*/*_conn_comp.tif
mintpy.load.demFile          = {data_dir}/{first_pair}/*_dem.tif
mintpy.load.incAngleFile     = {data_dir}/{first_pair}/*_lv_theta.tif
mintpy.load.waterMaskFile    = {data_dir}/{first_pair}/*_water_mask.tif

mintpy.reference.lalo        = {ref_lat},{ref_lon}
mintpy.network.tempBaseMax   = {temp_base_max}
mintpy.network.excludeIfgIndex = {exclude_ifg}

mintpy.networkInversion.weightFunc    = var
mintpy.networkInversion.maskDataset   = coherence
mintpy.networkInversion.maskThreshold = {coh_threshold}

mintpy.troposphericDelay.method = {tropo_method}
mintpy.troposphericDelay.weatherModel = ERA5
mintpy.deramp                = no
mintpy.topographicResidual   = yes
"""


def write_config(work_dir: str | Path, data_dir: str | Path, first_pair: str,
                 ref_lat: float, ref_lon: float, coh_threshold: float = 0.4,
                 temp_base_max: int = 48, tropo: bool = False,
                 exclude_ifg: str = "no") -> Path:
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cfg = MINTPY_CFG_TEMPLATE.format(
        data_dir=str(data_dir), first_pair=first_pair,
        ref_lat=ref_lat, ref_lon=ref_lon,
        temp_base_max=temp_base_max, coh_threshold=coh_threshold,
        tropo_method="pyaps" if tropo else "no",
        exclude_ifg=exclude_ifg,
    )
    path = work_dir / "rzecin.cfg"
    path.write_text(cfg)
    return path


def run(cfg_path: str | Path, work_dir: str | Path,
        dostep: str | None = None) -> int:
    """Lance smallbaselineApp.py ; retourne le code de sortie (log affiche)."""
    cmd = ["smallbaselineApp.py", str(cfg_path), "--work-dir", str(work_dir)]
    if dostep:
        cmd += ["--dostep", dostep]
    proc = subprocess.run(cmd, text=True)
    return proc.returncode


def load_timeseries(work_dir: str | Path):
    """Charge la serie temporelle MintPy (timeseries.h5) en xarray."""
    import h5py
    import numpy as np
    import pandas as pd
    import xarray as xr

    ts_file = Path(work_dir) / "timeseries.h5"
    with h5py.File(ts_file) as f:
        data = f["timeseries"][:]          # (n_dates, y, x), metres
        dates = pd.to_datetime([d.decode() for d in f["date"][:]])
        attrs = dict(f.attrs)
    x0 = float(attrs.get("X_FIRST", 0)); dx = float(attrs.get("X_STEP", 1))
    y0 = float(attrs.get("Y_FIRST", 0)); dy = float(attrs.get("Y_STEP", -1))
    ny, nx = data.shape[1:]
    return xr.DataArray(
        data * 1000.0,  # -> mm
        dims=("time", "y", "x"),
        coords={"time": dates,
                "y": y0 + dy * np.arange(ny),
                "x": x0 + dx * np.arange(nx)},
        name="los_displacement_mm", attrs={"units": "mm"},
    )
