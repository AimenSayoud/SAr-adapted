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
mintpy.load.demFile          = {data_dir}/{first_pair}/*_dem.tif
mintpy.load.incAngleFile     = {data_dir}/{first_pair}/*_lv_theta.tif
# NB: les produits HyP3 INSAR_ISCE_BURST ne fournissent ni conn_comp ni
# water_mask -> pas de connCompFile/waterMaskFile ici. Le masque d'eau vient
# de notre Phase 5 (dynamique), et la correction de deroulement utilise
# phase_closure (sans conn_comp).

mintpy.reference.lalo        = {ref_lat},{ref_lon}
mintpy.network.tempBaseMax   = {temp_base_max}
mintpy.network.excludeIfgIndex = {exclude_ifg}
# Selection de reseau alignee sur la Phase 3 : rejette les paires dont la
# coherence spatiale moyenne est sous le seuil (etes decorreles).
mintpy.network.coherenceBased = yes
mintpy.network.minCoherence   = {network_min_coherence}

# Correction des erreurs de deroulement : DESACTIVEE ('no') car les DEUX
# methodes MintPy (bridging ET phase_closure) exigent le dataset
# connectComponent, absent des produits HyP3 INSAR_ISCE_BURST. C'est une
# limitation documentee du SBAS standard sur ce type de produit ; la
# correction des sauts est faite en Phase 9 (ISBAS) par notre propre
# fermeture de triplets, independante de MintPy et sans conn_comp.
mintpy.unwrapError.method    = {unwrap_error_method}

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
                 exclude_ifg: str = "no",
                 network_min_coherence: float = 0.30,
                 unwrap_error_method: str = "no") -> Path:
    """Ecrit la config MintPy.

    unwrap_error_method : 'no' (defaut) car les methodes MintPy 'bridging' et
    'phase_closure' exigent toutes deux le dataset connectComponent, absent
    des produits HyP3 INSAR_ISCE_BURST. La correction des sauts de phase est
    faite en Phase 9 (ISBAS). Ne passer 'bridging'/'phase_closure' que si
    has_connected_component() est vrai (produits GAMMA, p.ex.).
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cfg = MINTPY_CFG_TEMPLATE.format(
        data_dir=str(data_dir), first_pair=first_pair,
        ref_lat=ref_lat, ref_lon=ref_lon,
        temp_base_max=temp_base_max, coh_threshold=coh_threshold,
        tropo_method="pyaps" if tropo else "no",
        exclude_ifg=exclude_ifg,
        network_min_coherence=network_min_coherence,
        unwrap_error_method=unwrap_error_method,
    )
    path = work_dir / "rzecin.cfg"
    path.write_text(cfg)
    return path


def has_connected_component(cropped_root: str | Path) -> bool:
    """Vrai si des fichiers conn_comp existent dans les produits croppes
    (condition necessaire pour la methode 'bridging')."""
    root = Path(cropped_root)
    return any(root.glob("*/*_conn_comp.tif"))


def run(cfg_path: str | Path, work_dir: str | Path,
        dostep: str | None = None, start: str | None = None,
        tail: int = 40) -> int:
    """Lance smallbaselineApp.py, capture le log et affiche sa fin.

    - Le log COMPLET (stdout+stderr) est ecrit dans work_dir/mintpy_run.log
      (persistant sur Drive) ;
    - les `tail` dernieres lignes sont affichees dans la cellule -> c'est la
      qu'apparait le Traceback/Error en cas d'echec ;
    - retourne le code de sortie (0 = succes).

    start='modify_network' relance apres un changement de config sans
    recharger les 346 GeoTIFF (le stack ifgramStack.h5 est reutilise).
    """
    cmd = ["smallbaselineApp.py", str(cfg_path), "--work-dir", str(work_dir)]
    if dostep:
        cmd += ["--dostep", dostep]
    if start:
        cmd += ["--start", start]

    proc = subprocess.run(cmd, text=True, capture_output=True)
    log = (proc.stdout or "") + "\n===== STDERR =====\n" + (proc.stderr or "")
    log_path = Path(work_dir) / "mintpy_run.log"
    log_path.write_text(log)

    lines = log.splitlines()
    print(f"EXIT CODE: {proc.returncode}  (log complet: {log_path})")
    print(f"----- {tail} dernieres lignes du log -----")
    print("\n".join(lines[-tail:]))
    return proc.returncode


def _grid_coords(attrs, ny, nx):
    import numpy as np

    x0 = float(attrs.get("X_FIRST", 0)); dx = float(attrs.get("X_STEP", 1))
    y0 = float(attrs.get("Y_FIRST", 0)); dy = float(attrs.get("Y_STEP", -1))
    return {"y": y0 + dy * np.arange(ny), "x": x0 + dx * np.arange(nx)}


def load_temporal_coherence(work_dir: str | Path):
    """Charge temporalCoherence.h5 (qualite d'inversion par pixel, 0-1).

    C'est le critere standard MintPy pour separer les pixels fiables du bruit
    de decorrelation : sans ce masque, la carte de vitesse montre tous les
    pixels, y compris ceux dont la phase est aleatoire.
    """
    import h5py
    import xarray as xr

    with h5py.File(Path(work_dir) / "temporalCoherence.h5") as f:
        data = f["temporalCoherence"][:]
        attrs = dict(f.attrs)
    return xr.DataArray(data, dims=("y", "x"),
                        coords=_grid_coords(attrs, *data.shape),
                        name="temporal_coherence")


def load_timeseries(work_dir: str | Path, coh_threshold: float | None = 0.7):
    """Charge la serie temporelle MintPy (timeseries.h5) en xarray (mm).

    Si coh_threshold est fourni, les pixels dont la coherence temporelle
    d'inversion est sous le seuil sont mis a NaN (recommande : 0.7 standard,
    0.5-0.6 acceptable sur tourbiere). Passer None pour la serie brute.
    """
    import h5py
    import numpy as np
    import pandas as pd
    import xarray as xr

    ts_file = Path(work_dir) / "timeseries.h5"
    with h5py.File(ts_file) as f:
        data = f["timeseries"][:]          # (n_dates, y, x), metres
        dates = pd.to_datetime([d.decode() for d in f["date"][:]])
        attrs = dict(f.attrs)
    ny, nx = data.shape[1:]
    ts = xr.DataArray(
        data * 1000.0,  # -> mm
        dims=("time", "y", "x"),
        coords={"time": dates, **_grid_coords(attrs, ny, nx)},
        name="los_displacement_mm", attrs={"units": "mm"},
    )
    if coh_threshold is not None:
        tcoh = load_temporal_coherence(work_dir)
        ts = ts.where(tcoh >= coh_threshold)
    return ts
