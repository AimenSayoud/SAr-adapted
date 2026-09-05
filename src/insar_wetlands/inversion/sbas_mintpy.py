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
# ISCE burst nomme les composantes connexes '_conncomp.tif' (sans underscore
# median) — indispensable aux corrections d'erreurs de deroulement.
mintpy.load.connCompFile     = {data_dir}/*/*_conncomp.tif
mintpy.load.demFile          = {data_dir}/{first_pair}/*_dem.tif
mintpy.load.incAngleFile     = {data_dir}/{first_pair}/*_lv_theta.tif

{reference_block}
mintpy.network.tempBaseMax   = {temp_base_max}
mintpy.network.excludeIfgIndex = {exclude_ifg}
# Selection de reseau alignee sur la Phase 3 : rejette les paires dont la
# coherence spatiale moyenne est sous le seuil (etes decorreles).
mintpy.network.coherenceBased = yes
mintpy.network.minCoherence   = {network_min_coherence}

# Correction des erreurs de deroulement : DESACTIVEE ('no'), et ce n'est
# plus une supposition mais un resultat mesure. 'bridging' echoue car aucun
# pixel n'appartient a la composante connexe de CHAQUE interferogramme
# (discussion MintPy #819). 'phase_closure', pourtant concu pour les
# reseaux redondants (346 paires ici), a lui-meme rapporte
# "number of common regions: 0" -> AUCUNE region n'est commune a tous les
# interferogrammes. C'est une mesure directe de la decorrelation du site
# (tourbiere + champs agricoles a 12j) : il n'y a rien a corriger avec ces
# methodes. Le SBAS reste donc epars mais non corrige des sauts de
# deroulement residuels -> role de comparaison/reference uniquement ; la
# mesure du fen est portee par l'ISBAS (Phase 9), qui n'a pas cette
# contrainte de composante connexe partagee.
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
                 ref_lat: float | None = None, ref_lon: float | None = None,
                 coh_threshold: float = 0.4,
                 temp_base_max: int = 48, tropo: bool = False,
                 exclude_ifg: str = "no",
                 network_min_coherence: float = 0.30,
                 unwrap_error_method: str = "no",
                 reference_min_coherence: float = 0.85,
                 ref_yx: tuple[int, int] | None = None) -> Path:
    """Ecrit la config MintPy.

    reference, par ordre de priorite :
      1. ref_yx=(row, col) : DETERMINISTE, recommande — utiliser
         best_reference_yx(work_dir) apres un premier run partiel ;
      2. ref_lat/ref_lon : impose un point geographique (echoue s'il est
         hors composante connexe) ;
      3. rien : selection auto MintPy (coherence > reference_min_coherence).
    On relit la reference effective avec read_reference() pour la partager
    avec l'ISBAS.

    unwrap_error_method : necessite les fichiers _conncomp.tif (verifier avec
    has_connected_component()). Passer 'no' pour desactiver.
    """
    # NB: chaque option doit etre EXPLICITE ('auto' compris) : MintPy
    # fusionne ce template dans son smallbaselineApp.cfg persistant — une
    # option omise garde sa vieille valeur d'un run precedent.
    if ref_yx is not None:
        # Reference deterministe (row,col) calculee par best_reference_yx()
        # -> garantie dans maskConnComp, pas de seuil arbitraire.
        reference_block = (
            f"mintpy.reference.yx           = {ref_yx[0]},{ref_yx[1]}\n"
            "mintpy.reference.lalo         = auto")
    elif ref_lat is not None and ref_lon is not None:
        reference_block = (
            f"mintpy.reference.lalo        = {ref_lat},{ref_lon}\n"
            "mintpy.reference.yx          = auto")
    else:
        # Selection auto MintPy : pixel avec coherence > minCoherence dans
        # maskConnComp. ATTENTION: echoue si aucun pixel ne depasse le seuil
        # (le defaut MintPy 0.85 est inatteignable sur tourbiere/champs).
        reference_block = (
            "mintpy.reference.lalo         = auto\n"
            "mintpy.reference.yx           = auto\n"
            f"mintpy.reference.minCoherence = {reference_min_coherence}")

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cfg = MINTPY_CFG_TEMPLATE.format(
        data_dir=str(data_dir), first_pair=first_pair,
        reference_block=reference_block,
        temp_base_max=temp_base_max, coh_threshold=coh_threshold,
        tropo_method="pyaps" if tropo else "no",
        exclude_ifg=exclude_ifg,
        network_min_coherence=network_min_coherence,
        unwrap_error_method=unwrap_error_method,
    )
    path = work_dir / "rzecin.cfg"
    path.write_text(cfg)
    return path


def reset_derived_products(work_dir: str | Path) -> list[str]:
    """Supprime les produits derives perimes du dossier de travail MintPy,
    en gardant inputs/ifgramStack.h5 (le stack deja charge depuis les 346
    GeoTIFF, couteux a reconstruire ~5 min) et rzecin.cfg.

    A utiliser quand un fichier intermediaire (maskConnComp.h5, etc.) est
    reste perime malgre la mise a jour des donnees source, a cause du cache
    incremental de MintPy ('output file already exists... skip').
    """
    work_dir = Path(work_dir)
    removed = []
    for pattern in ("*.h5", "*.png", "*.kmz", "*.txt", "*.cfg"):
        for f in work_dir.glob(pattern):
            # rzecin.cfg = notre template (source de verite, a garder) ;
            # smallbaselineApp.cfg = merge persistant genere par MintPy
            # (source du bug de valeurs perimees -> doit etre supprime).
            if f.name == "rzecin.cfg":
                continue
            f.unlink()
            removed.append(f.name)
    return removed


def inspect_h5(path: str | Path) -> None:
    """Diagnostic : structure REELLE d'un fichier HDF5 MintPy (cles, dtype,
    shape, min/max/nombre de valeurs distinctes). A utiliser au lieu de
    deviner le format d'un fichier (mask*, coherence...) avant de l'exploiter."""
    import h5py
    import numpy as np

    with h5py.File(path) as f:
        print(f"=== {path} ===")
        for k in f.keys():
            ds = f[k]
            if not hasattr(ds, "shape"):
                continue
            arr = ds[()]
            info = f"  {k}: dtype={arr.dtype} shape={arr.shape}"
            if arr.size and (np.issubdtype(arr.dtype, np.number)
                             or np.issubdtype(arr.dtype, np.bool_)):
                finite = arr[np.isfinite(arr)] if np.issubdtype(arr.dtype, np.floating) else arr
                uniq = np.unique(finite)
                info += f" min={finite.min():.4g} max={finite.max():.4g}"
                info += (f" n_unique={len(uniq)}"
                        + (f" values={uniq[:10].tolist()}" if len(uniq) <= 10 else ""))
            print(info)
        print("  attrs:", {k: v for k, v in list(f.attrs.items())[:10]})


def _first_dataset(h5file) -> str:
    """Nom du premier dataset 2D d'un fichier HDF5 MintPy."""
    for k in h5file.keys():
        if getattr(h5file[k], "ndim", 0) == 2:
            return k
    raise KeyError(f"aucun dataset 2D dans {h5file.filename}")


def best_reference_yx(work_dir: str | Path) -> dict:
    """Reference DETERMINISTE : pixel de coherence spatiale moyenne maximale
    a l'interieur de maskConnComp (donc garanti acceptable par MintPy).

    Se calcule apres l'etape load+quick_overview (avgSpatialCoh.h5 existe).
    Remplace la selection auto de MintPy dont le seuil par defaut (0.85) est
    inatteignable sur tourbiere/champs (~max 0.6-0.7 ici).
    """
    import h5py
    import numpy as np

    work_dir = Path(work_dir)
    with h5py.File(work_dir / "avgSpatialCoh.h5") as f:
        coh = f[_first_dataset(f)][:]
    mask_path = work_dir / "maskConnComp.h5"
    n_valid = coh.size
    if mask_path.exists():
        with h5py.File(mask_path) as f:
            raw = f[_first_dataset(f)][:]
        # Convention MintPy : 0 = hors composante connexe commune, tout
        # entier > 0 = ID de composante valide (jamais de valeurs negatives
        # en usage normal -> pas d'ambiguite avec un cast bool naif).
        m = raw > 0
        n_valid = int(m.sum())
        if n_valid == 0:
            raise RuntimeError(
                f"maskConnComp.h5 ({mask_path}) ne contient AUCUN pixel > 0 "
                f"(dtype={raw.dtype}, min={raw.min()}, max={raw.max()}) -> "
                "fichier probablement perime (genere avant l'ajout des "
                "_conncomp.tif) ou format inattendu. Diagnostiquer avec "
                "inspect_h5() avant de continuer ; eventuellement supprimer "
                "ce fichier + avgSpatialCoh.h5 pour forcer leur regeneration."
            )
        coh = np.where(m, coh, -np.inf)
    iy, ix = np.unravel_index(np.nanargmax(coh), coh.shape)
    best = float(coh[iy, ix])
    if not np.isfinite(best):
        raise RuntimeError("aucun pixel de coherence finie trouve dans le mask")
    return {"row": int(iy), "col": int(ix), "avg_spatial_coherence": best,
            "n_valid_pixels_in_mask": n_valid}


def read_reference(work_dir: str | Path) -> dict:
    """Relit le point de reference reellement utilise par MintPy.

    Lu depuis les attributs REF_* de timeseries.h5 -> a reutiliser tel quel
    pour l'ISBAS (Phase 9) afin que SBAS et ISBAS partagent le meme zero.
    """
    import h5py

    with h5py.File(Path(work_dir) / "timeseries.h5") as f:
        a = dict(f.attrs)

    def _get(k):
        v = a.get(k)
        return v.decode() if isinstance(v, bytes) else v

    x0 = float(a.get("X_FIRST", 0)); dx = float(a.get("X_STEP", 1))
    y0 = float(a.get("Y_FIRST", 0)); dy = float(a.get("Y_STEP", -1))
    ry, rx = int(_get("REF_Y")), int(_get("REF_X"))
    return {
        "y": y0 + dy * ry, "x": x0 + dx * rx,
        "row": ry, "col": rx,
        "lat": float(_get("REF_LAT")) if _get("REF_LAT") is not None else None,
        "lon": float(_get("REF_LON")) if _get("REF_LON") is not None else None,
    }


def has_connected_component(cropped_root: str | Path) -> bool:
    """Vrai si des fichiers de composantes connexes existent dans les
    produits croppes (necessaire aux corrections d'erreurs de deroulement).
    Les produits ISCE burst utilisent '_conncomp.tif', les GAMMA
    '_conn_comp.tif'."""
    root = Path(cropped_root)
    return (any(root.glob("*/*_conncomp.tif"))
            or any(root.glob("*/*_conn_comp.tif")))


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
