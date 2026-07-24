"""NISAR L-band GUNW -> même format que nos crops HyP3 (unw_phase + corr).

Le produit **GUNW** de NISAR (Geocoded Unwrapped Interferogram, L-band, λ≈24 cm)
est l'exact analogue de nos interférogrammes HyP3 : phase déroulée géocodée +
**coherenceMagnitude** + interférogramme enroulé + connected components, sur une
grille géographique (posting 80 m). On peut donc rebrancher TOUT le pipeline
existant (define_zones, coherence_by_zone, evd_phase_linking) sans ISCE.

Pourquoi ce module et pas un simple `load_layer` : le GUNW est un **HDF5**
(pas un GeoTIFF), en statut **beta V1** — les chemins internes peuvent bouger.
On ne code donc AUCUN chemin en dur : `explore_h5` imprime l'arbre, et
`find_gunw_layers` localise phase/cohérence/coordonnées **par nom** (robuste au
versionnage). Ainsi le notebook reste valide quand l'archive/le format évolue.

Intérêt scientifique : la bande L (24 cm) **pénètre la canopée** et voit la
surface du tapis, là où le C-band (5.5 cm) décorrèle par diffusion de volume
(verdict Phase D/E2). Si la temporal_coherence L-band remonte sur la zone A
(tapis) — là où le C-band est au plancher de bruit — le facteur limitant est la
longueur d'onde, et un déplacement du tapis redevient mesurable.

NB disponibilité : données NISAR publiques depuis 2026-07-20, observations à
partir de ~juin 2026 (record complet fin 2026). Test **prospectif**, pas
rétroactif sur 2022-2024.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

_TS = re.compile(r"(\d{8})T\d{6}Z")


def search_gunw(aoi_wkt: str, start: str, end: str,
                processing_level: str = "GUNW", **kwargs) -> list:
    """Cherche les GUNW NISAR intersectant l'AOI (asf_search).

    `processing_level`/kwargs sont passés tels quels à asf_search pour rester
    robuste au nommage beta (ex. collections=['NISAR_L2_GUNW_BETA_V1']). Retourne
    la liste brute des résultats asf_search. Peut être vide si l'archive ne
    couvre pas encore la Pologne — c'est normal en 2026, re-exécuter plus tard.
    """
    import asf_search as asf

    params = dict(dataset=["NISAR"], processingLevel=[processing_level],
                  intersectsWith=aoi_wkt, start=start, end=end)
    params.update(kwargs)
    try:
        results = list(asf.search(**params))
    except Exception as e:  # nommage beta pas encore stable dans asf_search
        print(f"  ! asf.search(dataset='NISAR', ...) a échoué ({e}).")
        print("    Essayez collections=['NISAR_L2_GUNW_BETA_V1'] ou vérifiez le "
              "nom de la collection dans ASF Vertex / le NISAR Data User Guide.")
        return []
    print(f"  {len(results)} GUNW trouvés sur {start}..{end}")
    return results


def gunw_pair_id(name: str) -> str:
    """'YYYYMMDD_YYYYMMDD' (réf_sec) depuis un nom de granule GUNW.

    Le nom contient les deux horodatages ref/sec (…_YYYYMMDDThhmmssZ_…). On prend
    les deux premiers et on ordonne (ancien_récent), format identique au reste du
    pipeline (dates_from_pairs, _pair_date_index)."""
    ts = _TS.findall(name)
    if len(ts) < 2:
        raise ValueError(f"impossible d'extraire ref/sec de: {name}")
    a, b = sorted(ts[:2])
    return f"{a}_{b}"


def explore_h5(path: str | Path, max_depth: int = 6) -> None:
    """Imprime l'arbre des datasets d'un GUNW (diagnostic runtime).

    À lancer UNE fois sur un fichier téléchargé pour voir les chemins réels
    (le format beta peut différer de la doc) avant de charger."""
    import h5py

    with h5py.File(path, "r") as f:
        def _walk(name, obj):
            depth = name.count("/")
            if depth <= max_depth and isinstance(obj, h5py.Dataset):
                print(f"  {name}  shape={obj.shape} dtype={obj.dtype}")
        f.visititems(_walk)


def find_gunw_layers(path: str | Path, freq: str = "frequencyA",
                     pol: str | None = None) -> dict:
    """Localise PAR NOM les datasets clés d'un GUNW (robuste au versionnage).

    Retourne un dict de chemins HDF5 : {'unw', 'coh', 'x', 'y', 'epsg', 'pol'}.
    Heuristique : on privilégie le groupe `freq` (frequencyA) et la polarisation
    `pol` si fournie, sinon la 1re polarisation disponible. Les noms cherchés :
    unwrapped phase, coherence magnitude, xCoordinates/yCoordinates, projection.
    """
    import h5py

    ds_paths = []
    with h5py.File(path, "r") as f:
        f.visititems(lambda n, o: ds_paths.append(n)
                     if isinstance(o, h5py.Dataset) else None)

        def pick(pred):
            cands = [p for p in ds_paths if pred(p.lower())]
            if freq.lower() in " ".join(cands).lower():
                cands = [p for p in cands if freq.lower() in p.lower()] or cands
            if pol:
                pol_c = [p for p in cands if f"/{pol.lower()}" in p.lower()
                         or p.lower().endswith(pol.lower())]
                cands = pol_c or cands
            return cands[0] if cands else None

        unw = pick(lambda p: "unwrapp" in p and "phase" in p)
        coh = pick(lambda p: "coherence" in p)
        x = pick(lambda p: p.endswith("xcoordinates") or p.endswith("/xcoordinate"))
        y = pick(lambda p: p.endswith("ycoordinates") or p.endswith("/ycoordinate"))
        proj = pick(lambda p: "projection" in p or p.endswith("epsg"))
        epsg = None
        if proj is not None:
            v = f[proj][()]
            epsg = int(np.ravel(v)[0]) if np.ndim(v) else int(v)
        # polarisation déduite du chemin unw (…/HH/…)
        pol_found = None
        if unw:
            m = re.search(r"/(HH|HV|VH|VV|RH|RV)/", unw + "/", re.I)
            pol_found = m.group(1).upper() if m else pol
    layers = {"unw": unw, "coh": coh, "x": x, "y": y, "epsg": epsg,
              "pol": pol_found}
    missing = [k for k in ("unw", "coh", "x", "y") if not layers[k]]
    if missing:
        raise KeyError(f"couches introuvables {missing} — lancez explore_h5() "
                       f"pour voir l'arbre réel. Trouvé: {layers}")
    return layers


def load_gunw(path: str | Path, freq: str = "frequencyA",
              pol: str | None = None) -> xr.Dataset:
    """Charge un GUNW en Dataset rioxarray {unw_phase, corr} géoréférencé.

    Mêmes noms de variables que nos crops HyP3 (unw_phase, corr) → compatible
    avec load_layer/define_zones/evd_phase_linking en aval."""
    import h5py

    lyr = find_gunw_layers(path, freq, pol)
    with h5py.File(path, "r") as f:
        unw = np.asarray(f[lyr["unw"]][()], "float32")
        coh = np.asarray(f[lyr["coh"]][()], "float32")
        x = np.asarray(f[lyr["x"]][()], "float64")
        y = np.asarray(f[lyr["y"]][()], "float64")
    # cohérence et phase peuvent avoir des postings différents (20 vs 80 m) :
    # on aligne la cohérence sur la grille de la phase si les tailles diffèrent.
    ds = xr.Dataset(
        {"unw_phase": (("y", "x"), unw)},
        coords={"y": y, "x": x},
    )
    if coh.shape != unw.shape:
        coh_da = xr.DataArray(coh, dims=("y", "x"),
                              coords={"y": np.linspace(y[0], y[-1], coh.shape[0]),
                                      "x": np.linspace(x[0], x[-1], coh.shape[1])})
        coh = coh_da.interp(y=ds.y, x=ds.x, method="nearest").values
    ds["corr"] = (("y", "x"), coh)
    ds = ds.rio.write_crs(f"EPSG:{lyr['epsg']}") if lyr["epsg"] else ds
    ds.attrs["polarization"] = lyr["pol"]
    return ds


def build_gunw_stack(files: list[str | Path], template: xr.DataArray | None = None,
                     freq: str = "frequencyA", pol: str | None = None
                     ) -> xr.Dataset:
    """Empile plusieurs GUNW en (pair, y, x) : {unw_phase, corr}.

    Si `template` est fourni (ex. la grille d'un crop C-band), chaque GUNW y est
    reprojeté (reproject_match) pour une comparaison pixel-à-pixel C vs L. Sinon
    on aligne tous les GUNW sur la grille du premier.
    """
    layers = []
    pairs = []
    ref = None
    for fp in sorted(files, key=lambda p: gunw_pair_id(Path(p).name)):
        ds = load_gunw(fp, freq, pol)
        if template is not None:
            ds = ds.rio.reproject_match(template)
        elif ref is None:
            ref = ds
        elif ds["unw_phase"].shape != ref["unw_phase"].shape:
            ds = ds.rio.reproject_match(ref["unw_phase"])
        pairs.append(gunw_pair_id(Path(fp).name))
        layers.append(ds)
    stack = xr.concat(layers, dim="pair")
    stack = stack.assign_coords(pair=pairs)
    return stack
