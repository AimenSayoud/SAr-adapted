"""Test du lecteur NISAR GUNW (parsing + auto-découverte HDF5 + stack).

On fabrique un HDF5 minimal aux chemins-type GUNW et on vérifie que :
(1) gunw_pair_id extrait 'réf_sec' du nom de granule ;
(2) find_gunw_layers localise phase/cohérence/coordonnées PAR NOM ;
(3) load_gunw rend {unw_phase, corr} sur la bonne grille ;
(4) build_gunw_stack empile en (pair, y, x) au format du pipeline.

Round-trip h5py réel (pas de mock). Aucune dépendance rioxarray : on omet la
couche projection -> load_gunw saute la branche .rio (testé sans Colab).

Execution : python tests/test_nisar_gunw.py
"""

import numpy as np

from insar_wetlands.nisar import build_gunw_stack, find_gunw_layers, gunw_pair_id, load_gunw

BASE = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram"


def _fake_gunw(path, ny=5, nx=7, seed=0, with_cc=True):
    import h5py

    rng = np.random.default_rng(seed)
    with h5py.File(path, "w") as f:
        g = f.create_group(f"{BASE}/HH")
        g.create_dataset("unwrappedPhase", data=rng.normal(0, 1, (ny, nx)).astype("f4"))
        g.create_dataset("coherenceMagnitude", data=rng.uniform(0, 1, (ny, nx)).astype("f4"))
        if with_cc:
            cc = np.ones((ny, nx), "i2")
            cc[0, 0] = 0                     # un pixel non fiable (îlot déroulement)
            g.create_dataset("connectedComponents", data=cc)
        f.create_dataset(f"{BASE}/xCoordinates", data=np.arange(nx) * 80.0 + 5e5)
        f.create_dataset(f"{BASE}/yCoordinates", data=np.arange(ny) * -80.0 + 6e6)


def test_pair_id():
    n = "NISAR_L2_PR_GUNW_001_A_20260617T161200Z_20260629T161200Z_v1.h5"
    assert gunw_pair_id(n) == "20260617_20260629"
    # ordre inversé dans le nom -> toujours ancien_récent
    n2 = "X_20260629T000000Z_20260617T000000Z_.h5"
    assert gunw_pair_id(n2) == "20260617_20260629"


def test_find_and_load(tmp_path=None):
    import os
    import tempfile

    d = tempfile.mkdtemp()
    fp = os.path.join(d, "g.h5")
    _fake_gunw(fp)
    lyr = find_gunw_layers(fp)
    assert lyr["unw"].endswith("unwrappedPhase")
    assert "coherence" in lyr["coh"].lower()
    assert lyr["cc"].endswith("connectedComponents")
    assert lyr["pol"] == "HH"
    assert lyr["epsg"] is None            # pas de projection -> branche rio sautée
    ds = load_gunw(fp)
    assert set(ds.data_vars) == {"unw_phase", "corr"}
    assert ds["unw_phase"].shape == (5, 7)
    assert float(ds["corr"].min()) >= 0 and float(ds["corr"].max()) <= 1
    # cc==0 -> pixel (0,0) masqué en NaN ; les autres restent finis
    assert np.isnan(ds["unw_phase"].values[0, 0])
    assert np.isfinite(ds["unw_phase"].values[1, 1])


def test_build_stack():
    import os
    import tempfile

    d = tempfile.mkdtemp()
    f1 = os.path.join(d, "A_20260617T000000Z_20260629T000000Z.h5")
    f2 = os.path.join(d, "B_20260629T000000Z_20260711T000000Z.h5")
    _fake_gunw(f1, seed=1)
    _fake_gunw(f2, seed=2)
    stack = build_gunw_stack([f2, f1])          # ordre volontairement inversé
    assert list(stack.pair.values) == ["20260617_20260629", "20260629_20260711"]
    assert stack["unw_phase"].dims == ("pair", "y", "x")
    assert stack["unw_phase"].shape == (2, 5, 7)


if __name__ == "__main__":
    test_pair_id()
    test_find_and_load()
    test_build_stack()
    print("ALL NISAR GUNW TESTS PASSED")
