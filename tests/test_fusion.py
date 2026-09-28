"""Ground truth for the fusion blocks (X-062): each routine recovers a known answer."""
import numpy as np
import pytest

from insar_wetlands import fusion as fu


def test_hamon_pet_is_small_in_frost_and_largest_in_summer():
    doy = np.array([15, 196])
    pet = fu.hamon_pet(np.array([-10.0, 20.0]), doy, 52.76)
    assert pet[0] < 0.5 and 2.5 < pet[1] < 6.0                  # Hamon never reaches exactly 0
    assert fu.day_length_h(172, 52.76) == pytest.approx(16.6, abs=0.4)


def test_bucket_recovers_its_parameters_and_forecasts():
    rng = np.random.default_rng(1)
    n = 1500
    rain = rng.gamma(0.4, 5.0, n) * (rng.random(n) < 0.4)
    pet = 2 + 1.5 * np.sin(2 * np.pi * np.arange(n) / 365)
    truth = {"a": 0.25, "b": 0.35, "c": 0.04, "w0": -5.0}
    w = fu.simulate_bucket(truth, -8.0, rain, pet) + rng.normal(0, 0.05, n)
    fit = fu.fit_bucket(w, rain, pet)
    for k in ("a", "b", "c"):
        assert fit[k] == pytest.approx(truth[k], rel=0.15), k
    sim = fu.simulate_bucket(fit, w[1000], rain[1000:], pet[1000:])
    assert fu.nse(w[1000:], sim) > 0.8


def test_kmeans_finds_the_planted_groups_and_the_index_flags_the_unrepresented():
    rng = np.random.default_rng(2)
    centres = np.array([[0, 0], [6, 0], [0, 6]], float)             # informative features only: an
    X = np.vstack([c + rng.normal(0, 0.5, (300, 2)) for c in centres])   # uninformative one dilutes the index
    Z, _, _ = fu.standardize(X)
    labels, k, _ = fu.choose_kmeans(Z, ks=range(2, 6))
    assert k == 3
    for g in range(3):
        assert len(np.unique(labels[g * 300:(g + 1) * 300])) == 1
    train = np.array([0, 1, 2, 300, 301, 302])                    # plots only in the first two groups
    di, thr = fu.dissimilarity_index(Z, train)
    outside = di > thr
    assert outside[600:].mean() > 0.95                             # the third group: nothing measured there
    assert outside[:600].mean() < 0.2


def test_gp_transfer_skill_is_positive_for_a_smooth_trait_and_not_for_noise():
    rng = np.random.default_rng(3)
    X = rng.uniform(0, 1, (12, 2))
    Xall = rng.uniform(0, 1, (200, 2))
    y = 2 * X[:, 0] + 0.02 * rng.normal(size=12)
    good = fu.gp_transfer(X, y, Xall)
    assert good["skill"] > 0.7
    assert np.corrcoef(good["mean"], 2 * Xall[:, 0])[0, 1] > 0.9
    noise = fu.gp_transfer(X, rng.normal(size=12), Xall)
    assert noise["skill"] < 0.2


def test_phase_noise_falls_with_coherence():
    s = fu.phase_sigma_mm([0.2, 0.5, 0.9])
    assert s[0] > s[1] > s[2] > 0


def test_bayes_linear_separates_motion_from_moisture():
    rng = np.random.default_rng(4)
    n = 60
    motion = rng.normal(0, 8, n)                                   # laser LOS change (mm)
    dwtd = 0.6 * motion / 3 + rng.normal(0, 2, n)                  # water change, correlated with motion
    sigma = fu.phase_sigma_mm(rng.uniform(0.3, 0.8, n))
    y = 1.0 * motion - 0.8 * dwtd + rng.normal(0, 1, n) * sigma
    fit = fu.bayes_linear(np.column_stack([motion, dwtd, np.ones(n)]), y, sigma)
    assert fit["beta"][0] == pytest.approx(1.0, abs=0.1)
    assert fit["beta"][1] == pytest.approx(-0.8, abs=0.3)
    assert fit["noise_scale"] == pytest.approx(1.0, abs=0.35)


def test_fusion_beats_both_inputs_and_rejects_an_unwrapping_error():
    rng = np.random.default_rng(5)
    n = 400
    truth = rng.normal(0, 8, n)
    Q, R = np.full(n, 25.0), np.full(n, 16.0)
    u = truth + rng.normal(0, 5, n)                                # model: σ 5 mm
    z = truth + rng.normal(0, 4, n)                                # radar: σ 4 mm
    z[10] += 27.7                                                  # one 2π unwrapping error
    out = fu.fuse_increments(u, Q, z, R)
    ok = np.arange(n) != 10
    rm = lambda a: np.sqrt(np.mean((a[ok] - truth[ok]) ** 2))      # noqa: E731
    assert out["rejected"][10]
    assert rm(out["inc"]) < min(rm(u), rm(z)) * 0.95
    assert out["cum"].shape[0] == n + 1 and out["cum"][0] == 0


def test_harmonic_amplitude_map_reads_each_pixel():
    t = np.linspace(0, 3, 90)
    Y = np.column_stack([3 * np.cos(2 * np.pi * t), 7 * np.sin(2 * np.pi * t) + t])
    assert fu.harmonic_amplitude_map(t, Y) == pytest.approx([3, 7], abs=0.05)
