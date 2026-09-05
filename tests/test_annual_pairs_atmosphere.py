"""Ground-truth tests for the annual-pair selection and the tropospheric model.

`annual_pairs` is the largest untested module (418 lines) and it decides which
interferograms get made — a bad selection is paid for in HyP3 credits and in a
noisier result, with nothing to say why. `atmosphere` converts water vapour into
millimetres of apparent displacement, so an error there moves a number in the
paper directly.

Each test constructs a case whose answer is known by construction."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from insar_wetlands.annual_pairs import (
    build_annual_pair_list, deramp, gather_candidate_dates, near_aoi_ring,
    pick_seasonal_dates, select_optimal_annual_pairs, select_topk_annual_pairs,
    surface_valid_mask,
)
from insar_wetlands.atmosphere import pair_delays_mm, zenith_wet_delay_mm


def inventory(dates) -> pd.DataFrame:
    return pd.DataFrame({"date": pd.to_datetime(dates),
                         "granule": [f"S1A_{d}" for d in dates]})


def scored(rows) -> pd.DataFrame:
    """rows: (year, 'YYYY-MM-DD', tcwv, atmo_score)"""
    return pd.DataFrame(
        [{"year": y, "date": pd.Timestamp(d), "granule": "g",
          "tcwv": t, "atmo_score": a} for y, d, t, a in rows])


# =========================================================================
# date selection
# =========================================================================

def test_the_date_nearest_the_target_is_picked_for_each_year():
    inv = inventory(["2022-04-02", "2022-04-14", "2022-05-20",
                     "2023-04-18", "2023-06-01"])
    got = pick_seasonal_dates(inv, target_month=4, target_day=15)
    assert list(got["year"]) == [2022, 2023]
    assert got.loc[0, "date"] == pd.Timestamp("2022-04-14")   # 1 day away
    assert got.loc[0, "gap_days"] == 1
    assert got.loc[1, "date"] == pd.Timestamp("2023-04-18")   # 3 days away


def test_a_year_with_only_distant_dates_still_yields_one_and_says_how_far():
    """Silently dropping the year would leave a gap in the series with no
    explanation; the gap is reported instead."""
    inv = inventory(["2022-04-15", "2023-11-30"])
    got = pick_seasonal_dates(inv, target_month=4, target_day=15)
    assert len(got) == 2
    assert got.loc[1, "gap_days"] > 200


def test_the_candidate_pool_keeps_everything_inside_the_window():
    """This is the difference from `pick_seasonal_dates`: a pool to optimise
    over, not one date taken on trust."""
    inv = inventory(["2022-04-01", "2022-04-14", "2022-04-27", "2022-05-30"])
    pool = gather_candidate_dates(inv, target_month=4, target_day=15,
                                  window_days=15)
    assert len(pool) == 3                       # 30 May is 45 days out
    assert pool["gap_days"].abs().max() <= 15


def test_the_window_is_symmetric_about_the_target():
    inv = inventory(["2022-03-31", "2022-04-30"])      # -15 and +15
    pool = gather_candidate_dates(inv, target_month=4, target_day=15,
                                  window_days=15)
    assert len(pool) == 2
    assert sorted(pool["gap_days"]) == [-15, 15]


# =========================================================================
# pair optimisation
# =========================================================================

def test_the_pair_minimising_differential_water_vapour_is_chosen():
    """Only the DIFFERENCE in vapour survives interferogram formation, so the
    selection must minimise |Δtcwv| — not pick the individually driest date."""
    cand = scored([
        (2022, "2022-04-10", 20.0, 0.0),
        (2022, "2022-04-14", 12.0, 0.0),     # driest in 2022
        (2023, "2023-04-12", 20.0, 0.0),     # matches the 20.0, Δ = 0
    ])
    got = select_optimal_annual_pairs(cand, max_gap_years=2)
    assert len(got) == 1
    assert got.loc[0, "pair"] == "20220410_20230412"
    assert got.loc[0, "d_tcwv_kg_m2"] == pytest.approx(0.0)


def test_an_individually_turbulent_date_is_penalised_even_when_delta_is_small():
    """`combined = |Δtcwv| + 0.5 * (score_ref + score_sec)` — the second term
    is what stops a matched-but-stormy pair winning."""
    cand = scored([
        (2022, "2022-04-10", 20.0, 10.0),    # matched, but both turbulent
        (2022, "2022-04-14", 21.0, 0.0),
        (2023, "2023-04-12", 20.0, 10.0),
        (2023, "2023-04-16", 21.0, 0.0),
    ])
    got = select_optimal_annual_pairs(cand)
    assert got.loc[0, "pair"] == "20220414_20230416"   # Δ=0 too, but calm


def test_transitions_beyond_the_year_gap_are_not_formed():
    cand = scored([(2020, "2020-04-15", 20.0, 0.0),
                   (2023, "2023-04-15", 20.0, 0.0)])
    assert select_optimal_annual_pairs(cand, max_gap_years=2).empty


def test_every_admissible_transition_gets_a_pair():
    cand = scored([(2022, "2022-04-15", 20.0, 0.0),
                   (2023, "2023-04-15", 20.0, 0.0),
                   (2024, "2024-04-15", 20.0, 0.0)])
    got = select_optimal_annual_pairs(cand, max_gap_years=2)
    assert set(got["pair"]) == {"20220415_20230415", "20230415_20240415",
                                "20220415_20240415"}


def test_topk_returns_k_ranked_pairs_per_transition():
    """The redundancy improvement over the source paper's single pair: an
    ensemble median over k pairs absorbs residual noise."""
    cand = scored([(2022, "2022-04-10", 20.0, 0.0), (2022, "2022-04-14", 21.0, 0.0),
                   (2023, "2023-04-12", 20.0, 0.0), (2023, "2023-04-16", 22.0, 0.0)])
    got = select_topk_annual_pairs(cand, k=3)
    assert len(got) == 3
    assert list(got["rank"]) == [1, 2, 3]
    assert got["combined_atmo_score"].is_monotonic_increasing
    assert got["transition"].unique().tolist() == ["2022-2023"]


def test_topk_never_repeats_the_same_pair():
    cand = scored([(2022, "2022-04-10", 20.0, 0.0), (2023, "2023-04-12", 20.0, 0.0)])
    got = select_topk_annual_pairs(cand, k=5)
    assert len(got) == 1                      # only one couple exists
    assert got["pair"].is_unique


def test_pair_identifiers_follow_the_hyp3_convention():
    cand = scored([(2022, "2022-04-15", 20.0, 0.0), (2023, "2023-04-15", 20.0, 0.0)])
    pair = select_optimal_annual_pairs(cand).loc[0, "pair"]
    assert pair == "20220415_20230415"
    # the simpler path through build_annual_pair_list must agree
    built = build_annual_pair_list(pick_seasonal_dates(
        inventory(["2022-04-15", "2023-04-15"])))
    assert built.loc[0, "pair"] == pair
    assert list(built.columns) == ["ref_date", "sec_date", "dt_days", "pair"]


# =========================================================================
# deramp
# =========================================================================

def test_a_known_plane_is_removed_exactly():
    """Ground truth: a field that IS a plane must become zero."""
    y = np.arange(10.0)
    x = np.arange(12.0)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    field = xr.DataArray(3.0 * xx + 2.0 * yy + 7.0,
                         dims=("y", "x"), coords={"y": y, "x": x})
    mask = xr.ones_like(field, dtype=bool)

    out = deramp(field, mask)

    assert np.allclose(out["corrected"].values, 0.0, atol=1e-9)
    assert out.attrs["ramp_coef_x_per_m"] == pytest.approx(3.0)
    assert out.attrs["ramp_coef_y_per_m"] == pytest.approx(2.0)
    assert out.attrs["ramp_offset"] == pytest.approx(7.0)


def test_signal_inside_the_aoi_survives_a_ramp_fitted_outside_it():
    """The reason the fit is restricted to stable ground: real deformation must
    not be absorbed into the plane."""
    y = np.arange(10.0); x = np.arange(10.0)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    ramp = 0.5 * xx + 0.25 * yy
    signal = np.zeros_like(ramp)
    signal[4:6, 4:6] = -20.0                       # subsidence patch
    field = xr.DataArray(ramp + signal, dims=("y", "x"),
                         coords={"y": y, "x": x})
    outside = xr.DataArray(np.ones_like(ramp, dtype=bool),
                           dims=("y", "x"), coords={"y": y, "x": x})
    outside[3:7, 3:7] = False                      # fit away from the patch

    out = deramp(field, outside)
    assert out["corrected"].values[4:6, 4:6].mean() == pytest.approx(-20.0, abs=1e-6)
    assert abs(out["corrected"].values[0, 0]) < 1e-6


def test_deramp_refuses_to_fit_a_plane_on_too_few_pixels():
    """Three points define a plane exactly; fitting one would produce a
    confident, meaningless correction. It raises instead."""
    y = np.arange(5.0); x = np.arange(5.0)
    field = xr.DataArray(np.zeros((5, 5)), dims=("y", "x"),
                         coords={"y": y, "x": x})
    mask = xr.zeros_like(field, dtype=bool)
    mask[0, 0] = True
    with pytest.raises(ValueError, match="pas assez"):
        deramp(field, mask)


# =========================================================================
# masks
# =========================================================================

def test_permanent_water_is_excluded_from_the_interpretable_surface():
    """Class 5 has no coherent C-band phase; leaving it in produced pixels
    showing a rate that is physically impossible."""
    cls = xr.DataArray(np.array([[1, 5], [3, 4]]), dims=("y", "x"),
                       coords={"y": [0, 1], "x": [0, 1]})
    keep = surface_valid_mask(cls)
    assert keep.values.tolist() == [[True, False], [True, True]]


def test_the_coherence_floor_is_applied_as_well():
    cls = xr.DataArray(np.array([[1, 1]]), dims=("y", "x"),
                       coords={"y": [0], "x": [0, 1]})
    corr = xr.DataArray(np.array([[0.5, 0.1]]), dims=("y", "x"),
                        coords={"y": [0], "x": [0, 1]})
    keep = surface_valid_mask(cls, corr, corr_min=0.20)
    assert keep.values.tolist() == [[True, False]]


def test_dropping_the_intermittent_class_is_opt_in():
    cls = xr.DataArray(np.array([[4]]), dims=("y", "x"),
                       coords={"y": [0], "x": [0]})
    assert surface_valid_mask(cls).item() is np.True_ or surface_valid_mask(cls).item()
    assert not surface_valid_mask(cls, drop_classes=(4, 5)).item()


def test_the_ring_grows_the_aoi_by_the_requested_distance():
    """Fitting the ramp on a near ring rather than distant corners: a plane
    calibrated 2 km away models the local turbulent atmosphere badly."""
    aoi = xr.DataArray(np.zeros((7, 7), dtype=bool), dims=("y", "x"),
                       coords={"y": np.arange(7), "x": np.arange(7)})
    aoi[3, 3] = True
    ring = near_aoi_ring(aoi, max_dist_px=2)
    assert ring[3, 3].item()
    assert ring[3, 5].item()          # 2 px away, inside
    assert not ring[3, 6].item()      # 3 px away, outside
    assert ring.sum().item() > aoi.sum().item()


# =========================================================================
# atmosphere
# =========================================================================

def era5_tcwv(times, values) -> xr.Dataset:
    return xr.Dataset(
        {"tcwv": (("time", "latitude", "longitude"),
                  np.asarray(values, dtype=float).reshape(-1, 1, 1))},
        coords={"time": pd.to_datetime(times),
                "latitude": [52.7632], "longitude": [16.3098]})


def test_zenith_delay_is_proportional_to_water_vapour():
    era5 = era5_tcwv(["2024-04-15T17:00", "2024-04-16T17:00"], [10.0, 20.0])
    zwd = zenith_wet_delay_mm(era5, lon=16.3098, lat=52.7632)
    assert zwd.values[1] == pytest.approx(2 * zwd.values[0])
    assert (zwd.values > 0).all()


def test_the_pair_delay_is_the_difference_between_its_two_dates():
    """Only the differential delay contaminates the interferogram."""
    era5 = era5_tcwv(["2022-04-15T17:00", "2023-04-15T17:00"], [10.0, 30.0])
    zwd = zenith_wet_delay_mm(era5, lon=16.3098, lat=52.7632)
    d = pair_delays_mm(zwd, ["20220415_20230415"])
    assert d["20220415_20230415"] == pytest.approx(
        float(zwd.values[1] - zwd.values[0]))


def test_an_identical_atmosphere_on_both_dates_gives_no_correction():
    era5 = era5_tcwv(["2022-04-15T17:00", "2023-04-15T17:00"], [18.0, 18.0])
    zwd = zenith_wet_delay_mm(era5, lon=16.3098, lat=52.7632)
    assert pair_delays_mm(zwd, ["20220415_20230415"]).iloc[0] == pytest.approx(0.0)


def test_the_default_acquisition_hour_matches_the_afternoon_overpass():
    """The Rzecin track passes at ~16:36 UTC. The old 05:00 default sampled the
    MORNING atmosphere — a systematic bias, not noise. Nearest ERA5 hour is 17:00.

    Built so the morning and afternoon values differ: picking the wrong hour
    changes the answer, which is what makes this test meaningful."""
    era5 = era5_tcwv(["2022-04-15T05:00", "2022-04-15T17:00",
                      "2023-04-15T05:00", "2023-04-15T17:00"],
                     [5.0, 25.0, 5.0, 15.0])
    zwd = zenith_wet_delay_mm(era5, lon=16.3098, lat=52.7632)

    afternoon = pair_delays_mm(zwd, ["20220415_20230415"]).iloc[0]
    morning = pair_delays_mm(zwd, ["20220415_20230415"],
                             acq_time_utc="05:00").iloc[0]

    assert afternoon < 0 and morning == pytest.approx(0.0)
    assert afternoon != pytest.approx(morning)
