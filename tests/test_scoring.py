"""The common scorer (C-061) on synthetic truths with known answers."""
from __future__ import annotations

import subprocess

import numpy as np
import pandas as pd
import pytest

from insar_wetlands import scoring as S


def _pairs(start="2021-04-01", end="2022-12-31", every=6, truth=None, est=None):
    t = pd.date_range(start, end, freq=f"{every}D", tz="UTC")
    lev = pd.Series(truth(t) if truth else 10 * np.sin(2 * np.pi * (t.dayofyear / 365.25)), index=t)
    p = pd.DataFrame({"t1": t[:-1], "t2": t[1:]})
    p["truth_inc_mm"] = np.diff(lev.to_numpy())
    p["est_mm"] = est(p) if est else p.truth_inc_mm
    hourly = lev.reindex(pd.date_range(t[0], t[-1], freq="h", tz="UTC")).interpolate()
    return p, hourly


def test_a_perfect_estimate_scores_perfectly():
    p, lev = _pairs()
    r = S.score(p, lev, S.Protocol(min_epochs=5))
    full = r["levels"].set_index("view").loc["full span"]
    assert full.pearson_detrended == pytest.approx(1.0) and full.rmse_detrended_mm == pytest.approx(0.0, abs=1e-9)
    inc = r["increments"].iloc[0]
    assert inc.slope_through_origin == pytest.approx(1.0) and inc.within_quarter_wave == 1.0
    assert (r["winters"].difference_mm.abs() < 1e-9).all()


def test_a_half_amplitude_estimate_is_seen_in_slope_and_amplitude_ratio():
    p, lev = _pairs(est=lambda p: 0.5 * p.truth_inc_mm)
    r = S.score(p, lev, S.Protocol(min_epochs=5))
    assert r["increments"].iloc[0].slope_through_origin == pytest.approx(0.5)
    assert r["levels"].set_index("view").loc["full span"].amplitude_ratio == pytest.approx(0.5, rel=1e-6)


def test_a_planted_winter_loss_is_measured_across_the_winter_only():
    def est(p):
        e = p.truth_inc_mm.copy()
        e[(p.t1 >= "2021-12-15") & (p.t1 < "2021-12-21")] -= 7.0      # one winter pair loses 7 mm
        return e
    p, lev = _pairs(est=est)
    w = S.score(p, lev, S.Protocol(min_epochs=5))["winters"]
    assert w.set_index("winter").loc["2021/22"].difference_mm == pytest.approx(-7.0)


def test_a_missing_pair_starts_a_new_segment():
    p, lev = _pairs()
    p = p.drop(index=30).reset_index(drop=True)
    c = S.chain(p)
    assert c.seg.nunique() == 2


def test_a_relation_that_flips_sign_between_years_is_flagged():
    p, lev = _pairs(est=lambda p: np.where(p.t1.dt.year == 2022, -1.0, 1.0) * p.truth_inc_mm)
    y = S.score(p, lev, S.Protocol(min_epochs=5))["years"].set_index("year")
    assert y.loc[2021, "sign"] == 1 and y.loc[2022, "sign"] == -1


def test_datum_breaks_split_the_views():
    p, lev = _pairs()
    pr = S.Protocol(datum_breaks=[(pd.Timestamp("2022-07-05", tz="UTC"), pd.Timestamp("2022-08-08", tz="UTC"))],
                    min_epochs=5)
    views = set(S.score(p, lev, pr)["levels"].view)
    assert {"datum A", "datum B"} <= views


def test_the_lockbox_refuses_unless_asked():
    p, lev = _pairs(start="2020-04-01", end="2021-10-01")
    pr = S.Protocol(lockbox_years=(2020,), min_epochs=5)
    with pytest.raises(S.LockboxError):
        S.score(p, lev, pr)
    assert len(S.score(p, lev, pr, lockbox=True)["levels"])


def test_preregistration_must_be_committed_before_the_outputs(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "config", "user.email", "t@t"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "config", "user.name", "t"], check=True)
    pre = tmp_path / "X-999.yaml"
    pre.write_text("ticket: X-999\nquestion: q\nprimary: rmse\ncomparison: plain chain\nyears: [2021]\nbetter_if: lower\n")
    with pytest.raises(RuntimeError, match="not committed"):
        S.assert_preregistered(pre)
    subprocess.run(["git", "-C", str(tmp_path), "add", "."], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "commit", "-qm", "prereg"], check=True)
    assert S.assert_preregistered(pre)["ticket"] == "X-999"
    old = tmp_path / "out.csv"
    old.write_text("x")
    import os
    os.utime(old, (0, 0))
    with pytest.raises(RuntimeError, match="older than"):
        S.assert_preregistered(pre, [old])
    (tmp_path / "bad.yaml").write_text("ticket: X-1\n")
    with pytest.raises(ValueError):
        S.load_preregistration(tmp_path / "bad.yaml")
