"""C-033: HyP3 credits are spent only on an explicit, count-exact confirmation."""
import sys
import types

import pandas as pd

from insar_wetlands.hyp3 import jobs


class FakeJob:
    def __init__(self, i):
        self.job_id = f"job{i}"


class FakeHyP3:
    submitted: list = []

    def __init__(self):
        pass

    def find_jobs(self, name=None):
        return []

    def prepare_insar_isce_burst_job(self, ref, sec, **kw):
        return {"ref": ref, "sec": sec}

    def submit_prepared_jobs(self, prepared_jobs):
        FakeHyP3.submitted += prepared_jobs
        return [FakeJob(i) for i in range(len(prepared_jobs))]


def _setup(monkeypatch):
    FakeHyP3.submitted = []
    monkeypatch.setitem(sys.modules, "hyp3_sdk", types.SimpleNamespace(HyP3=FakeHyP3))
    monkeypatch.setattr(jobs, "_existing_granule_pairs", lambda hyp3, name: set())
    dates = pd.to_datetime(["2022-01-01", "2022-01-13", "2022-01-25"])
    pairs = pd.DataFrame({"ref_date": dates[:2], "sec_date": dates[1:], "pair": ["p1", "p2"]})
    granules = {d: f"G{i}" for i, d in enumerate(dates)}
    return pairs, granules


def test_no_confirmation_spends_nothing(monkeypatch):
    pairs, granules = _setup(monkeypatch)
    out = jobs.submit_pairs(pairs, granules, "rzecin_test")
    assert FakeHyP3.submitted == []
    assert set(out.status) == {"NEEDS_CONFIRMATION"}


def test_a_stale_confirmation_for_another_batch_spends_nothing(monkeypatch):
    pairs, granules = _setup(monkeypatch)
    jobs.submit_pairs(pairs, granules, "rzecin_test", confirm="spend 5 credits on rzecin_test")
    jobs.submit_pairs(pairs, granules, "rzecin_test", confirm="spend 2 credits on other_name")
    assert FakeHyP3.submitted == []


def test_the_exact_phrase_submits(monkeypatch):
    pairs, granules = _setup(monkeypatch)
    out = jobs.submit_pairs(pairs, granules, "rzecin_test", confirm=jobs.confirmation_phrase(2, "rzecin_test"))
    assert len(FakeHyP3.submitted) == 2
    assert set(out.status) == {"PENDING"}
