"""C-028: ctx.archive() attaches this run's executed notebook without being told to."""
import logging
import os
import time
from pathlib import Path
from types import SimpleNamespace

from insar_wetlands import run_archive
from insar_wetlands.bootstrap import Context

PHASES = """phaseX:
  title: t
  notebook: notebooks/04_hypotheses/phaseX_test.ipynb
  group: hypotheses
  question: q
  status: current
  inputs: []
  outputs: []
  paper: "-"
"""


def _ctx(repo: Path) -> Context:
    (repo / "config").mkdir(parents=True)
    (repo / "config" / "phases.yaml").write_text(PHASES)
    (repo / "notebooks" / "04_hypotheses").mkdir(parents=True)
    paths = SimpleNamespace(repo=repo, drive=repo / "drive", outputs=repo / "outputs" / "phaseX")
    return Context(phase="phaseX", cfg={}, paths=paths, log=logging.getLogger("t"))


def test_finds_the_output_notebook_written_during_this_run(tmp_path):
    ctx = _ctx(tmp_path)
    out = tmp_path / "notebooks/04_hypotheses/phaseX_test_output.ipynb"
    out.write_text("{}")
    assert ctx.executed_notebook() == out


def test_ignores_an_output_notebook_from_an_earlier_run(tmp_path):
    ctx = _ctx(tmp_path)
    out = tmp_path / "notebooks/04_hypotheses/phaseX_test_output.ipynb"
    out.write_text("{}")
    old = ctx.started - 3600
    os.utime(out, (old, old))
    assert ctx.executed_notebook() is None


def test_unknown_phase_gives_none(tmp_path):
    ctx = _ctx(tmp_path)
    ctx.phase = "not_declared"
    assert ctx.executed_notebook() is None


def test_archive_forwards_it(tmp_path, monkeypatch):
    ctx = _ctx(tmp_path)
    out = tmp_path / "notebooks/04_hypotheses/phaseX_test_output.ipynb"
    time.sleep(0.01)
    out.write_text("{}")
    seen = {}
    monkeypatch.setattr(run_archive, "archive_run", lambda *a, **k: seen.update(k) or tmp_path / "run")
    ctx.archive()
    assert seen["executed_notebook"] == out
