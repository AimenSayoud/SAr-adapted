"""The register (config/register.yaml): loads, validates, finds; the real one is sound."""
from __future__ import annotations

from pathlib import Path

from insar_wetlands.register import find, load_register, producers_consumers, untracked_scripts, validate

REPO = Path(__file__).resolve().parents[1]

SAMPLE = """
version: 1
entries:
  phase08:
    kind: phase
    status: finished
    runs_on: colab
    code: notebooks/gone.ipynb
    code_ref: some-tag
    outputs: [ts_sbas.nc]
    tries: [T13]
  x100_slc:
    kind: experiment
    status: current
    runs_on: mac
    code: experiments/x100.py
    command: python experiments/x100.py
    inputs: [phase08.ts_sbas.nc]
    tickets: [X-100]
    chains: [slc-pilot]
  broken:
    kind: notebook
    status: done
    runs_on: cloud
    code: experiments/missing.py
"""


def _repo(tmp_path):
    (tmp_path / "experiments").mkdir()
    (tmp_path / "experiments" / "x100.py").write_text("")
    (tmp_path / "experiments" / "stray.py").write_text("")
    p = tmp_path / "register.yaml"
    p.write_text(SAMPLE)
    return p


def test_validate_reports_bad_enums_and_missing_code(tmp_path):
    reg = load_register(_repo(tmp_path))
    bad = validate(reg, repo=tmp_path, check_git=False)
    assert any("broken: kind" in b for b in bad)
    assert any("broken: status" in b for b in bad)
    assert any("broken: runs_on" in b for b in bad)
    assert any("missing.py is not in the tree" in b for b in bad)
    assert not any(b.startswith(("phase08", "x100_slc")) for b in bad)


def test_find_and_provenance(tmp_path):
    reg = load_register(_repo(tmp_path))
    assert [e.id for e in find(reg, ticket="X-100")] == ["x100_slc"]
    assert [e.id for e in find(reg, chain="slc-pilot")] == ["x100_slc"]
    assert [e.id for e in find(reg, try_id="T13")] == ["phase08"]
    prod, cons, _ = producers_consumers(reg)
    assert prod["ts_sbas.nc"] == ["phase08"] and cons["ts_sbas.nc"] == ["x100_slc"]
    assert untracked_scripts(reg, repo=tmp_path) == ["experiments/stray.py"]


def test_the_real_register_is_sound():
    reg = load_register()
    assert validate(reg, check_git=False) == []
    assert untracked_scripts(reg) == []
