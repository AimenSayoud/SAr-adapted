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


def test_a_run_leaves_a_record_with_commit_and_outputs(tmp_path, monkeypatch):
    import json
    import time

    from insar_wetlands import paths, register
    hub = tmp_path / "hub"
    (hub / "08_deliverables" / "x").mkdir(parents=True)
    monkeypatch.setattr(paths, "hub_root", lambda: hub)
    e = register.Entry(id="x_demo", kind="experiment", status="current", runs_on="mac", code="experiments/x.py",
                       outputs=["hub:08_deliverables/x"], tickets=["X-999"])
    t0 = time.time()
    (hub / "08_deliverables" / "x" / "table.csv").write_text("a\n1\n")
    rec = register.record_run(e, ["--flag"], 0, t0)
    m = json.loads((rec / "manifest.json").read_text())
    assert rec.parent == hub / "05_code" / "runs" / "x_demo"
    assert m["exit_code"] == 0 and m["args"] == ["--flag"] and m["tickets"] == ["X-999"]
    assert [f["path"].endswith("table.csv") for f in m["outputs_written"]] == [True]
    assert "commit" in m["git"] and "python" in m["environment"]


def test_a_colab_launcher_runs_its_entry_through_the_register(monkeypatch):
    from insar_wetlands import colab, register
    fake = {"x_heavy": register.Entry(id="x_heavy", kind="experiment", status="open", runs_on="colab",
                                      code="experiments/x_heavy.py", title="t", tickets=["X-100"])}
    monkeypatch.setattr(colab, "load_register", lambda: fake)
    nb = colab.launcher("x_heavy", ref="main")
    src = "".join("".join(c["source"]) for c in nb["cells"])
    assert 'REF = "main"' in src and "insar_wetlands.register run x_heavy" in src and "drive.mount" in src
