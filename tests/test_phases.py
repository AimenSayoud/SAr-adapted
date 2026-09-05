"""Tests for the declared pipeline.

Two jobs. The unit tests below build small synthetic graphs and check that each
kind of inconsistency is caught. The last test runs the validator against the
project's real `config/phases.yaml`, which is what stops the declaration drifting
away from the notebooks on disk — the failure that left the README documenting
14 of 38 phases."""

from __future__ import annotations

from pathlib import Path

import yaml

from insar_wetlands.phases import (
    GROUPS,
    execution_order,
    find_cycles,
    load_phases,
    missing_inputs,
    readme_table,
    summary,
    validate,
)

REPO = Path(__file__).resolve().parents[1]


def write(tmp_path: Path, spec: dict) -> Path:
    p = tmp_path / "phases.yaml"
    p.write_text(yaml.safe_dump(spec))
    return p


# --- parsing ---------------------------------------------------------------

def test_a_phase_reads_its_declared_fields(tmp_path):
    p = load_phases(write(tmp_path, {
        "phaseG": {"title": "Aggregation", "notebook": "notebooks/g.ipynb",
                   "question": "Does it work?", "status": "current",
                   "inputs": ["phaseD.coh.csv", "water_mask.nc"],
                   "outputs": ["series_AC.csv"], "paper": "§4.3"},
    }))["phaseG"]
    assert p.title == "Aggregation"
    assert p.upstream() == ["phaseD"]        # bare artefacts are not phases
    assert p.is_live is True


def test_superseded_phases_are_not_live(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phase13": {"title": "old", "notebook": "n.ipynb", "status": "superseded"},
        "phaseI": {"title": "new", "notebook": "n.ipynb", "status": "current"},
    }))
    assert ps["phase13"].is_live is False
    assert ps["phaseI"].is_live is True


# --- validation ------------------------------------------------------------

def test_unknown_upstream_phase_is_caught(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "notebook": "n.ipynb", "inputs": ["phaseX.out.csv"]},
    }))
    assert any("unknown phase" in p for p in validate(ps))


def test_a_live_phase_may_not_read_from_a_dead_one(tmp_path):
    """The point of marking a phase superseded is that nothing still uses it."""
    ps = load_phases(write(tmp_path, {
        "phase13": {"title": "old", "notebook": "n.ipynb", "status": "superseded",
                    "outputs": ["hydro.csv"]},
        "phaseG": {"title": "g", "notebook": "n.ipynb", "status": "current",
                   "inputs": ["phase13.hydro.csv"]},
    }))
    assert any("dead one" in p for p in validate(ps))


def test_claiming_to_supersede_a_live_phase_is_caught(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phaseE": {"title": "old", "notebook": "n.ipynb", "status": "current"},
        "phaseE2": {"title": "new", "notebook": "n.ipynb", "supersedes": ["phaseE"]},
    }))
    assert any("still" in p and "phaseE" in p for p in validate(ps))


def test_unknown_status_is_caught(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "notebook": "n.ipynb", "status": "maybe"},
    }))
    assert any("unknown status" in p for p in validate(ps))


def test_a_missing_notebook_is_caught(tmp_path):
    (tmp_path / "notebooks").mkdir()
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "notebook": "notebooks/absent.ipynb"},
    }))
    assert any("notebook not found" in p for p in validate(ps, repo=tmp_path))


def test_an_undeclared_notebook_is_caught(tmp_path):
    """The failure that left the README covering 14 of 38."""
    (tmp_path / "notebooks").mkdir()
    (tmp_path / "notebooks" / "declared.ipynb").write_text("{}")
    (tmp_path / "notebooks" / "forgotten.ipynb").write_text("{}")
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "notebook": "notebooks/declared.ipynb"},
    }))
    problems = validate(ps, repo=tmp_path)
    assert any("undeclared notebook" in p and "forgotten" in p for p in problems)


def test_a_cycle_is_reported_not_hung_on(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phaseA": {"title": "a", "notebook": "n.ipynb", "inputs": ["phaseB.x.csv"]},
        "phaseB": {"title": "b", "notebook": "n.ipynb", "inputs": ["phaseA.y.csv"]},
    }))
    assert find_cycles(ps)
    assert any("cycle" in p for p in validate(ps))


# --- groups ----------------------------------------------------------------

def test_a_phase_without_a_group_is_caught(tmp_path):
    """Group is what the two filename schemes never carried."""
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "notebook": "notebooks/04_hypotheses/g.ipynb"},
    }))
    assert any("no group declared" in p for p in validate(ps))


def test_an_unknown_group_is_caught(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "group": "miscellaneous",
                   "notebook": "notebooks/04_hypotheses/g.ipynb"},
    }))
    assert any("unknown group" in p for p in validate(ps))


def test_a_notebook_filed_outside_its_group_directory_is_caught(tmp_path):
    """The declaration and the filesystem are two statements of one fact. If
    they may disagree, the flat layout grows back one notebook at a time."""
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "group": "hypotheses",
                   "notebook": "notebooks/02_inversion/g.ipynb"},
    }))
    assert any("belongs in notebooks/04_hypotheses/" in p for p in validate(ps))


def test_group_is_orthogonal_to_status(tmp_path):
    """A superseded phase keeps the group of the question it asked, so it is
    filed next to the phase that replaced it rather than in a graveyard."""
    ps = load_phases(write(tmp_path, {
        "phaseC2": {"title": "old", "group": "inversion", "status": "superseded",
                    "notebook": "notebooks/02_inversion/c2.ipynb"},
        "phaseE2": {"title": "new", "group": "inversion", "status": "current",
                    "notebook": "notebooks/02_inversion/e2.ipynb",
                    "supersedes": ["phaseC2"]},
    }))
    assert ps["phaseC2"].group == ps["phaseE2"].group
    assert ps["phaseC2"].is_live is False


def test_every_real_phase_declares_a_known_group():
    phases = load_phases(REPO / "config" / "phases.yaml")
    for name, p in phases.items():
        assert p.group in GROUPS, f"{name} has group {p.group!r}"


def test_the_readme_table_is_grouped_and_still_covers_everything():
    phases = load_phases(REPO / "config" / "phases.yaml")
    table = readme_table(phases)
    for directory in GROUPS.values():
        assert f"### {directory}" in table
    for name in phases:
        assert f"`{name}`" in table


# --- ordering --------------------------------------------------------------

def test_execution_order_respects_dependencies(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phaseG": {"title": "g", "notebook": "n.ipynb", "inputs": ["phaseD.c.csv"]},
        "phaseD": {"title": "d", "notebook": "n.ipynb", "inputs": ["phase02.h.nc"]},
        "phase02": {"title": "h", "notebook": "n.ipynb"},
    }))
    order = execution_order(ps)
    assert order.index("phase02") < order.index("phaseD") < order.index("phaseG")


def test_execution_order_is_stable(tmp_path):
    """An unstable order would make the generated README churn on every build."""
    spec = {f"phase{i:02d}": {"title": str(i), "notebook": "n.ipynb"} for i in range(6)}
    ps = load_phases(write(tmp_path, spec))
    assert execution_order(ps) == execution_order(ps)


def test_superseded_phases_are_excluded_from_the_live_order(tmp_path):
    ps = load_phases(write(tmp_path, {
        "phase13": {"title": "old", "notebook": "n.ipynb", "status": "superseded"},
        "phaseI": {"title": "new", "notebook": "n.ipynb"},
    }))
    assert execution_order(ps) == ["phaseI"]
    assert "phase13" in execution_order(ps, live_only=False)


# --- inputs on disk --------------------------------------------------------

def test_missing_inputs_are_listed_before_a_run(tmp_path):
    drive = tmp_path / "drive"
    (drive).mkdir()
    (drive / "water_mask.nc").write_bytes(b"")
    ps = load_phases(write(tmp_path, {
        "phaseD": {"title": "d", "notebook": "n.ipynb", "outputs": ["coh.csv"]},
        "phaseG": {"title": "g", "notebook": "n.ipynb",
                   "inputs": ["water_mask.nc", "phaseD.coh.csv"]},
    }))
    missing = missing_inputs(ps["phaseG"], ps, drive)
    assert missing == ["phaseD.coh.csv"]     # the present one is not reported


# --- the real file ---------------------------------------------------------

def test_the_projects_own_declaration_is_sound():
    """Runs against config/phases.yaml. Fails if a notebook is added without
    being declared, or if a dependency points at a superseded phase."""
    phases = load_phases(REPO / "config" / "phases.yaml")
    problems = validate(phases, repo=REPO)
    assert not problems, "\n".join(problems)


def test_every_notebook_on_disk_is_declared():
    phases = load_phases(REPO / "config" / "phases.yaml")
    declared = {p.notebook for p in phases.values()}
    on_disk = {str(p.relative_to(REPO))
               for p in (REPO / "notebooks").glob("**/*.ipynb")}
    assert on_disk == declared


def test_readme_table_covers_every_phase():
    phases = load_phases(REPO / "config" / "phases.yaml")
    table = readme_table(phases)
    for name in phases:
        assert f"`{name}`" in table
    assert summary(phases)["n_phases"] == len(phases)
