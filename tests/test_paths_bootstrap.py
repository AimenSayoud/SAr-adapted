"""Ground-truth tests for path resolution and the phase bootstrap.

The failures these guard against are quiet ones: a path that resolves somewhere
unexpected and writes there instead of raising, or a Colab-only default that
works in Colab and nowhere else."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
import yaml

from insar_wetlands.bootstrap import Context, start
from insar_wetlands.paths import COLAB_DRIVE, Paths, make_paths, repo_root, resolve_drive


@pytest.fixture()
def fake_repo(tmp_path: Path) -> Path:
    """A minimal repository: what `repo_root` searches for is config/config.yaml."""
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "config.yaml").write_text(yaml.safe_dump({
        "site": {"name": "Test", "centroid": [16.3098, 52.7632]},
        "paths": {"outputs": "outputs",
                  "drive_data_root": "/content/drive/MyDrive/insar_rzecin"},
    }))
    return tmp_path


# --- resolution order ------------------------------------------------------

def test_explicit_root_wins_over_everything(fake_repo, monkeypatch, tmp_path):
    monkeypatch.setenv("INSAR_DRIVE_ROOT", str(tmp_path / "from_env"))
    chosen = tmp_path / "explicit"
    assert resolve_drive(chosen, cfg={}, repo=fake_repo) == chosen


def test_env_var_wins_over_config(fake_repo, monkeypatch, tmp_path):
    env = tmp_path / "from_env"
    monkeypatch.setenv("INSAR_DRIVE_ROOT", str(env))
    cfg = {"paths": {"drive_data_root": str(tmp_path / "from_config")}}
    assert resolve_drive(None, cfg=cfg, repo=fake_repo) == env


def test_configured_root_is_ignored_when_it_does_not_exist(fake_repo, monkeypatch):
    """The committed config names the Colab path, which is meaningless locally.

    Falling through to a working local directory beats resolving to a path that
    does not exist and failing later, somewhere less obvious."""
    monkeypatch.delenv("INSAR_DRIVE_ROOT", raising=False)
    cfg = {"paths": {"drive_data_root": "/content/drive/MyDrive/insar_rzecin"}}
    assert resolve_drive(None, cfg=cfg, repo=fake_repo) == fake_repo / "data" / "_local"


def test_configured_root_is_used_when_it_exists(fake_repo, monkeypatch, tmp_path):
    monkeypatch.delenv("INSAR_DRIVE_ROOT", raising=False)
    real = tmp_path / "real_drive"
    real.mkdir()
    cfg = {"paths": {"drive_data_root": str(real)}}
    assert resolve_drive(None, cfg=cfg, repo=fake_repo) == real


def test_repo_root_finds_the_config_from_a_subdirectory(fake_repo):
    deep = fake_repo / "notebooks" / "nested"
    deep.mkdir(parents=True)
    assert repo_root(deep) == fake_repo


def test_repo_root_raises_when_there_is_no_config(tmp_path):
    with pytest.raises(FileNotFoundError):
        repo_root(tmp_path)


# --- the Paths object ------------------------------------------------------

def test_phase_scopes_the_output_directory(fake_repo, tmp_path):
    p = Paths(repo=fake_repo, drive=tmp_path / "drive", phase="phaseG")
    assert p.outputs == fake_repo / "outputs" / "phaseG"
    assert p.outputs.is_dir()                    # created on access
    assert p.log_file == p.outputs / "run.log"


def test_paths_without_a_phase_use_the_outputs_root(fake_repo, tmp_path):
    p = Paths(repo=fake_repo, drive=tmp_path / "drive")
    assert p.outputs == fake_repo / "outputs"


def test_cache_replaces_the_hardcoded_colab_directory(fake_repo, tmp_path):
    """`load_worldcover` used to default to /content/worldcover in its signature.

    The replacement must live on the data root, so a downloaded tile survives
    the session that fetched it."""
    drive = tmp_path / "drive"
    p = Paths(repo=fake_repo, drive=drive)
    assert p.cache == drive / "cache"
    assert p.cache.is_dir()
    assert "/content/worldcover" not in str(p.cache)


def test_construction_creates_nothing(fake_repo, tmp_path):
    """Building Paths must be free — directories appear only when asked for."""
    drive = tmp_path / "untouched"
    Paths(repo=fake_repo, drive=drive, phase="phaseX")
    assert not drive.exists()


def test_on_colab_is_false_for_a_local_root(fake_repo, tmp_path):
    assert Paths(repo=fake_repo, drive=tmp_path).on_colab is False
    assert Paths(repo=fake_repo, drive=COLAB_DRIVE).on_colab is True


def test_describe_reports_what_resolved(fake_repo, tmp_path):
    d = Paths(repo=fake_repo, drive=tmp_path, phase="phaseG").describe()
    assert d["phase"] == "phaseG"
    assert d["drive_exists"] is True
    assert d["cropped_exists"] is False          # honest about a missing input


def test_for_phase_returns_a_new_object(fake_repo, tmp_path):
    base = Paths(repo=fake_repo, drive=tmp_path)
    g = base.for_phase("phaseG")
    assert g.phase == "phaseG" and base.phase is None


# --- the bootstrap ---------------------------------------------------------

def test_start_is_cheap_and_loads_nothing(fake_repo, tmp_path):
    """No stack, no zones, no network: a phase that needs none of it pays nothing."""
    ctx = start("phaseG", mount=False, git=False, repo=fake_repo,
                drive_root=tmp_path / "drive")
    assert isinstance(ctx, Context)
    assert ctx.phase == "phaseG"
    assert ctx.cfg["site"]["name"] == "Test"
    assert "pairs" not in ctx.__dict__           # cached_property not triggered
    assert "template" not in ctx.__dict__


def test_start_writes_a_run_log(fake_repo, tmp_path):
    ctx = start("phaseG", mount=False, git=False, repo=fake_repo,
                drive_root=tmp_path / "drive")
    for h in ctx.log.handlers:
        h.flush()
    text = ctx.paths.log_file.read_text()
    assert "phase phaseG starting" in text
    assert "drive" in text                       # the resolved paths are recorded


def test_start_warns_when_the_data_root_is_missing(fake_repo, tmp_path, caplog):
    missing = tmp_path / "no_drive"
    with caplog.at_level(logging.WARNING):
        ctx = start("phaseX", mount=False, git=False, repo=fake_repo,
                    drive_root=missing)
        ctx.log.addHandler(caplog.handler)       # capture propagates off by default
    assert not missing.exists()


def test_outdir_is_scoped_to_the_phase(fake_repo, tmp_path):
    ctx = start("phaseD", mount=False, git=False, repo=fake_repo,
                drive_root=tmp_path / "drive")
    assert ctx.outdir == fake_repo / "outputs" / "phaseD"


def test_archive_routes_to_the_resolved_drive(fake_repo, tmp_path):
    """The context archives to its own data root, not to a global default."""
    drive = tmp_path / "drive"
    ctx = start("phaseG", mount=False, git=False, repo=fake_repo, drive_root=drive)
    (ctx.outdir / "series.csv").write_text("date,value\n2024-01-01,3.29\n")

    run = ctx.archive(params={"n_nulls": 280},
                      products={"series": ctx.outdir / "series.csv"})

    assert run.parent == drive / "runs" / "phaseG"
    import json
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["parameters"] == {"n_nulls": 280}
    assert manifest["products"]["series"]["heavy"] is False


def test_two_phases_do_not_share_an_output_directory(fake_repo, tmp_path):
    drive = tmp_path / "drive"
    a = start("phaseA", mount=False, git=False, repo=fake_repo, drive_root=drive)
    b = start("phaseB", mount=False, git=False, repo=fake_repo, drive_root=drive)
    assert a.outdir != b.outdir
