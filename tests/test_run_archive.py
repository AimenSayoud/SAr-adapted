"""Ground-truth tests for run archiving.

Each test constructs a known situation and asserts the archive reports it
correctly — the same shape as the other test modules here."""

from __future__ import annotations

from pathlib import Path

import pytest

from insar_wetlands.run_archive import (
    archive_run,
    compare_runs,
    drive_root,
    latest_run,
    list_runs,
    load_manifest,
    run_id,
)


@pytest.fixture()
def workspace(tmp_path: Path) -> tuple[Path, Path]:
    """An output directory holding one light and one heavy product."""
    outdir = tmp_path / "outputs" / "phaseX"
    outdir.mkdir(parents=True)
    (outdir / "series.csv").write_text("date,value\n2024-01-01,3.29\n")
    (outdir / "stack.nc").write_bytes(b"\x00" * 2048)   # heavy by suffix
    return outdir, tmp_path / "drive"


def test_run_id_is_sortable_and_carries_the_commit():
    from datetime import datetime, timezone
    early = run_id(datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc), "aaaaaaa")
    late = run_id(datetime(2026, 9, 5, 17, 12, 0, tzinfo=timezone.utc), "b4995f9")
    assert early < late                      # chronological by string sort
    assert late.endswith("_b4995f9")
    assert late.startswith("20260905T171200Z")


def test_heavy_products_are_described_but_not_copied(workspace):
    outdir, drive = workspace
    run = archive_run("phaseX", outdir, root=drive,
                      products={"stack": outdir / "stack.nc",
                                "series": outdir / "series.csv"})
    manifest = load_manifest(run)

    assert manifest["products"]["stack"]["heavy"] is True
    assert manifest["products"]["series"]["heavy"] is False

    # the light file is copied, the heavy one is not
    assert (run / "products" / "series.csv").exists()
    assert not (run / "products" / "stack.nc").exists()

    # but the heavy file is still identified well enough to detect substitution
    assert manifest["products"]["stack"]["bytes"] == 2048
    assert manifest["products"]["stack"]["sha256_head"]


def test_a_missing_product_is_recorded_rather_than_raising(workspace):
    outdir, drive = workspace
    run = archive_run("phaseX", outdir, root=drive,
                      products={"absent": outdir / "never_written.csv"})
    assert load_manifest(run)["products"]["absent"]["exists"] is False


def test_runs_accumulate_and_are_never_overwritten(workspace, monkeypatch):
    outdir, drive = workspace
    import insar_wetlands.run_archive as ra

    monkeypatch.setattr(ra, "run_id", lambda **kw: "20260101T000000Z_aaaaaaa")
    ra.archive_run("phaseX", outdir, root=drive)
    with pytest.raises(FileExistsError):
        ra.archive_run("phaseX", outdir, root=drive)

    monkeypatch.setattr(ra, "run_id", lambda **kw: "20260102T000000Z_bbbbbbb")
    ra.archive_run("phaseX", outdir, root=drive)

    runs = list_runs("phaseX", root=drive)
    assert len(runs) == 2
    assert runs[0].name < runs[1].name              # oldest first
    assert latest_run("phaseX", root=drive) == runs[1]


def test_environment_and_git_state_are_captured(workspace):
    outdir, drive = workspace
    manifest = load_manifest(archive_run("phaseX", outdir, root=drive))
    assert manifest["environment"]["python"]
    assert "numpy" in manifest["environment"]["packages"]
    assert set(manifest["git"]) == {"commit", "short", "branch", "dirty"}


def test_compare_runs_names_what_changed(workspace, monkeypatch):
    outdir, drive = workspace
    import insar_wetlands.run_archive as ra

    monkeypatch.setattr(ra, "run_id", lambda **kw: "20260101T000000Z_aaaaaaa")
    ra.archive_run("phaseX", outdir, root=drive, params={"n_nulls": 92})
    monkeypatch.setattr(ra, "run_id", lambda **kw: "20260102T000000Z_bbbbbbb")
    ra.archive_run("phaseX", outdir, root=drive, params={"n_nulls": 280})

    diff = compare_runs("phaseX", root=drive)
    assert diff["comparable"] is True
    assert diff["parameters_changed"] == {"n_nulls": (92, 280)}


def test_compare_runs_is_honest_about_a_single_run(workspace):
    outdir, drive = workspace
    archive_run("phaseX", outdir, root=drive)
    assert compare_runs("phaseX", root=drive) == {"comparable": False, "n_runs": 1}


def test_drive_root_falls_back_when_drive_is_not_mounted(monkeypatch):
    monkeypatch.delenv("INSAR_DRIVE_ROOT", raising=False)
    monkeypatch.setattr("insar_wetlands.run_archive.DEFAULT_DRIVE_ROOT",
                        "/definitely/not/mounted")
    assert drive_root() == Path("outputs") / "_runs"


def test_env_var_overrides_the_colab_mount(monkeypatch, tmp_path):
    monkeypatch.setenv("INSAR_DRIVE_ROOT", str(tmp_path))
    assert drive_root() == tmp_path


def test_describe_file_flags_a_large_light_file_as_heavy(tmp_path):
    import insar_wetlands.run_archive as ra
    big = tmp_path / "huge.csv"
    big.write_bytes(b"x" * 128)
    monkeypatch_limit = 64
    original = ra.MAX_COPY_BYTES
    try:
        ra.MAX_COPY_BYTES = monkeypatch_limit
        assert ra.describe_file(big)["heavy"] is True
    finally:
        ra.MAX_COPY_BYTES = original
