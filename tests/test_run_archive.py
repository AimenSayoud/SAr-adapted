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
    extract_notebook_images,
    latest_run,
    list_runs,
    load_manifest,
    run_id,
)

# A real, minimal 1x1 transparent PNG, base64-encoded -- exactly the shape a
# `colab exec`-produced `_output.ipynb` embeds for a `plt.show()`ed figure.
ONE_PIXEL_PNG_B64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk"
    "+A8AAQUBAScY42YAAAAASUVORK5CYII="
)


def _notebook_with_outputs(cell_outputs: list[list[dict]]) -> dict:
    """A minimal nbformat structure with the given per-cell outputs."""
    return {
        "cells": [
            {"cell_type": "code", "source": [], "outputs": outputs}
            for outputs in cell_outputs
        ],
        "nbformat": 4,
        "nbformat_minor": 5,
    }


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


def test_scalar_products_do_not_crash(workspace):
    outdir, drive = workspace
    run = archive_run("phaseX", outdir, root=drive,
                      products={"count": 12, "ratio": 0.75, "active": True})
    manifest = load_manifest(run)
    assert manifest["products"]["count"]["value"] == 12
    assert manifest["products"]["count"]["is_scalar"] is True


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


# --- extracting images from an executed notebook ---------------------------

def test_extract_notebook_images_pulls_embedded_png(tmp_path: Path):
    import base64
    import json

    nb = _notebook_with_outputs([
        [{"output_type": "display_data",
          "data": {"image/png": ONE_PIXEL_PNG_B64}}],
    ])
    nb_path = tmp_path / "phaseX_output.ipynb"
    nb_path.write_text(json.dumps(nb))

    out_dir = tmp_path / "images"
    written = extract_notebook_images(nb_path, out_dir)

    assert len(written) == 1
    assert written[0] == out_dir / "cell000_output00.png"
    assert written[0].read_bytes() == base64.b64decode(ONE_PIXEL_PNG_B64)


def test_extract_notebook_images_ignores_text_only_outputs(tmp_path: Path):
    import json

    nb = _notebook_with_outputs([
        [{"output_type": "stream", "name": "stdout", "text": ["356 pairs\n"]}],
        [{"output_type": "execute_result",
          "data": {"text/plain": ["<DataFrame>"]}}],
    ])
    nb_path = tmp_path / "phaseX_output.ipynb"
    nb_path.write_text(json.dumps(nb))

    assert extract_notebook_images(nb_path, tmp_path / "images") == []


def test_extract_notebook_images_numbers_multiple_figures_in_order(tmp_path: Path):
    import json

    nb = _notebook_with_outputs([
        [{"output_type": "display_data", "data": {"image/png": ONE_PIXEL_PNG_B64}}],
        [{"output_type": "stream", "name": "stdout", "text": ["...\n"]},
         {"output_type": "display_data", "data": {"image/png": ONE_PIXEL_PNG_B64}}],
    ])
    nb_path = tmp_path / "phaseX_output.ipynb"
    nb_path.write_text(json.dumps(nb))

    written = extract_notebook_images(nb_path, tmp_path / "images")
    names = sorted(p.name for p in written)
    assert names == ["cell000_output00.png", "cell001_output01.png"]


def test_archive_run_with_executed_notebook_extracts_images(workspace):
    import json

    outdir, drive = workspace
    nb = _notebook_with_outputs([
        [{"output_type": "display_data", "data": {"image/png": ONE_PIXEL_PNG_B64}}],
    ])
    nb_path = outdir / "phaseX_output.ipynb"
    nb_path.write_text(json.dumps(nb))

    run = archive_run("phaseX", outdir, root=drive, executed_notebook=nb_path)
    manifest = load_manifest(run)

    assert manifest["images"] == ["images/cell000_output00.png"]
    assert (run / "images" / "cell000_output00.png").is_file()


def test_archive_run_reports_a_missing_executed_notebook(workspace):
    outdir, drive = workspace
    run = archive_run("phaseX", outdir, root=drive,
                      executed_notebook=outdir / "does_not_exist.ipynb")
    manifest = load_manifest(run)
    assert manifest["images"] == [f"MISSING: {outdir / 'does_not_exist.ipynb'}"]


def test_archive_run_without_executed_notebook_has_empty_images(workspace):
    outdir, drive = workspace
    run = archive_run("phaseX", outdir, root=drive)
    assert load_manifest(run)["images"] == []
