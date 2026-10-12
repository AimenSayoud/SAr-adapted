"""Data roots (paths.cube_root, s1_local_root, s1_crops): env first, then the hub, then the Drive backup on Colab."""
from __future__ import annotations

from insar_wetlands import paths


def test_env_wins(monkeypatch, tmp_path):
    monkeypatch.setenv(paths.CUBE_ENV, str(tmp_path / "c"))
    monkeypatch.setenv(paths.S1_ENV, str(tmp_path / "s"))
    assert paths.cube_root() == tmp_path / "c"
    assert paths.s1_crops("2025", "descending") == tmp_path / "s" / "s1_2025" / "hyp3_cropped_descending"


def test_without_a_hub_the_drive_backup_is_used(monkeypatch, tmp_path):
    for v in (paths.CUBE_ENV, paths.S1_ENV, paths.HUB_ENV):
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setattr(paths, "hub_root", lambda: None)
    monkeypatch.setattr(paths, "MAC_BACKUP", tmp_path / "mac_backup")
    assert paths.cube_root() == tmp_path / "mac_backup" / "06_data" / "cube"
    assert paths.s1_local_root() == tmp_path / "mac_backup" / "05_code" / "local"


def test_2022_2024_crops_are_in_the_drive_working_area(tmp_path):
    assert paths.s1_crops("2022_2024", "ascending", drive=tmp_path) == tmp_path / "hyp3_cropped"
    assert paths.s1_crops("2022_2024", "descending", drive=tmp_path) == tmp_path / "hyp3_cropped_descending"
