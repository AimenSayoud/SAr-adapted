"""Publication checks reject actual proposed leaks, including force-added ignored files."""
from __future__ import annotations

import importlib.util
import json
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("public_boundary", REPO / "scripts/check_public_tree.py")
boundary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(boundary)


def test_ignored_runtime_and_credentials_can_still_be_detected_if_force_added(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / ".gitignore").write_text((REPO / ".gitignore").read_text())
    names = [".env.local", "private.env", ".netrc", ".cdsapirc", "08_deliverables/private/table.csv"]
    for name in names:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder")
        assert subprocess.run(["git", "check-ignore", "--quiet", name], cwd=tmp_path).returncode == 0
    assert boundary.problems(tmp_path) == []
    subprocess.run(["git", "add", "-f", *names], cwd=tmp_path, check=True)
    errors = boundary.problems(tmp_path)
    assert all(any(name in error for error in errors) for name in names)


def test_notebook_outputs_and_secret_literals_are_reported_without_values(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "example.ipynb").write_text(json.dumps({"cells": [{"outputs": [{"text": "sample"}]}]}))
    token = "ghp_" + "x" * 36
    (tmp_path / "sample.txt").write_text(token)
    errors = boundary.problems(tmp_path)
    assert any("notebook outputs" in error for error in errors)
    assert any("credential-shaped" in error for error in errors)
    assert all(token not in error for error in errors)
