"""Every name a live notebook uses must be defined somewhere in that notebook.

Why this test exists. The 2026-09-05 bootstrap migration (3f685a5) replaced each
notebook's setup cell with `ctx = start(...)`, but deleted imports and aliases
(`CROPPED`, `list_pairs`, `load_worldcover`, `json`, …) that later cells still
used. `make phases` checked the preamble and passed; nothing checked the body.
Twenty-two notebooks could not run end to end until the web-atlas audit
(2026-09-24) ran them locally and restored the names. ruff's F821 sees the
whole notebook as one module, which is exactly the question.
"""

import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]

# `bootstrap.start()` injects these into the caller's globals at run time.
INJECTED = {"np", "pd", "plt", "xr"}

# Known gaps that need a person, not a restored import. Each carries its reason;
# an entry here that no longer fires fails the test, so the list cannot rot.
KNOWN = {
    ("phase01", "coverage"): "logic gap: a coverage table the migration removed; acquisition phase, needs Aymen",
    ("phase01", "dates"): "logic gap: same cell as `coverage`",
    ("phase01", "pairs"): "archive cell refers to a variable this phase never defines",
    ("phase02b", "pairs"): "the topk->pairs selection cell was removed; HyP3 submission phase",
    ("phaseB", "EGMS"): "an EGMS product constant never defined; exploratory phase",
}


def _ruff() -> str | None:
    return shutil.which("ruff") or next(
        (str(p) for p in Path.home().glob("Library/Python/*/bin/ruff")), None)


@pytest.mark.skipif(_ruff() is None, reason="ruff not installed")
def test_live_notebooks_define_every_name_they_use():
    phases = yaml.safe_load((REPO / "config" / "phases.yaml").read_text())
    live = {k: v["notebook"] for k, v in phases.items() if v.get("status") != "superseded"}
    out = subprocess.run([_ruff(), "check", "--select", "F821", "--output-format", "json",
                          "--no-cache", *live.values()],
                         cwd=REPO, capture_output=True, text=True)
    import json
    found = set()
    by_path = {v: k for k, v in live.items()}
    for d in json.loads(out.stdout or "[]"):
        name = d["message"].split("`")[1]
        if name in INJECTED:
            continue
        rel = str(Path(d["filename"]).resolve().relative_to(REPO))
        found.add((by_path[rel], name))
    unexpected = sorted(found - set(KNOWN))
    assert not unexpected, f"undefined names in live notebooks: {unexpected}"
    stale = sorted(set(KNOWN) - found)
    assert not stale, f"KNOWN entries that no longer fire — remove them: {stale}"
