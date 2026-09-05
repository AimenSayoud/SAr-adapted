"""Archive every phase execution to Drive, so a run can be reconstructed later.

Why this module exists. `utils.manifest.write_manifest` records what a phase
produced, but it writes to the Colab disk — which is erased between sessions —
and it overwrites, so only the most recent run survives. When the numbers in the
paper moved between the July draft and the round-2 corrections, there was no
record of *which execution* had produced the earlier values, and the Drive copy
of the data silently fell behind the repository.

An archived run answers three questions that a manifest alone cannot:

  * *When was this produced, and from which commit?* — run id and git SHA.
  * *In what environment?* — Python and key package versions, captured at run time.
  * *What did the previous run give?* — runs are never overwritten.

Usage is one line at the end of a phase notebook::

    from insar_wetlands.run_archive import archive_run
    archive_run("phaseG", outdir, params=PARAMS, products=PRODUCTS)

Heavy products (``.nc``, ``.tif``, ``.h5``) are **recorded but not copied**: they
already live on Drive and duplicating them per run would exhaust the quota. They
are listed with their path, size and modification time, which is enough to tell
whether a later run read the same input.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

# Where the Drive is mounted in Colab. Overridable for local runs and tests.
DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/insar_rzecin"
ENV_VAR = "INSAR_DRIVE_ROOT"

# Extensions that stay where they are. Copying them per run would multiply
# gigabytes for no gain: the manifest records enough to identify them.
HEAVY_SUFFIXES = {".nc", ".tif", ".tiff", ".h5", ".zip", ".pkl", ".safe"}

# A light product larger than this is treated as heavy. Guards against a CSV
# that turns out to be 400 MB.
MAX_COPY_BYTES = 25 * 1024 * 1024

# Packages worth pinning down, because a change in any of them can move a number.
TRACKED_PACKAGES = ("numpy", "pandas", "xarray", "scipy", "sklearn",
                    "rioxarray", "rasterio", "shapely", "h5py")


def drive_root(root: str | Path | None = None) -> Path:
    """Resolve the archive root.

    Order: explicit argument, then ``INSAR_DRIVE_ROOT``, then the Colab mount
    point. Falls back to ``./outputs/_runs`` when Drive is not mounted, so a
    local run archives *somewhere* rather than failing — the manifest records
    which root was used."""
    if root is not None:
        return Path(root)
    env = os.environ.get(ENV_VAR)
    if env:
        return Path(env)
    colab = Path(DEFAULT_DRIVE_ROOT)
    if colab.is_dir():
        return colab
    return Path("outputs") / "_runs"


def git_state(repo: str | Path = ".") -> dict:
    """Commit SHA and whether the tree was dirty at run time.

    A dirty tree means the archived run cannot be reproduced from the commit
    alone; recording it is more useful than pretending otherwise."""
    def _run(*args: str) -> str | None:
        try:
            out = subprocess.run(args, cwd=str(repo), capture_output=True,
                                 text=True, timeout=10)
            return out.stdout.strip() if out.returncode == 0 else None
        except (OSError, subprocess.SubprocessError):
            return None

    sha = _run("git", "rev-parse", "HEAD")
    status = _run("git", "status", "--porcelain")
    return {
        "commit": sha,
        "short": sha[:7] if sha else None,
        "branch": _run("git", "rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(status) if status is not None else None,
    }


def environment_snapshot() -> dict:
    """Python and package versions, captured now rather than assumed later."""
    versions = {}
    for name in TRACKED_PACKAGES:
        try:
            module = __import__(name)
            versions[name] = getattr(module, "__version__", "unknown")
        except ImportError:
            continue
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": versions,
    }


def run_id(when: datetime | None = None, sha: str | None = None) -> str:
    """A sortable, unique run identifier: ``20260905T171200Z_b4995f9``.

    Timestamp first so a directory listing is chronological."""
    when = when or datetime.now(timezone.utc)
    stamp = when.strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}_{sha}" if sha else stamp


def _digest(path: Path, limit: int = 8 * 1024 * 1024) -> str | None:
    """SHA-256 of a file's first `limit` bytes. Enough to detect substitution
    without reading gigabytes."""
    try:
        h = hashlib.sha256()
        with path.open("rb") as fh:
            h.update(fh.read(limit))
        return h.hexdigest()[:16]
    except OSError:
        return None


def describe_file(path: str | Path) -> dict:
    """Identify a product without necessarily copying it."""
    path = Path(path)
    if not path.exists():
        return {"path": str(path), "exists": False}
    stat = path.stat()
    return {
        "path": str(path),
        "exists": True,
        "bytes": stat.st_size,
        "modified_utc": datetime.fromtimestamp(
            stat.st_mtime, timezone.utc).isoformat(timespec="seconds"),
        "sha256_head": _digest(path),
        "heavy": path.suffix.lower() in HEAVY_SUFFIXES
                 or stat.st_size > MAX_COPY_BYTES,
    }


def archive_run(phase: str,
                outdir: str | Path,
                params: dict | None = None,
                products: dict | None = None,
                root: str | Path | None = None,
                repo: str | Path = ".",
                copy_light: bool = True) -> Path:
    """Archive one phase execution under ``<root>/runs/<phase>/<run_id>/``.

    Parameters
    ----------
    phase
        Phase label, e.g. ``"phaseG"``. Becomes a directory.
    outdir
        The phase's working output directory. Light files in it are copied.
    params
        The inputs that define the run — thresholds, date ranges, zone
        definitions. Whatever you would need to justify the result.
    products
        Named outputs, ``{"series": "outputs/phaseG/series_AC.csv", ...}``.
        Heavy files are recorded by path and hash instead of copied.
    copy_light
        Set False to record only, copying nothing.

    Returns the run directory. Runs are never overwritten: a second call in the
    same second with the same commit raises rather than clobbering.
    """
    outdir = Path(outdir)
    git = git_state(repo)
    rid = run_id(sha=git.get("short"))
    base = drive_root(root) / "runs" / phase / rid
    if base.exists():
        raise FileExistsError(f"run directory already exists: {base}")
    base.mkdir(parents=True)

    described = {name: describe_file(p) for name, p in (products or {}).items()}

    copied = []
    if copy_light and outdir.is_dir():
        for src in sorted(outdir.rglob("*")):
            if not src.is_file():
                continue
            info = describe_file(src)
            if info.get("heavy"):
                continue
            dest = base / "products" / src.relative_to(outdir)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            copied.append(str(src.relative_to(outdir)))

    manifest = {
        "phase": phase,
        "run_id": rid,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": git,
        "environment": environment_snapshot(),
        "parameters": params or {},
        "products": described,
        "copied": copied,
        "outdir": str(outdir),
        "archive_root": str(drive_root(root)),
    }
    (base / "manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8")

    # Keep the existing per-phase manifest working for callers that expect it.
    try:
        from .utils.manifest import write_manifest
        write_manifest(phase, outdir, params or {}, described)
    except Exception:  # pragma: no cover - never fail a run over bookkeeping
        pass

    return base


def list_runs(phase: str, root: str | Path | None = None) -> list[Path]:
    """Every archived run of a phase, oldest first (run ids sort naturally)."""
    d = drive_root(root) / "runs" / phase
    return sorted(p for p in d.glob("*") if p.is_dir()) if d.is_dir() else []


def latest_run(phase: str, root: str | Path | None = None) -> Path | None:
    """The most recent archived run, or None if the phase has never run."""
    runs = list_runs(phase, root)
    return runs[-1] if runs else None


def load_manifest(run_dir: str | Path) -> dict:
    """Read one run's manifest."""
    return json.loads((Path(run_dir) / "manifest.json").read_text(encoding="utf-8"))


def compare_runs(phase: str, root: str | Path | None = None) -> dict:
    """What changed between the two most recent runs of a phase.

    Answers the question that went unanswered when the July numbers moved: did
    the code change, the environment, or the parameters?"""
    runs = list_runs(phase, root)
    if len(runs) < 2:
        return {"comparable": False, "n_runs": len(runs)}
    prev, curr = load_manifest(runs[-2]), load_manifest(runs[-1])
    changed_params = {k: (prev["parameters"].get(k), curr["parameters"].get(k))
                      for k in set(prev["parameters"]) | set(curr["parameters"])
                      if prev["parameters"].get(k) != curr["parameters"].get(k)}
    prev_pkg = prev["environment"]["packages"]
    curr_pkg = curr["environment"]["packages"]
    changed_pkgs = {k: (prev_pkg.get(k), curr_pkg.get(k))
                    for k in set(prev_pkg) | set(curr_pkg)
                    if prev_pkg.get(k) != curr_pkg.get(k)}
    return {
        "comparable": True,
        "previous": prev["run_id"],
        "current": curr["run_id"],
        "commit_changed": prev["git"].get("commit") != curr["git"].get("commit"),
        "parameters_changed": changed_params,
        "packages_changed": changed_pkgs,
    }
