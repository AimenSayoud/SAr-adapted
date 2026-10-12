"""The one place a path is constructed.

Why this module exists. Paths were being built ad hoc in every notebook and in a
few module defaults — including ``cache_dir="/content/worldcover"`` sitting in a
function signature, which is silently wrong anywhere but Colab and, being a
default, does not raise: it writes somewhere unexpected instead.

Thirty-eight notebooks each carrying their own idea of where the Drive is means
thirty-eight places to edit when the layout changes, and one that gets missed.
After this module, no other code contains a path literal.

Resolution order for the data root, most explicit first:

1. the ``root`` argument;
2. the ``INSAR_DRIVE_ROOT`` environment variable — how tests and local runs override;
3. ``paths.drive_data_root`` in ``config/config.yaml``, if that directory exists;
4. the Colab mount point, if it exists;
5. on the Mac, the hub's read-only Drive snapshot ``05_code/local/drive_pristine`` (never the live mount);
6. ``<repo>/data/_local`` as a last resort, so a local run works rather than failing.

Step 3 checks existence deliberately: the committed config names the Colab path,
which is meaningless off Colab, and falling through to a working local directory
is friendlier than a confusing error.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

ENV_VAR = "INSAR_DRIVE_ROOT"
COLAB_DRIVE = Path("/content/drive/MyDrive/insar_rzecin")


def repo_root(start: str | Path | None = None) -> Path:
    """The directory containing ``config/config.yaml``, searching upward."""
    p = Path(start or Path.cwd()).resolve()
    for cand in [p, *p.parents]:
        if (cand / "config" / "config.yaml").exists():
            return cand
    raise FileNotFoundError(
        "config/config.yaml not found — run from inside the cloned repository."
    )


def resolve_drive(root: str | Path | None = None, cfg: dict | None = None,
                  repo: Path | None = None) -> Path:
    """Resolve the data root. See the module docstring for the order."""
    if root is not None:
        return Path(root).expanduser()

    env = os.environ.get(ENV_VAR)
    if env:
        return Path(env).expanduser()

    if cfg:
        configured = (cfg.get("paths") or {}).get("drive_data_root")
        if configured and Path(configured).is_dir():
            return Path(configured)

    if COLAB_DRIVE.is_dir():
        return COLAB_DRIVE

    pristine = _hub_snapshot(repo)
    if pristine is not None:
        return pristine

    return (repo or repo_root()) / "data" / "_local"


def _hub_snapshot(repo: Path | None = None) -> Path | None:
    """On the Mac: the read-only snapshot of the Drive (``<hub>/05_code/local/drive_pristine``), never the live
    CloudStorage mount — a run must not be able to overwrite a Drive product in place."""
    here = Path(repo).resolve() if repo else Path(__file__).resolve().parents[2]
    snap = here.parent / "local" / "drive_pristine" if here.parent.name == "05_code" else None
    return snap if snap is not None and snap.is_dir() else None


# ------------------------------------------------------------------------------------- the hub's data roots
# The cube and the HyP3 crops outside 2022–2024 are made on the Mac (hub ``06_data/cube``, ``05_code/local/s1_*``) and
# copied to Drive ``insar_rzecin/mac_backup/`` (05_code/local/scripts/backup_mac_to_drive.py). On Colab those copies
# are the only way to reach them; they are read there, never written (a write would be overwritten by the next backup
# and would not exist on the Mac). Each root can be forced by an environment variable.
CUBE_ENV, S1_ENV, HUB_ENV = "RZECIN_CUBE_ROOT", "RZECIN_S1_ROOT", "RZECIN_HUB"
MAC_BACKUP = COLAB_DRIVE / "mac_backup"
S1_PERIODS_ON_DRIVE = ("2022_2024",)          # the Colab working area holds only these crops (hyp3_cropped*/)


def on_colab() -> bool:
    return COLAB_DRIVE.is_dir() and Path("/content").is_dir()


def hub_root() -> Path | None:
    """The private hub around this repository (``<hub>/05_code/SAr-adapted``), or ``RZECIN_HUB``; None on Colab."""
    if os.environ.get(HUB_ENV):
        return Path(os.environ[HUB_ENV]).expanduser()
    here = Path(__file__).resolve().parents[2]
    return here.parents[1] if here.parent.name == "05_code" else None


def cube_root() -> Path:
    """The data cube: ``RZECIN_CUBE_ROOT``, else the hub's ``06_data/cube``, else (Colab) its Drive backup copy."""
    if os.environ.get(CUBE_ENV):
        return Path(os.environ[CUBE_ENV]).expanduser()
    hub = hub_root()
    if hub is not None and (hub / "06_data" / "cube").is_dir():
        return hub / "06_data" / "cube"
    return MAC_BACKUP / "06_data" / "cube"


def s1_local_root() -> Path:
    """Where the HyP3 crops outside 2022–2024 are: ``RZECIN_S1_ROOT``, the hub's ``05_code/local``, or (Colab) the
    Drive backup copy of it."""
    if os.environ.get(S1_ENV):
        return Path(os.environ[S1_ENV]).expanduser()
    hub = hub_root()
    if hub is not None and (hub / "05_code" / "local").is_dir():
        return hub / "05_code" / "local"
    return MAC_BACKUP / "05_code" / "local"


def s1_crops(period: str, track: str, drive: Path | None = None) -> Path:
    """The folder of HyP3 crops (one sub-folder per pair) for a period (``2017_2019`` … ``2026``, ``annual``) and track.
    2022–2024 lives in the Drive working area (``hyp3_cropped`` / ``hyp3_cropped_descending``); the rest under
    ``s1_local_root()/s1_<period>/hyp3_cropped_<track>``."""
    if period in S1_PERIODS_ON_DRIVE:
        root = Path(drive) if drive else resolve_drive()
        return root / ("hyp3_cropped" if track == "ascending" else f"hyp3_cropped_{track}")
    return s1_local_root() / f"s1_{period}" / f"hyp3_cropped_{track}"


@dataclass(frozen=True)
class Paths:
    """Every directory the pipeline uses, resolved once.

    Directories are created on access rather than at construction, so building a
    ``Paths`` is free and safe in a test."""

    repo: Path
    drive: Path
    phase: str | None = None
    track: str | None = None

    # --- inputs -------------------------------------------------------------
    @property
    def cropped(self) -> Path:
        """Cropped HyP3 burst interferograms — the pipeline's main input."""
        if self.track and self.track.lower() != "ascending":
            return self.drive / f"hyp3_cropped_{self.track.lower()}"
        return self.drive / "hyp3_cropped"

    @property
    def cache(self) -> Path:
        """Downloaded auxiliary rasters (WorldCover tiles, DEM extracts).

        Replaces the hardcoded ``/content/worldcover``. Lives on the Drive so a
        tile survives the Colab session that fetched it."""
        return self._made(self.drive / "cache")

    @property
    def artifacts(self) -> Path:
        """Small inter-phase artefacts (``pair_selection.csv``, reference-pixel
        JSON), namespaced by track so a non-default-track run cannot overwrite
        the default track's copy."""
        if self.track and self.track.lower() != "ascending":
            return self._made(self.drive / f"artifacts_{self.track.lower()}")
        return self._made(self.drive / "artifacts")

    def track_file(self, name: str) -> Path:
        """A Drive-root product (e.g. ``ts_sbas_ref_only.nc``), suffixed by
        track when non-default so it cannot collide with the default track's
        file of the same name. Identical to ``drive_file(name)`` for the
        default track — no change to the existing ascending pipeline."""
        if self.track and self.track.lower() != "ascending":
            stem = Path(name)
            suffix = stem.suffix
            return self.drive / f"{stem.stem}_{self.track.lower()}{suffix}"
        return self.drive / name

    # --- outputs ------------------------------------------------------------
    @property
    def outputs(self) -> Path:
        """``outputs/<phase>/`` in the repo, or ``outputs/<phase>_<track>/`` when non-default."""
        base = self.repo / "outputs"
        if self.phase:
            suffix = f"_{self.track.lower()}" if (self.track and self.track.lower() != "ascending") else ""
            return self._made(base / f"{self.phase}{suffix}")
        return self._made(base)

    @property
    def runs(self) -> Path:
        """Archived executions, on the Drive so they outlive the session."""
        if self.track and self.track.lower() != "ascending":
            return self._made(self.drive / f"runs_{self.track.lower()}")
        return self._made(self.drive / "runs")

    @property
    def figures(self) -> Path:
        """Exported figures (PNG) of the write-up."""
        return self.repo / "docs" / "paper" / "figures"

    @property
    def tables(self) -> Path:
        """Results tables ``T*.csv`` and the aggregated series — the numbers' source of truth."""
        return self.repo / "results" / "tables"

    @property
    def robustness(self) -> Path:
        """Robustness tables and figures (K*, L*, X*) of the main results."""
        return self.repo / "results" / "robustness"

    @property
    def paper(self) -> Path:
        return self.repo / "docs" / "paper"

    @property
    def log_file(self) -> Path:
        return self.outputs / "run.log"

    # --- helpers ------------------------------------------------------------
    @cached_property
    def on_colab(self) -> bool:
        return self.drive == COLAB_DRIVE

    def drive_file(self, name: str) -> Path:
        """A named artefact at the Drive root, e.g. ``water_mask.nc``."""
        return self.drive / name

    def for_phase(self, phase: str, track: str | None = None) -> Paths:
        return Paths(repo=self.repo, drive=self.drive, phase=phase,
                     track=track or self.track)

    @staticmethod
    def _made(p: Path) -> Path:
        p.mkdir(parents=True, exist_ok=True)
        return p

    def describe(self) -> dict:
        """What resolved to what — logged at the start of every run, so a
        surprising result can be traced to a surprising path."""
        return {
            "repo": str(self.repo),
            "drive": str(self.drive),
            "phase": self.phase,
            "track": self.track or "ascending",
            "on_colab": self.on_colab,
            "drive_exists": self.drive.is_dir(),
            "cropped_exists": self.cropped.is_dir(),
        }


def make_paths(phase: str | None = None, cfg: dict | None = None,
               root: str | Path | None = None,
               repo: str | Path | None = None,
               track: str | None = None) -> Paths:
    """Build the `Paths` for a phase and optional track."""
    r = Path(repo) if repo else repo_root()
    return Paths(repo=r, drive=resolve_drive(root, cfg, r), phase=phase, track=track)
