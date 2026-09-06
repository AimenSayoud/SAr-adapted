"""One entry point for every phase notebook.

Why this module exists. Measured across the 38 notebooks, **108 pairs shared more
than 35 % of their substantive lines**, and the overlap was almost entirely
setup: mount Drive, pull the repo, load the config, resolve the AOI, open the
coherence template, rebuild the zones. Thirty-eight copies of thirty lines is
about eleven hundred lines of preamble, and every copy is an independent chance
for a path or a threshold to drift out of step with the others.

A notebook now opens with two lines::

    from insar_wetlands.bootstrap import start
    ctx = start("phaseG")

and ends with one::

    ctx.archive(params=PARAMS, products=PRODUCTS)

Everything expensive is lazy. Asking for ``ctx.zones`` loads the water mask,
WorldCover and the Sentinel-2 features; not asking costs nothing. That matters
because an acquisition phase needs none of it and should not pay for it.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Any

from .paths import Paths, make_paths

LOG_FORMAT = "%(asctime)s  %(levelname)-7s %(name)s: %(message)s"


def _in_colab() -> bool:
    return "google.colab" in sys.modules or Path("/content").is_dir()


def _configure_logging(phase: str, log_file: Path) -> logging.Logger:
    """Log to the notebook and to ``outputs/<phase>/run.log``.

    The file is what `archive_run` captures, so a run that behaved oddly can be
    read back months later instead of being reconstructed from memory."""
    logger = logging.getLogger(f"insar.{phase}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(logging.Formatter("%(levelname)-7s %(message)s"))
    logger.addHandler(stream)

    try:
        fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        fh.setFormatter(logging.Formatter(LOG_FORMAT))
        logger.addHandler(fh)
    except OSError:          # read-only filesystem in a test; the stream is enough
        pass

    return logger


@dataclass
class Context:
    """Everything a phase needs, resolved once and shared.

    Attributes are lazy: the first access loads, later accesses reuse."""

    phase: str
    cfg: dict
    paths: Paths
    log: logging.Logger
    _extra: dict[str, Any] = field(default_factory=dict)

    # --- inputs -------------------------------------------------------------
    @cached_property
    def pairs(self) -> list[str]:
        """Interferogram pair identifiers found under ``hyp3_cropped/``."""
        from .stack import list_pairs
        pairs = list_pairs(self.paths.cropped)
        self.log.info("pairs: %d", len(pairs))
        return pairs

    @cached_property
    def template(self):
        """One coherence layer, used as the reference grid for everything else."""
        from .stack import load_layer
        t = load_layer(self.paths.cropped, "corr", [self.pairs[0]]).isel(pair=0)
        self.log.info("template grid: %s", dict(t.sizes))
        return t

    @cached_property
    def dem(self):
        from .stack import load_static_layer
        return load_static_layer(self.paths.cropped, "dem")

    @cached_property
    def flooded_fraction(self):
        import xarray as xr

        from .masking.water_mask import flooded_fraction as ff
        path = self.paths.drive_file("water_mask.nc")
        return self.to_grid(ff(xr.open_dataset(path)))

    @cached_property
    def zones(self) -> dict:
        """Zone stratification A–D on the radar grid.

        WorldCover and the Sentinel-2 features are optional refinements: if
        either is unavailable the zones are still defined, and the log says so
        rather than the notebook failing halfway through."""
        import xarray as xr

        from .stratify import define_zones, load_worldcover, s2_landcover_features

        wc = None
        try:
            wc = load_worldcover(self.template, self.cfg, cache_dir=self.paths.cache)
        except Exception as e:                       # noqa: BLE001 - reported, not swallowed
            self.log.warning("WorldCover unavailable, continuing without it: %s", e)

        s2f = None
        try:
            s2f = s2_landcover_features(
                xr.load_dataset(self.paths.drive_file("s2_stack.nc")))
        except Exception as e:                       # noqa: BLE001
            self.log.warning("Sentinel-2 features unavailable: %s", e)

        z = define_zones(self.template, self.cfg, self.flooded_fraction,
                         worldcover=wc, s2_feat=s2f, dem=self.dem)
        counts = {k: z["info"].get(f"n_{k}") for k in "ABCD"}
        self.log.info("zones: %s", counts)
        return z

    # --- helpers ------------------------------------------------------------
    def layer(self, name: str, pairs: list[str] | None = None):
        """Load an interferometric layer onto the template grid."""
        from .stack import load_layer
        return load_layer(self.paths.cropped, name, pairs or self.pairs)

    def to_grid(self, obj):
        """Align anything to the template grid, reprojecting when necessary.

        The plain `align_grid` fails when the CRS is not written on the object;
        the reproject_match fallback handles that case."""
        from .stack import to_grid
        return to_grid(obj, self.template)

    def cache_df(self, tag: str, fn, index_col=None, **kw):
        """Cache a DataFrame under <drive>/figures_cache/."""
        cache_dir = self.paths.drive / "figures_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_df(cache_dir, tag, fn, index_col=index_col, **kw)

    @property
    def outdir(self) -> Path:
        return self.paths.outputs

    def archive(self, params: dict | None = None,
                products: dict | None = None) -> Path:
        """Archive this execution to the Drive. Call at the end of every phase."""
        from .run_archive import archive_run
        run = archive_run(self.phase, self.outdir, params=params,
                          products=products, root=self.paths.drive,
                          repo=self.paths.repo)
        self.log.info("archived run: %s", run)
        return run


def cache_df(cache_dir: str | Path, tag: str, fn, index_col=None, **kw):
    """Cache a DataFrame to CSV in cache_dir.

    `index_col` matters: a frame keyed by its index (like pair_hydro_change
    keyed by `pair`) loses that key if written with index=False, and the next
    run silently merges against a RangeIndex, matching nothing.
    """
    import pandas as pd

    f = Path(cache_dir) / f"{tag}.csv"
    if f.exists():
        return pd.read_csv(f, index_col=index_col, **kw)
    df = fn()
    df.to_csv(f, index=index_col is not None)
    return df


def start(phase: str,
          *,
          mount: bool = True,
          git: bool = True,
          config_path: str | Path | None = None,
          drive_root: str | Path | None = None,
          repo: str | Path | None = None) -> Context:
    """Prepare a phase and return its `Context`.

    Parameters
    ----------
    phase
        Label such as ``"phaseG"``. Names the output directory, the logger and
        the archived run.
    mount, git
        Colab-only steps. Both are skipped silently off Colab, so the same
        notebook runs locally.

    Nothing here is expensive: stacks load on first use.
    """
    from .config import load_config
    from .paths import repo_root

    if mount and _in_colab() and not Path("/content/drive/MyDrive").is_dir():
        try:
            from google.colab import drive as _drive  # type: ignore
            _drive.mount("/content/drive")
        except Exception as e:                          # noqa: BLE001
            logging.getLogger("insar").warning("Drive mount skipped: %s", e)

    # Resolve the repository BEFORE loading the config. `config.load_config`
    # searches upward from the working directory, so passing `repo` would
    # otherwise be silently ignored and the wrong config.yaml read — which is
    # exactly what a test caught.
    root = Path(repo) if repo else repo_root()
    cfg = load_config(config_path or root / "config" / "config.yaml")
    paths = make_paths(phase, cfg=cfg, root=drive_root, repo=root)
    log = _configure_logging(phase, paths.log_file)

    log.info("phase %s starting", phase)
    for k, v in paths.describe().items():
        log.info("  %-16s %s", k, v)

    if not paths.drive.is_dir():
        log.warning("data root does not exist: %s — inputs will not be found",
                    paths.drive)

    if git and _in_colab():
        try:
            from .config import setup_git
            setup_git(paths.repo)
        except Exception as e:                          # noqa: BLE001
            log.warning("git setup skipped: %s", e)

    return Context(phase=phase, cfg=cfg, paths=paths, log=log)
