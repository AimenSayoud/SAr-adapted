# Google Colab Execution Guide

Compute-intensive entries (OPERA COMPASS SLC coregistration, full-stack SLC phase linking, MintPy-class network inversions, and large Monte Carlo permutation nulls) can be executed on Google Colab using thin generated launcher notebooks. All other entries run locally via `make run E=<entry>`.

## Design Rules

- **All scientific logic lives in `.py` modules** under `src/insar_wetlands/`, `experiments/`, or `pipelines/`, declared in `config/register.yaml` with `runs_on: colab`.
- **Launcher notebooks contain no inline science code.** `colab/<entry>.ipynb` is generated automatically and never edited by hand:

  ```bash
  PYTHONPATH=src python -m insar_wetlands.colab <entry> [--ref main] [--extra acquisition] [--pip "<heavy tools>"]
  ```

## Running in the Browser

1. Commit and push the target entry's `.py` script (the launcher checks out `REF` from GitHub).
2. Open `colab/<entry>.ipynb` from GitHub in Google Colab (*File → Open notebook → GitHub*, repository `AimenSayoud/SAr-adapted`, branch `main`), then select *Runtime → Run all*.
3. Execution provenance (`manifest.json`: git SHA, package versions, runtime duration, exit status, and output checksums) is written to `insar_rzecin/runs/<entry>/<UTC>_<sha>/manifest.json` on Google Drive.

## Data Paths on Colab

`insar_wetlands.paths` resolves all storage locations automatically:

| Dataset | Colab Path | Access |
|---|---|---|
| HyP3 crops (2022–2024), phase products, ERA5, `runs/` | `/content/drive/MyDrive/insar_rzecin/` | Read / write working area |
| Harmonised data cube (`raw / silver / gold`) | `/content/drive/MyDrive/insar_rzecin/mac_backup/06_data/cube/` | **Read-only** mirror |
| HyP3 crops (2017–2021, 2025, 2026, annual pairs) | `/content/drive/MyDrive/insar_rzecin/mac_backup/05_code/local/s1_<period>/` | **Read-only** mirror |

## Credentials

Store credentials in Colab's *Secrets* panel (`EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` or `EARTHDATA_TOKEN`, and `CDSAPI_KEY`). Never commit credentials to notebooks, scripts, or run manifests.

## Headless CLI Route (`insar_wetlands.colab_cli`)

Entries can also be executed from a local terminal via `google-colab-cli`:

```bash
python -m insar_wetlands.colab_cli new slc_cslc_x101 -s x101
colab drivemount -s x101 /content/drive                                # Interactive OAuth approval in browser
python -m insar_wetlands.colab_cli bootstrap slc_cslc_x101 -s x101 --conda "compass isce3 s1reader burst2safe sentineleof sardem gdal libgdal-netcdf"
python -m insar_wetlands.colab_cli netrc -s x101                       # Copy Earthdata entry from ~/.netrc to VM
python -m insar_wetlands.colab_cli run slc_cslc_x101 -s x101 --timeout 21600
python -m insar_wetlands.colab_cli stop -s x101                        # Release cloud VM when complete
```

Operational notes (`TROUBLESHOOTING.md`):
- Conda binaries (`COMPASS`, `ISCE3`, `GDAL`) are installed into an isolated `micromamba` prefix rather than `condacolab` so headless kernel sessions are not dropped by a runtime restart.
- `colab_cli run` exports `PROJ_DATA`, `GDAL_DATA`, `PROJ_NETWORK=ON`, and `MPLBACKEND=Agg` automatically.
