# Colab — how heavy entries run

Only entries that cannot run on the Mac (2 cores, 8 GB) run here: SLC coregistration, phase linking on SLCs,
MintPy-class inversions, null realisations. Everything else is `make run E=<entry>` on the Mac.

## The rule

- **The logic lives in a `.py`**, in `experiments/` (it answers a ticket) or `pipelines/` (it makes data), declared
  in `config/register.yaml` with `runs_on: colab`.
- **The notebook only launches it.** `colab/<entry>.ipynb` is generated, never edited by hand:

  ```bash
  PYTHONPATH=src python -m insar_wetlands.colab <entry> [--ref <branch>] [--extra acquisition] [--pip "<heavy tools>"]
  ```

## Running one

1. Commit and push the entry's code. The launcher checks out `REF` from GitHub, not your working copy.
2. Open `colab/<entry>.ipynb` from GitHub in Colab (*File → Open notebook → GitHub*, repository `AimenSayoud/SAr-adapted`,
   branch `REF`), then *Runtime → Run all*.
3. The run record is written to Drive `insar_rzecin/runs/<entry>/<UTC>_<sha>/manifest.json`. It holds the commit and
   whether the tree was dirty, the package versions, the duration, the exit code and every output file written.

## Where the data are on Colab

`insar_wetlands.paths` resolves them; nothing is typed in a notebook.

| Data | Colab path | Access |
|---|---|---|
| HyP3 crops 2022–2024, phase products, ERA5, `runs/` | `/content/drive/MyDrive/insar_rzecin/` | read and write (the working area) |
| The cube | `…/insar_rzecin/mac_backup/06_data/cube/` | **read only** |
| HyP3 crops 2017–2021, 2025, 2026, annual pairs | `…/insar_rzecin/mac_backup/05_code/local/s1_<period>/` | **read only** |

`mac_backup/` is the Mac's backup, refreshed from the Mac by `05_code/local/scripts/backup_mac_to_drive.py`. A file
written there would not exist on the Mac and would be overwritten by the next backup. Write new products into the
working area, under a folder named after the entry.

## Secrets

Use Colab's key manager (*Secrets*): `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD` or `EARTHDATA_TOKEN`, and `CDSAPI_KEY`.
Never put a secret in a notebook, a commit or a run record.

## When it breaks

`TROUBLESHOOTING.md` lists ten failure modes seen on this project: Drive-mount deadlocks, session locks, git
identity, secrets in headless runs and others, each with its fix.

## Terminal route (headless, from the Mac) — `insar_wetlands.colab_cli`

The same run without a browser, through the `colab` CLI (google-colab-cli; `colab skill` prints its guide). An agent
can do every step except mounting Drive, which Colab reserves for a person:

```bash
python -m insar_wetlands.colab_cli new slc_cslc_x101 -s x101           # allocate a CPU VM
colab drivemount -s x101 /content/drive                                # THE PERSON: open the URL, approve, Enter
python -m insar_wetlands.colab_cli bootstrap slc_cslc_x101 -s x101 --conda "compass isce3 s1reader burst2safe sentineleof sardem gdal libgdal-netcdf"
python -m insar_wetlands.colab_cli netrc -s x101                       # Earthdata line of ~/.netrc → VM, as a file
python -m insar_wetlands.colab_cli run slc_cslc_x101 -s x101 --timeout 21600
python -m insar_wetlands.colab_cli stop -s x101                        # always: an idle VM burns compute units
```

Credentials never appear in executed code (the CLI keeps a history of it); conda tools go into a micromamba
environment on the VM instead of condacolab (whose kernel restart breaks a headless session).
