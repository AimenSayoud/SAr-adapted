# Google Colab & Drive Troubleshooting Reference

Practical solutions for common failure modes encountered when running heavy SAR processing entries (`COMPASS` CSLC coregistration, ISCE3/GDAL stacks, and large NetCDF/GeoTIFF streaming) on Google Colab.

---

## 1. Drive FUSE Streaming Timeouts (`KernelTimeoutError`)
- **Symptom:** Large stack reads across `/content/drive/MyDrive/...` time out after 5–10 minutes in headless CLI runs.
- **Fix:** Pass an explicit `--timeout` (e.g., `--timeout 21600` in `python -m insar_wetlands.colab_cli run ...`) and set `REQUEST_TIMEOUT=3600` in the host environment when streaming hundreds of interferograms over Drive FUSE.

## 2. Headless Drive Mount Deadlock (`dfs_ephemeral`)
- **Symptom:** Calling `google.colab.drive.mount('/content/drive')` inside a non-interactive terminal session hangs indefinitely waiting for browser OAuth confirmation.
- **Fix:** Check `Path('/content/drive/MyDrive').is_dir()` before calling `drive.mount(...)`, and mount Drive interactively beforehand via `colab drivemount -s <session> /content/drive`.

## 3. Stale CLI Session Lock (`sessions.json`)
- **Symptom:** After an interrupted terminal run, subsequent `colab` CLI commands report the session as busy.
- **Fix:** Clear stale `"running"` flags in `~/.config/colab-cli/sessions.json`, or allocate a fresh session and terminate orphaned runtimes via *Runtime → Manage sessions* in the Colab browser UI.

## 4. Headless Credential Access (`userdata.get` Timeout)
- **Symptom:** `google.colab.userdata.get(...)` times out in headless runs because browser confirmation cannot be clicked.
- **Fix:** Either grant secret access interactively once in the browser UI, or inject `~/.netrc` / environment variables via `python -m insar_wetlands.colab_cli netrc -s <session>`. `insar_wetlands.config.get_secret()` falls back automatically to environment variables and `~/.netrc` / `~/.cdsapirc`.

## 5. Avoiding Kernel Restarts with `micromamba`
- **Symptom:** Installing `ISCE3` / `COMPASS` via `condacolab` forces a kernel restart that drops headless CLI connections.
- **Fix:** `insar_wetlands.colab_cli bootstrap` installs binary SAR dependencies into `/root/micromamba/envs/sar` and invokes that environment directly with `PROJ_DATA`, `GDAL_DATA`, and `MPLBACKEND=Agg` set.

## 6. Scoped Python Dependencies (`colab/requirements.txt`)
- **Symptom:** Full-system `pip freeze` lockfiles fail to compile OS packages (`dbus-python`, `python-apt`) on fresh Colab VMs.
- **Fix:** `colab/setup.sh` installs only the scoped scientific dependencies in `colab/requirements.txt` followed by `pip install -e .` (or `pip install --no-deps -e .` when paired with a `micromamba` SAR environment).
