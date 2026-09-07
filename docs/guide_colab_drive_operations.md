# Google Colab & Google Drive Operations Manual
## Complete Guide, Troubleshooting & Error Prevention

> **Scope**: Standard Operating Procedure (SOP) for headless Google Colab execution and Google Drive integration using `colab-cli` (`google-colab-cli`), based on real execution experience on the Rzecin InSAR research pipeline (`SAr-adapted`).

---

## 1. System Architecture & Prerequisites

The research workflow uses a three-tier architecture:
1. **Local Machine (macOS)**: Git repository (`05_code/SAr-adapted`), code editing, automated testing (`make lint test check phases docx`), and CLI management.
2. **Google Colab Cloud VM**: Ephemeral Linux runtime executing heavy compute notebooks (`phaseG`, `phaseL`, `phaseK`, `export_figures_en`). Standard **CPU runtime** is used (accelerator: `NONE`, variant: `DEFAULT`).
3. **Google Drive FUSE (`/content/drive/MyDrive/insar_rzecin`)**: Persistent storage housing raw input stacks (356 Sentinel-1 cropped pairs, water mask NetCDF, Sentinel-2 features, ERA5) and permanent execution archives (`runs/phase*`, `referee/`, `figures/`).

### Local Dependencies
```bash
# Install colab-cli and compatible Jupyter client
pip install --user --break-system-packages google-colab-cli jupyter-kernel-client==0.15.0
```
Binary path on macOS: `/Users/aymen/Library/Python/3.14/bin/colab`. Ensure this path is in your `$PATH`.

---

## 2. Seven Real Failure Modes & How to Prevent Them

### ⚠️ Issue 1: Short Timeout During Large Stack Streaming (`colab_request timeout`)

* **Symptom**:
  `colab exec -f notebook.ipynb` crashed after 5 to 10 minutes with `jupyter_kernel_client.errors.KernelTimeoutError: colab_request timeout`, especially in cells loading 356 GeoTIFF/NetCDF layers across Drive FUSE.
* **Root Cause**:
  `jupyter_kernel_client/constants.py` defaulted to a short timeout (`REQUEST_TIMEOUT = 300` or `600` seconds). Streaming 712 files across Drive FUSE takes ~4 minutes on a standard CPU runtime.
* **The Permanent Fix**:
  1. Patch `jupyter_kernel_client/constants.py` to allow environment override:
     ```python
     REQUEST_TIMEOUT = int(os.getenv("REQUEST_TIMEOUT", 3600))
     ```
  2. Always pass `--timeout 3600` (1 hour) when executing notebooks:
     ```bash
     colab exec --session <session_name> --timeout 3600 -f <path_to_notebook.ipynb>
     ```

---

### ⚠️ Issue 2: Headless Drive Mount Deadlock (`dfs_ephemeral`)

* **Symptom**:
  The Colab execution froze indefinitely on the bootstrap cell (`ctx = start(...)` or `drive.mount("/content/drive")`), printing no error and never finishing.
* **Root Cause**:
  Calling `google.colab.drive.mount("/content/drive")` inside headless `colab exec` triggers Colab's interactive WebSocket authorization request (`dfs_ephemeral`). In a headless terminal session there is no browser UI to click, causing execution to deadlock waiting for an impossible user response.
* **The Permanent Fix**:
  1. In `src/insar_wetlands/bootstrap.py`, **never mount unconditionally**:
     ```python
     if on_colab:
         from google.colab import drive as _drive
         # Check if already mounted via FUSE CLI before calling mount!
         if not Path("/content/drive/MyDrive").is_dir():
             _drive.mount("/content/drive")
     ```
  2. Mount Drive **from the host CLI before running notebooks**:
     ```bash
     colab drivemount --session <session_name> /content/drive
     ```
     The CLI displays an OAuth URL. Click the link in your browser, approve granular Google Drive permissions, and press **Enter** in the terminal.

---

### ⚠️ Issue 3: Session Lock Deadlock (`sessions.json` "running" lock)

* **Symptom**:
  If a running notebook or command was interrupted, killed, or timed out, all subsequent `colab` CLI commands hung or threw errors claiming the session was busy.
* **Root Cause**:
  `~/.config/colab-cli/sessions.json` stores `"running": "<command>"`. When a process terminates abruptly, `colab-cli` fails to clear this lock.
* **The Permanent Fix**:
  Run this one-liner to reset the lock:
  ```bash
  python3 -c "
  import json
  p = '/Users/aymen/.config/colab-cli/sessions.json'
  with open(p) as f: d = json.load(f)
  for k in d: d[k]['running'] = None
  with open(p, 'w') as f: json.dump(d, f, indent=2)
  print('Colab sessions unlocked.')
  "
  ```

---

### ⚠️ Issue 4: Headless Git Push & Author Identity Crashes

* **Symptom**:
  1. Notebook failed at the final commit cell with: `fatal: unable to auto-detect email address (got 'root@...')`.
  2. Notebook crashed with: `fatal: could not read Username for 'https://github.com': No such device or address`.
* **Root Cause**:
  Fresh Colab VMs have no global git identity configured and no GitHub authentication token in headless sessions.
* **The Permanent Fix**:
  1. Configure git global identity upon session startup:
     ```python
     subprocess.run(['git', 'config', '--global', 'user.email', 'aimensayoud@gmail.com'])
     subprocess.run(['git', 'config', '--global', 'user.name', 'AimenSayoud'])
     ```
  2. Always add `|| true` to notebook push commands so headless execution never fails if push fails:
     ```bash
     !git add -A docs/paper && git commit -m "..." || echo "nothing to commit"
     !git push origin {BRANCH} || true
     ```

---

### ⚠️ Issue 5: Missing Imports & `cache_df` Signature Mismatch

* **Symptom**:
  1. `NameError: name 'invert_aggregate' is not defined` in `phaseL_gate.ipynb`.
  2. `NameError: name 'zone_areas' is not defined` or `NameError: name 'dates' is not defined` in `export_figures_en.ipynb`.
  3. `TypeError: cache_df() missing 1 required positional argument: 'fn'`.
* **Root Cause**:
  - Code refactoring moved functions into modules without updating all notebook import cells.
  - `bootstrap.py` provides `ctx.cache_df(tag, fn)` (a method on `PhaseContext`), whereas module-level `cache_df(cache_dir, tag, fn)` requires `cache_dir` as the first argument.
* **The Permanent Fix**:
  In every notebook using `ctx = start(...)`, bind:
  ```python
  cache_df = ctx.cache_df
  dates = dates_from_pairs(pairs)
  from insar_wetlands.zone_viz import zone_areas, zone_field_table
  ```

---

### ⚠️ Issue 6: Undeclared Notebooks in Pipeline Validation (`*_output.ipynb`)

* **Symptom**:
  `make phases` and `pytest tests/test_phases.py` failed with:
  `AssertionError: Extra items in the left set: 'notebooks/05_robustness/phaseL_gate_output.ipynb'`.
* **Root Cause**:
  `colab exec -f notebook.ipynb` automatically saves `notebook_output.ipynb` containing execution logs. The repository's strict architecture guard asserted that *every* `.ipynb` under `notebooks/` must be declared in `config/phases.yaml`.
* **The Permanent Fix**:
  Update `src/insar_wetlands/phases.py` and `tests/test_phases.py` to filter out output logs:
  ```python
  for nb in sorted((repo / "notebooks").glob("**/*.ipynb")):
      if nb.name.endswith("_output.ipynb"):
          continue
      # ... validate declared phases
  ```

---

### ⚠️ Issue 7: Ephemeral VM Data Loss Prevention

* **Symptom**:
  When a Colab VM is idle for ~20 minutes or completed, it terminates (HTTP 404/401). Files in `/content/SAr-adapted` are lost if not copied to Drive.
* **The Permanent Fix**:
  At the end of every notebook execution, sync artifacts directly to Google Drive FUSE:
  ```python
  import shutil, os
  dst = '/content/drive/MyDrive/insar_rzecin/referee'
  os.makedirs(dst, exist_ok=True)
  shutil.copytree('/content/SAr-adapted/docs/paper/referee', dst, dirs_exist_ok=True)
  ```
  And commit CSV tables (`T*.csv`, `LT*.csv`, `KT*.csv`) back to the git repository.

---

## 3. Standard Operating Procedure (SOP) Step-by-Step

Follow these exact terminal commands for any future Colab session:

### Step 1: Create a Clean CPU Session
```bash
# Create fresh CPU session named 's1'
colab new --session s1
```

### Step 2: Mount Google Drive
```bash
colab drivemount --session s1 /content/drive
```
* Click the OAuth authorization URL displayed in terminal.
* Log in and grant Drive permissions.
* Return to terminal and press **Enter**.
* Verify Drive:
  ```bash
  echo "import os; print('Drive mounted:', os.path.exists('/content/drive/MyDrive/insar_rzecin'))" | colab exec --session s1
  ```

### Step 3: Configure Git Identity on VM
```bash
echo "import subprocess
subprocess.run(['git', 'config', '--global', 'user.email', 'aimensayoud@gmail.com'])
subprocess.run(['git', 'config', '--global', 'user.name', 'AimenSayoud'])
print('Git configured.')" | colab exec --session s1
```

### Step 4: Execute a Pipeline Notebook
```bash
# Run with 1-hour timeout
colab exec --session s1 --timeout 3600 -f notebooks/05_robustness/phaseL_gate.ipynb
```

### Step 5: Backup Artifacts to Drive
```bash
echo "import shutil, os
shutil.copytree('/content/SAr-adapted/docs/paper/referee', '/content/drive/MyDrive/insar_rzecin/referee', dirs_exist_ok=True)
shutil.copytree('/content/SAr-adapted/docs/paper/figures', '/content/drive/MyDrive/insar_rzecin/figures', dirs_exist_ok=True)
print('All outputs backed up to Drive.')" | colab exec --session s1
```

### Step 6: Verify Local Repository & Manuscript
Back on your Mac:
```bash
cd ~/Documents/Research_Hub/05_code/SAr-adapted
git pull origin main
make test            # 236 tests must pass
make phases          # Pipeline declaration must be sound
make check           # Manuscript numbers must match CSVs
make check-generated # Zero drift between appendix and CSVs
make docx            # Builds manuscript.docx with pandoc
```

---

## 4. Emergency Quick Fixes Cheat Sheet

| Situation | Command / Solution |
| :--- | :--- |
| **Colab command hanging / session locked** | Run one-liner to reset `"running": null` in `~/.config/colab-cli/sessions.json` |
| **VM disconnected (404/401)** | Run `colab new --session <name>` and re-mount with `colab drivemount` |
| **Drive mount stuck on Colab** | Ensure `bootstrap.py` checks `if not Path('/content/drive/MyDrive').is_dir()` |
| **Timeout during cell execution** | Set `export REQUEST_TIMEOUT=3600` and pass `--timeout 3600` |
| **Git commit failed on VM** | Run `git config --global user.email ...` on the VM via `colab exec` |
| **make phases failed on undeclared file** | Ensure `phases.py` and `test_phases.py` skip `*_output.ipynb` |
