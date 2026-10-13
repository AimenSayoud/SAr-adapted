"""Run a register entry on Colab from this terminal with the ``colab`` CLI (google-colab-cli), headless.

The browser launchers (``colab/<entry>.ipynb``) remain for runs started by hand. This is the same run driven from the
Mac, step by step, so an agent can do everything except the one step Colab reserves for a person — mounting Drive:

    python -m insar_wetlands.colab_cli new <entry> -s <session>        # allocate a CPU VM; prints the mount command
    colab drivemount -s <session> /content/drive                        # THE PERSON: open the URL, approve, Enter
    python -m insar_wetlands.colab_cli bootstrap <entry> -s <session> [--conda "pkgs"]   # clone, install
    python -m insar_wetlands.colab_cli netrc -s <session> [--machine urs.earthdata.nasa.gov]  # credentials file
    python -m insar_wetlands.colab_cli run <entry> -s <session> [--timeout 21600] [-- args]
    python -m insar_wetlands.colab_cli stop -s <session>

Rules from ``colab/TROUBLESHOOTING.md`` and the CLI's own guide (``colab skill``):
- ``drive.mount()`` deadlocks headless (issue 2): Drive is mounted by ``colab drivemount``, by the person.
- Colab Secrets cannot be read headless (issue 8), and every executed code string lands in the CLI's history logs:
  credentials go up as a file with ``colab upload`` (only the needed ``~/.netrc`` line, from a private temp file
  deleted afterwards); no secret is ever part of executed code.
- condacolab restarts the kernel, which breaks a headless session: conda tools go into a separate micromamba
  environment (``/content/env``) and the entry runs with that environment's Python.
- Long runs: one blocking ``colab exec`` with a long ``--timeout`` keeps the kernel busy; output streams back.
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_URL = "https://github.com/AimenSayoud/SAr-adapted.git"
VM_REPO = "/content/SAr-adapted"
VM_ENV = "/content/env"
MAMBA = "/content/mm/bin/micromamba"
CONDA_BASE = ("python=3.11 xarray netcdf4 h5py rioxarray rasterio pyproj shapely pyyaml pandas scipy matplotlib "
              "requests")


def colab(*args: str, stdin: str | None = None, timeout: int | None = None) -> int:
    cmd = ["colab", *args]
    print("$", " ".join(shlex.quote(a) for a in cmd), flush=True)
    return subprocess.run(cmd, input=stdin, text=True, timeout=timeout).returncode


def vm_python(code: str, session: str, timeout: int = 3600) -> int:
    return colab("exec", "-s", session, "--timeout", str(timeout), stdin=code, timeout=timeout + 300)


def cmd_new(entry: str, session: str) -> int:
    rc = colab("new", "-s", session)
    if rc == 0:
        print(f"\n>>> Now, in a terminal, the person runs:\n    colab drivemount -s {session} /content/drive\n"
              "    (open the URL it prints, approve Drive access, come back and press Enter)\n", flush=True)
    return rc


def cmd_bootstrap(entry: str, session: str, ref: str, conda: str) -> int:
    code = f'''
import os, subprocess, sys
assert os.path.isdir("/content/drive/MyDrive/insar_rzecin"), "Drive is not mounted: run colab drivemount first"
def sh(c):
    print("+", c, flush=True); subprocess.run(c, shell=True, check=True)
if not os.path.isdir("{VM_REPO}"):
    sh("git clone -q {REPO_URL} {VM_REPO}")
sh("cd {VM_REPO} && git fetch -q origin && git checkout -q {ref} && git pull -q --ff-only || true")
print(subprocess.run(["git", "-C", "{VM_REPO}", "log", "-1", "--oneline"], capture_output=True, text=True).stdout)
'''
    if conda:
        code += f'''
if not os.path.exists("{MAMBA}"):
    sh("mkdir -p /content/mm && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /content/mm bin/micromamba")
if not os.path.exists("{VM_ENV}/bin/python"):
    sh("{MAMBA} create -q -y -p {VM_ENV} -c conda-forge {CONDA_BASE} {conda}")
sh("{VM_ENV}/bin/pip install -q --no-deps -e {VM_REPO}")
if not os.path.exists("{VM_ENV}/share/proj/us_nga_egm08_25.tif"):   # EGM2008 grid: geoid → ellipsoid heights
    sh("PROJ_DATA={VM_ENV}/share/proj {VM_ENV}/bin/projsync --file us_nga_egm08_25.tif --target-dir {VM_ENV}/share/proj")
sh("cd {VM_REPO} && PATH={VM_ENV}/bin:$PATH {VM_ENV}/bin/python -c 'from insar_wetlands import paths; print(\\"drive\\", paths.resolve_drive()); print(\\"cube\\", paths.cube_root())'")
'''
    else:
        code += f'''
sh("cd {VM_REPO} && bash colab/setup.sh")
'''
    code += 'print(">>> bootstrap done", flush=True)\n'
    return vm_python(code, session, timeout=5400)


def cmd_netrc(session: str, machine: str) -> int:
    """Upload only the ``machine`` line of the local ~/.netrc to the VM's /root/.netrc. Nothing is printed."""
    src = Path.home() / ".netrc"
    lines = [ln for ln in src.read_text().splitlines() if ln.split()[:2] == ["machine", machine]] if src.exists() else []
    if not lines:
        raise SystemExit(f"no 'machine {machine}' line in ~/.netrc")
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "netrc"
        f.touch(mode=0o600)
        f.write_text("\n".join(lines) + "\n")
        rc = colab("upload", "-s", session, str(f), "/content/.netrc_upload")
    if rc == 0:      # uploads land under /content; move it home, private (the code carries no secret)
        rc = vm_python('import os, shutil; d = os.path.expanduser("~/.netrc"); '
                       'shutil.move("/content/.netrc_upload", d); os.chmod(d, 0o600); '
                       'print(">>> netrc in place (0600):", d)', session, 300)
    return rc


def cmd_run(entry: str, session: str, timeout: int, args: list[str]) -> int:
    """Run the entry through the register on the VM (micromamba environment if bootstrap made one), streaming."""
    code = f'''
import os, subprocess
conda = os.path.exists("{VM_ENV}/bin/python")
env = dict(os.environ, PYTHONPATH="{VM_REPO}/src", PYTHONUNBUFFERED="1", MPLBACKEND="Agg",  # kernel's inline backend
           # is unknown to the environment's matplotlib and breaks any import that pulls it in
           **({{"PROJ_DATA": "{VM_ENV}/share/proj", "PROJ_LIB": "{VM_ENV}/share/proj", "GDAL_DATA": "{VM_ENV}/share/gdal",
               "PROJ_NETWORK": "ON"}} if conda else {{}}),   # the env is not activated: point PROJ/GDAL at its data
           PATH=("{VM_ENV}/bin:" if conda else "") + os.environ["PATH"])
py = "{VM_ENV}/bin/python" if conda else "python"
p = subprocess.Popen([py, "-m", "insar_wetlands.register", "run", "{entry}", *{["--", *args] if args else []!r}],
                     cwd="{VM_REPO}", env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in p.stdout:
    print(line, end="", flush=True)
print(">>> exit", p.wait(), flush=True)
'''
    return vm_python(code, session, timeout=timeout)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m insar_wetlands.colab_cli", description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["new", "bootstrap", "netrc", "run", "stop"])
    ap.add_argument("entry", nargs="?", default="")
    ap.add_argument("-s", "--session", required=True)
    ap.add_argument("--ref", default="fusion-x062")
    ap.add_argument("--conda", default="", help="conda-forge packages for a micromamba environment on the VM")
    ap.add_argument("--machine", default="urs.earthdata.nasa.gov")
    ap.add_argument("--timeout", type=int, default=21600)
    argv = list(sys.argv[1:] if argv is None else argv)
    rest = argv[argv.index("--") + 1:] if "--" in argv else []          # arguments for the entry itself
    a = ap.parse_args(argv[:argv.index("--")] if "--" in argv else argv)
    if a.step == "new":
        return cmd_new(a.entry, a.session)
    if a.step == "bootstrap":
        return cmd_bootstrap(a.entry, a.session, a.ref, a.conda)
    if a.step == "netrc":
        return cmd_netrc(a.session, a.machine)
    if a.step == "run":
        return cmd_run(a.entry, a.session, a.timeout, rest)
    return colab("stop", "-s", a.session)


if __name__ == "__main__":
    raise SystemExit(main())
