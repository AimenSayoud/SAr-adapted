"""Check proposed public files for credential names, local paths and notebook outputs.

Checks tracked files plus non-ignored additions. It does not certify that the Git
history, ignored files or arbitrary binary content are safe to publish.
"""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIVATE_DIRS = {"06_data", "08_deliverables", "runs"}
TOKEN = re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{25,}|github_pat_[A-Za-z0-9_]{30,}|"
                   rb"AKIA[0-9A-Z]{16}|sk-[A-Za-z0-9]{30,}|"
                   rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----)")
LOCAL_PATH = re.compile(rb"/(?:Users|home)/[A-Za-z0-9_.-]+/")


def problems(root: Path = ROOT) -> list[str]:
    names = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=root,
    ).decode().split("\0")
    errors = []
    for name in sorted(set(names) - {""}):
        path = root / name
        if not path.is_file():
            continue
        basename = path.name
        credential = (basename in {".env", ".netrc", ".cdsapirc"}
                      or basename.startswith(".env.") or basename.endswith((".env", ".pem", ".key")))
        if credential and basename != ".env.example":
            errors.append(f"{name}: credential filename")
        if Path(name).parts[0] in PRIVATE_DIRS:
            errors.append(f"{name}: private runtime directory")
        data = path.read_bytes()
        if b"\0" in data[:8192]:
            continue
        if TOKEN.search(data):
            errors.append(f"{name}: credential-shaped literal (value redacted)")
        if LOCAL_PATH.search(data):
            errors.append(f"{name}: personal filesystem path")
        if path.suffix == ".ipynb":
            try:
                notebook = json.loads(data)
                if any(cell.get("outputs") for cell in notebook.get("cells", [])):
                    errors.append(f"{name}: saved notebook outputs")
            except (ValueError, UnicodeDecodeError):
                errors.append(f"{name}: invalid notebook JSON")
    return errors


if __name__ == "__main__":
    found = problems()
    for problem in found:
        print(f"[FAIL] {problem}")
    if not found:
        print("[OK] proposed public files: boundary checks passed")
    raise SystemExit(bool(found))
