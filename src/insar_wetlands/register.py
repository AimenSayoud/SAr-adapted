"""The register of everything that runs or has run: ``config/register.yaml`` (the public, technical half).

One entry per notebook phase (historical, Colab), experiment (answers a ticket, writes a hub deliverable),
pipeline (makes data) or Colab launcher. It says what each entry is, where its code is (``code`` at ``code_ref``:
``HEAD`` for the working tree, a tag for code removed from it), how to run it, what it reads and writes, and which
tickets and tries it belongs to. What an entry *found* is not here — this repo is public — but in the hub's tries
map and ledger; the hub joins both halves into ``05_code/REGISTER.md``.

    from insar_wetlands.register import load_register, validate, find
    reg = load_register()
    find(reg, ticket="X-082")          # entries that answered or serve a ticket
    find(reg, chain="fusion-v3")       # entries of a work chain
    find(reg, try_id="T13")            # entries a try rests on
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
REGISTER = REPO / "config" / "register.yaml"
KINDS = ("phase", "experiment", "pipeline", "colab")
STATUSES = ("current", "finished", "superseded", "open", "removed")
RUNS_ON = ("mac", "colab")


@dataclass
class Entry:
    id: str
    kind: str
    status: str
    runs_on: str
    code: str | None
    code_ref: str = "HEAD"
    title: str = ""
    question: str = ""
    command: str | None = None
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    tickets: list[str] = field(default_factory=list)
    chains: list[str] = field(default_factory=list)
    tries: list[str] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    @property
    def in_tree(self) -> bool:
        return self.code_ref == "HEAD"


KNOWN = {f for f in Entry.__dataclass_fields__ if f not in ("id", "extra")}


def load_register(path: str | Path | None = None) -> dict[str, Entry]:
    raw = yaml.safe_load(Path(path or REGISTER).read_text(encoding="utf-8")) or {}
    out = {}
    for name, e in (raw.get("entries") or {}).items():
        e = dict(e or {})
        known = {k: e.pop(k) for k in list(e) if k in KNOWN}
        out[name] = Entry(id=name, extra=e, **{"code": None, "kind": "", "status": "", "runs_on": "", **known})
    return out


def _git_has(repo: Path, ref: str, path: str) -> bool | None:
    """Does ``path`` exist at ``ref``? None when git cannot answer (no repo, ref not fetched)."""
    try:
        r = subprocess.run(["git", "-C", str(repo), "cat-file", "-e", f"{ref}:{path}"],
                           capture_output=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return None
    if r.returncode == 0:
        return True
    probe = subprocess.run(["git", "-C", str(repo), "rev-parse", "--verify", "--quiet", ref], capture_output=True)
    return False if probe.returncode == 0 else None


def validate(reg: dict[str, Entry], repo: str | Path | None = None, check_git: bool = True) -> list[str]:
    """Problems with the register, empty when sound."""
    repo = Path(repo or REPO)
    bad = []
    for name, e in reg.items():
        if e.kind not in KINDS:
            bad.append(f"{name}: kind {e.kind!r} not in {KINDS}")
        if e.status not in STATUSES:
            bad.append(f"{name}: status {e.status!r} not in {STATUSES}")
        if e.runs_on not in RUNS_ON:
            bad.append(f"{name}: runs_on {e.runs_on!r} not in {RUNS_ON}")
        if e.status == "open":
            continue
        if not e.code:
            bad.append(f"{name}: no code path")
            continue
        if e.in_tree:
            if not (repo / e.code).exists():
                bad.append(f"{name}: {e.code} is not in the tree (removed? set code_ref to the tag that has it)")
        elif check_git and _git_has(repo, e.code_ref, e.code) is False:
            bad.append(f"{name}: {e.code} not found at {e.code_ref}")
        if e.in_tree and e.status in ("current", "finished") and e.kind != "phase" and not e.command:
            bad.append(f"{name}: runnable entry without a command")
        for t in e.tickets:
            if not (len(t) == 5 and t[1] == "-" and t[0] in "XCDRP" and t[2:].isdigit()):
                bad.append(f"{name}: ticket id {t!r} malformed")
        for t in e.tries:
            if not (t.startswith("T") and t[1:].isdigit()):
                bad.append(f"{name}: try id {t!r} malformed")
    return bad


def untracked_scripts(reg: dict[str, Entry], repo: str | Path | None = None) -> list[str]:
    """Runnable files in the tree that no entry declares (experiments/, pipelines/, colab/)."""
    repo = Path(repo or REPO)
    declared = {e.code for e in reg.values() if e.in_tree}
    files = [p for d in ("experiments", "pipelines") for p in (repo / d).glob("*.py")]
    files += list((repo / "colab").glob("*.ipynb")) if (repo / "colab").exists() else []
    launchers = {f"colab/{n}.ipynb" for n, e in reg.items() if e.runs_on == "colab"}
    return sorted(str(p.relative_to(repo)) for p in files
                  if str(p.relative_to(repo)) not in declared | launchers and not p.name.startswith("_"))


def find(reg: dict[str, Entry], ticket: str | None = None, chain: str | None = None,
         try_id: str | None = None, kind: str | None = None) -> list[Entry]:
    out = []
    for e in reg.values():
        if ticket and ticket not in e.tickets:
            continue
        if chain and chain not in e.chains:
            continue
        if try_id and try_id not in e.tries:
            continue
        if kind and e.kind != kind:
            continue
        out.append(e)
    return out


def producers_consumers(reg: dict[str, Entry]) -> tuple[dict, dict, dict]:
    """(output basename → producing entries, input basename → consuming entries, entry → status) — the provenance
    hints the hub's data catalog reads."""
    prod: dict[str, list[str]] = {}
    cons: dict[str, list[str]] = {}
    for name, e in reg.items():
        for o in e.outputs:
            prod.setdefault(Path(str(o).split(":", 1)[-1].rstrip("/")).name, []).append(name)
        for i in e.inputs:
            ref = str(i)
            head, _, tail = ref.partition(".")
            if tail and head in reg:          # phase-style reference 'phase03.pair_selection.csv'
                ref = tail
            cons.setdefault(Path(ref.split(":", 1)[-1].rstrip("/")).name, []).append(name)
    return prod, cons, {n: e.status for n, e in reg.items()}


def summary(reg: dict[str, Entry]) -> str:
    from collections import Counter
    c = Counter((e.kind, e.status) for e in reg.values())
    return "\n".join(f"  {k:10s} {s:10s} {n}" for (k, s), n in sorted(c.items()))


# ------------------------------------------------------------------------------------------------- command line
def run_entry(name: str, args: list[str] | None = None, reg: dict[str, Entry] | None = None) -> int:
    """Run one register entry on this machine. Colab-only and removed entries are refused with the reason."""
    import os
    import sys
    reg = reg or load_register()
    e = reg.get(name)
    if e is None:
        raise SystemExit(f"no register entry {name!r} (python -m insar_wetlands.register list)")
    from .paths import on_colab
    if e.runs_on == "colab" and not on_colab():
        raise SystemExit(f"{name} runs on Colab: open colab/{name}.ipynb (see colab/README.md)")
    if not e.in_tree:
        raise SystemExit(f"{name} is not in the tree: its code is at {e.code_ref}:{e.code}")
    if e.status == "superseded":
        raise SystemExit(f"{name} is superseded: run what replaced it, not this")
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    cmd = [sys.executable, str(REPO / e.code), *(args or [])]
    print(">>", " ".join(cmd), flush=True)
    import time
    t0 = time.time()
    code = subprocess.call(cmd, env=env, cwd=REPO)
    rec = record_run(e, args or [], code, t0)
    print(f"run record: {rec}", flush=True)
    return code


def runs_root() -> Path:
    """Where run records go: the hub's ``05_code/runs`` on the Mac (private: outputs hold field values), the Drive's
    ``runs/`` on Colab."""
    from .paths import hub_root, on_colab, resolve_drive
    hub = hub_root()
    if hub is not None and not on_colab():
        return hub / "05_code" / "runs"
    return resolve_drive() / "runs"


def record_run(e: Entry, args: list[str], exit_code: int, started: float) -> Path:
    """One manifest per execution: entry, command, commit (and whether the tree was dirty), environment, duration,
    exit code, and every file under the entry's declared output folders written during the run."""
    import json
    import time
    from datetime import datetime, timezone

    from .paths import hub_root
    from .run_archive import describe_file, environment_snapshot, git_state, run_id
    git = git_state(REPO)
    rid = run_id(datetime.fromtimestamp(started, timezone.utc), git.get("short"))
    out_dir = runs_root() / e.id / rid
    out_dir.mkdir(parents=True, exist_ok=True)
    hub = hub_root()
    written = []
    for o in e.outputs:
        root, _, rel = str(o).partition(":")
        base = (hub / rel) if (root == "hub" and hub is not None) else None
        if base is None or not base.exists():
            continue
        files = [base] if base.is_file() else [f for f in base.rglob("*") if f.is_file()]
        written += [describe_file(f) for f in files if f.stat().st_mtime >= started - 1]
    manifest = {"entry": e.id, "kind": e.kind, "tickets": e.tickets, "tries": e.tries, "code": e.code,
                "args": args, "exit_code": exit_code, "started_utc": datetime.fromtimestamp(started, timezone.utc)
                .isoformat(timespec="seconds"), "duration_s": round(time.time() - started, 1), "git": git,
                "environment": environment_snapshot(), "outputs_written": written}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
    return out_dir


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(prog="python -m insar_wetlands.register",
                                 description="The code register: check it, find entries, run one.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="validate the register and find undeclared scripts")
    ls = sub.add_parser("list", help="list entries, optionally filtered")
    ls.add_argument("--ticket"), ls.add_argument("--chain"), ls.add_argument("--try", dest="try_id"), ls.add_argument("--kind")
    sh = sub.add_parser("show", help="one entry in full")
    sh.add_argument("name")
    rn = sub.add_parser("run", help="run one entry on this machine")
    rn.add_argument("name")
    rn.add_argument("args", nargs=argparse.REMAINDER)
    a = ap.parse_args(argv)
    reg = load_register()
    if a.cmd == "check":
        bad = validate(reg) + [f"undeclared: {u}" for u in untracked_scripts(reg)]
        print(summary(reg))
        print("\n".join("  " + b for b in bad) if bad else "  register is sound")
        return 1 if bad else 0
    if a.cmd == "list":
        for e in find(reg, ticket=a.ticket, chain=a.chain, try_id=a.try_id, kind=a.kind):
            print(f"{e.id:30s} {e.kind:10s} {e.status:10s} {e.runs_on:6s} {','.join(e.tickets):18s} {e.title[:70]}")
        return 0
    if a.cmd == "show":
        e = reg.get(a.name)
        if e is None:
            raise SystemExit(f"no register entry {a.name!r}")
        print(yaml.safe_dump({e.id: {k: v for k, v in e.__dict__.items() if v not in (None, [], {}, "") and k != "id"}},
                             sort_keys=False, allow_unicode=True, width=120))
        return 0
    return run_entry(a.name, [x for x in a.args if x != "--"], reg)


if __name__ == "__main__":
    raise SystemExit(main())
