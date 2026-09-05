"""Read and validate the declared pipeline in ``config/phases.yaml``.

Why this module exists. The execution order, what each phase consumes and
produces, and which of the ``b`` / ``bis`` / ``ter`` variants is still current
lived only in filenames and in README prose that covered 14 of 38 notebooks. A
filename cannot say *what changed*, and prose cannot be checked.

With the pipeline declared, four things stop being maintained by hand:

* the README phase table is generated, so it cannot fall out of date again;
* ``status: superseded`` records which notebook is dead, in the place a reader
  looks rather than in a suffix;
* ``paper:`` maps each phase to the claim it supports — the reverse direction of
  ``traceability.md``;
* a phase whose inputs are missing is reported before a Colab session is spent
  discovering it halfway through.

The graph is deliberately validated rather than trusted: a typo in an ``inputs``
entry is exactly the sort of thing that would otherwise be discovered by a
confusing failure three phases later.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

STATUSES = {"current", "superseded", "exploratory", "tooling"}

# Groups answer "what question is this phase part of", which is the axis the
# filenames never carried: two naming schemes grew side by side (phase01..15 by
# chronology, phaseA..L by hypothesis) and neither says that phaseE2 is the
# decisive H1 test while phase08 and phase09 are the inversions it supersedes.
#
# A group is orthogonal to `status`: phaseC1 is exploratory but belongs with the
# inversions, and phaseC2 is superseded but is filed next to the phaseE2 that
# replaced it — which is where a reader looks for it.
#
# The value is also the directory under `notebooks/`, so the declaration and the
# filesystem cannot disagree: `validate` checks that every notebook actually
# lives in its group's directory. The numeric prefixes give reading order to a
# listing that was otherwise 38 files in no order at all.
GROUPS = {
    "data":        "01_data",
    "inversion":   "02_inversion",
    "corrections": "03_corrections",
    "hypotheses":  "04_hypotheses",
    "robustness":  "05_robustness",
    "manuscript":  "06_manuscript",
}


@dataclass(frozen=True)
class Phase:
    """One declared phase."""

    name: str
    title: str
    notebook: str
    # Defaulted rather than required so the synthetic graphs in the tests stay
    # readable; `validate` reports an entry that omits it, which is what keeps
    # the real declaration honest.
    group: str = ""
    question: str = "-"
    status: str = "current"
    supersedes: list[str] = field(default_factory=list)
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    paper: str = "-"

    @property
    def is_live(self) -> bool:
        """Should this phase be re-run? Superseded phases are kept for
        provenance and must not be."""
        return self.status in {"current", "tooling"}

    def upstream(self) -> list[str]:
        """Phase names this one depends on.

        An input written ``phaseD.coh_by_zone.csv`` names another phase; one
        written ``water_mask.nc`` is a bare artefact with no declared producer."""
        out = []
        for i in self.inputs:
            head = i.split(".")[0]
            if head.startswith("phase") or head in {"export_figures_en",
                                                    "build_manuscript_docx"}:
                out.append(head)
        return out


def load_phases(path: str | Path | None = None) -> dict[str, Phase]:
    """Parse ``config/phases.yaml`` into `Phase` objects."""
    if path is None:
        from .paths import repo_root
        path = repo_root() / "config" / "phases.yaml"
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    return {name: Phase(name=name, **body) for name, body in raw.items()}


def validate(phases: dict[str, Phase], repo: str | Path | None = None) -> list[str]:
    """Every problem found, as readable lines. Empty means the graph is sound.

    Returns problems rather than raising, so one run reports all of them
    instead of stopping at the first."""
    problems: list[str] = []
    repo = Path(repo) if repo else None

    for name, p in phases.items():
        if p.status not in STATUSES:
            problems.append(f"{name}: unknown status {p.status!r} "
                            f"(expected one of {sorted(STATUSES)})")

        if repo and not (repo / p.notebook).exists():
            problems.append(f"{name}: notebook not found — {p.notebook}")

        # The group and the directory are two statements of the same fact, so a
        # notebook filed in the wrong place is a declaration error, not a
        # cosmetic one: it is how the flat layout would grow back.
        if not p.group:
            problems.append(f"{name}: no group declared "
                            f"(expected one of {sorted(GROUPS)})")
        elif p.group not in GROUPS:
            problems.append(f"{name}: unknown group {p.group!r} "
                            f"(expected one of {sorted(GROUPS)})")
        else:
            want = f"notebooks/{GROUPS[p.group]}"
            if str(Path(p.notebook).parent) != want:
                problems.append(
                    f"{name}: group {p.group!r} means the notebook belongs in "
                    f"{want}/, but it is declared at {p.notebook}")

        for dep in p.upstream():
            if dep not in phases:
                problems.append(f"{name}: input refers to unknown phase {dep!r}")
            elif not phases[dep].is_live and p.is_live:
                problems.append(
                    f"{name} is {p.status} but depends on {dep}, which is "
                    f"{phases[dep].status} — a live phase must not read from a dead one")

        for s in p.supersedes:
            if s not in phases:
                problems.append(f"{name}: supersedes unknown phase {s!r}")
            elif phases[s].status != "superseded":
                problems.append(
                    f"{name} claims to supersede {s}, but {s} is still "
                    f"marked {phases[s].status!r}")

    # A notebook on disk that nobody declared is how the README fell behind.
    if repo:
        declared = {p.notebook for p in phases.values()}
        for nb in sorted((repo / "notebooks").glob("**/*.ipynb")):
            rel = str(nb.relative_to(repo))
            if rel not in declared:
                problems.append(f"undeclared notebook: {rel}")

    problems.extend(f"dependency cycle: {' -> '.join(c)}" for c in find_cycles(phases))
    return problems


def find_cycles(phases: dict[str, Phase]) -> list[list[str]]:
    """Cycles in the dependency graph, each as the path that closes it."""
    cycles, state = [], {}

    def walk(name: str, trail: list[str]) -> None:
        if state.get(name) == "done":
            return
        if state.get(name) == "open":
            cycles.append(trail[trail.index(name):] + [name])
            return
        state[name] = "open"
        for dep in phases[name].upstream():
            if dep in phases:
                walk(dep, trail + [name])
        state[name] = "done"

    for n in phases:
        walk(n, [])
    return cycles


def execution_order(phases: dict[str, Phase], live_only: bool = True) -> list[str]:
    """Phases in an order that satisfies their dependencies.

    Ties are broken by name so the order is stable between runs — an unstable
    order would make the generated README churn."""
    todo = {n: p for n, p in phases.items() if p.is_live or not live_only}
    done, order = set(), []

    def rank(n: str) -> tuple:
        # Tooling last: export and build steps have no upstream phase, so a
        # plain alphabetical tie-break would otherwise put them first, which
        # reads as though the pipeline begins by building the manuscript.
        return (phases[n].status == "tooling", n)

    while todo:
        ready = sorted((n for n, p in todo.items()
                        if all(d in done or d not in todo for d in p.upstream())),
                       key=rank)
        if not ready:                     # a cycle; emit the rest deterministically
            order.extend(sorted(todo))
            break
        for n in ready:
            order.append(n); done.add(n); todo.pop(n)
    return order


def missing_inputs(phase: Phase, phases: dict[str, Phase],
                   drive: str | Path) -> list[str]:
    """Declared inputs that are not on disk.

    Checked before a phase runs, so a missing artefact costs a second rather
    than forty minutes of Colab."""
    drive = Path(drive)
    missing = []
    for item in phase.inputs:
        head = item.split(".")[0]
        artefact = item[len(head) + 1:] if head in phases else item
        if not artefact:
            continue
        if not (drive / artefact).exists():
            missing.append(item)
    return missing


GROUP_TITLES = {
    "data":        "Data preparation — acquisition, network, masks, products",
    "inversion":   "Inversion — every attempt to recover per-pixel displacement (H1)",
    "corrections": "Corrections — atmosphere and viewing geometry",
    "hypotheses":  "Hypothesis tests — zone contrast, aggregation, hydrology (H2–H4)",
    "robustness":  "Robustness — falsification, external controls, referee rounds",
    "manuscript":  "Manuscript — figures, tables and the assembled document",
}


def readme_table(phases: dict[str, Phase]) -> str:
    """The phase table, generated. Replaces the hand-written one that
    documented 14 of 38.

    One sub-table per group. A single 38-row table was sorted by dependency,
    which is the right order to *run* the pipeline and the wrong one to read it:
    it interleaved the H1 inversions with the data preparation and left no way
    to see that phaseD, phaseG and phaseI are one argument."""
    marks = {"current": "✅", "superseded": "⛔ superseded",
             "exploratory": "🔍 exploratory", "tooling": "🔧 tooling"}
    order = execution_order(phases, live_only=False)
    out: list[str] = []
    for group, heading in GROUP_TITLES.items():
        members = [n for n in order if phases[n].group == group]
        if not members:
            continue
        out += [f"### {GROUPS[group]} — {heading}", "",
                "| Phase | Question | Status | Paper |", "|---|---|---|---|"]
        for name in members:
            ph = phases[name]
            q = ph.question if ph.question != "-" else "—"
            out.append(f"| [`{name}`]({ph.notebook}) — {ph.title} | {q} | "
                       f"{marks.get(ph.status, ph.status)} | {ph.paper} |")
        out.append("")
    # A phase whose group is unknown would otherwise vanish from the table
    # silently; `validate` reports it, and it still gets listed here.
    rest = [n for n in order if phases[n].group not in GROUP_TITLES]
    if rest:
        out += ["### Ungrouped", "",
                "| Phase | Question | Status | Paper |", "|---|---|---|---|"]
        for name in rest:
            ph = phases[name]
            q = ph.question if ph.question != "-" else "—"
            out.append(f"| [`{name}`]({ph.notebook}) — {ph.title} | {q} | "
                       f"{marks.get(ph.status, ph.status)} | {ph.paper} |")
        out.append("")
    return "\n".join(out).rstrip()


def summary(phases: dict[str, Phase]) -> dict:
    counts: dict[str, int] = {}
    for p in phases.values():
        counts[p.status] = counts.get(p.status, 0) + 1
    return {"n_phases": len(phases), "by_status": counts,
            "order": execution_order(phases)}
