# Environment

Two files, deliberately.

| File | Use |
|---|---|
| `requirements.txt` | Loose. Use when you intend to move the stack forward. |
| `requirements-lock.txt` | Exact versions that produced the published results. **Use this to reproduce.** |

## Why the lock file matters

`colab_setup.sh` previously ran `pip install --upgrade` with no constraints, so every
Colab session installed whatever was current that day. A rerun months later ran on a
different stack than the one behind the published numbers, and nothing recorded which
stack that had been.

For a paper whose contribution is partly methodological, that is the reproducibility
gap a referee is most entitled to ask about.

## Generating the lock file

Once, from a working Colab session that has produced good results:

```bash
!bash environment/colab_setup.sh          # install the loose set
!pip freeze > environment/requirements-lock.txt
```

Then commit it, and note the date in `DECISIONS.md`. Regenerate deliberately — never
as a side effect.

## Per-run capture

Every archived run also records the versions actually in use at run time
(`run_archive.environment_snapshot`). The lock file says what *should* be installed;
the run manifest says what *was*. When a number moves, `compare_runs` tells you whether
the environment changed.
