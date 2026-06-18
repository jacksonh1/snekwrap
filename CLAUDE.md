
## Working style

### Push back when something doesn't make sense

If a request is scientifically or technically wrong — the wrong GROMACS flag, a misunderstanding of T-REMD, a workflow that would silently produce bad results — say so clearly before implementing it. A wrong simulation or analysis run wastes real compute time and can produce results that look plausible but are wrong. See the **Known gotchas** section below for specific discovered pitfalls.

### Suggest better alternatives

If there's a cleaner, faster, or more correct approach than what was asked for, say so and explain the tradeoff. Don't just implement what was asked if something clearly better exists. Include enough context for an informed decision.

### Prefer simple and scalable solutions

Solutions should work for 8, 48, or 128 replicas without special-casing. Prefer shell/Python idioms that stay readable as the codebase grows. If a task has a five-line solution and a fifty-line solution, understand why the complexity is or isn't justified before recommending it. Don't add abstractions or generality beyond what the current task requires.

---

## Code style: fail loudly

This project follows a "fail loudly" philosophy. Bugs that crash immediately are strongly preferred over bugs that silently produce wrong results.

### Core philosophy

- **Crashes are cheap; silent bugs are expensive.** Prefer code that crashes obviously when assumptions are violated over code that "handles" the violation by producing degraded output.
- **Don't paper over uncertainty.** If you're unsure whether a value can be None, empty, or wrong-typed, either ask, add an assertion, or leave a clearly-marked comment — never add a default to make the question go away.
- **Make illegal states unrepresentable.** Prefer types and structures where the invalid case can't be expressed, over runtime checks for the invalid case.

### Error handling

- **No bare `except:` or `except Exception:` clauses** unless the exception is logged AND re-raised, or the recovery is documented and intentional.
- **Don't catch exceptions just to log and continue.** If the operation failed, the caller needs to know.
- **No `.get(key, default)` patterns** unless the default is semantically meaningful, not just a way to avoid a `KeyError`.
- **No `value or fallback` shortcuts** (`x or []`, `x or {}`, `x or 0`) unless `None`/empty/zero is genuinely interchangeable with the fallback. These hide bugs where `x` was unexpectedly empty.
- **Don't add defensive `if x is not None:` checks** unless `None` is a real expected case. If `None` would be a bug, let it crash.

### Indexing and iteration

- **Prefer iteration over indexing.** Use `for item in items`, not `for i in range(len(items))`. When you need the index too, use `enumerate`.
- **Use `zip(strict=True)`** (Python 3.10+) so mismatched-length iterables crash instead of silently truncating.
- **Assert invariants before code that relies on them.** E.g., `assert len(a) == len(b)` before zipping when the lengths must match.

### Types and structure

- **Use dataclasses or TypedDicts, not raw dicts**, when the shape matters and is fixed.
- **Parse, don't validate, at boundaries.** Convert untrusted input into typed structures at the edge; the rest of the code should be able to assume the data is valid.

### When in doubt

- **Ask before adding error handling.** If tempted to wrap something in try/except, ask what the intended behavior is when it fails.
- **Flag assumptions explicitly.** If making an assumption about input shape, range, or type that isn't enforced by the types, leave a comment like `# ASSUMES: items is non-empty`.

---

## Known gotchas

Specific pitfalls discovered in this environment. **When a new pitfall is discovered during work — a surprising behavior, a cluster quirk, a wrong assumption that caused a failure — add it here immediately without waiting to be asked.** The goal is that the same mistake is never made twice.

### SLURM: sbatch scripts are copied to a temp path at execution time

SLURM copies the `.sbatch` script to a temporary location before running it, so `BASH_SOURCE[0]` inside an sbatch script does not point to the repo. Code like `$(dirname "${BASH_SOURCE[0]}")/../site_config.sh` will silently resolve to the wrong place or fail.

**Fix (two options):**
- Use `$SLURM_SUBMIT_DIR` — a SLURM-provided environment variable that always holds the directory from which `sbatch` was called. Works as long as the job is submitted from the repo root (the normal case).
- Or pass the repo path explicitly via `--export` (e.g. `GROMACS_SCRIPTS_DIR`) and use that variable inside the sbatch.

Regular shell scripts called from *outside* SLURM (e.g. `analysis_scripts/`) can use `BASH_SOURCE[0]` reliably.

