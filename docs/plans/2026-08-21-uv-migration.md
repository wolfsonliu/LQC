# uv Package-Management Migration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `setup.py` packaging with a hatchling-built `pyproject.toml`, expose the CLI via a `lqc.cli:main` entry point, single-source the version, and adopt `uv sync`/`uv build` with a committed `uv.lock`.

**Architecture:** Packaging metadata moves to PEP 621 `pyproject.toml` (hatchling backend). The CLI body moves verbatim from `bin/lqc` into `lqc/cli.py:main()` (the existing 4-space-indented body becomes the function body with no re-indent). The version lives once in `lqc/__init__.py` and is read by hatchling (dynamic version) and the CLI (via `from lqc import __version__ as VERSION`).

**Tech Stack:** Python 3.11, uv 0.12, hatchling, ruff (transient via `uvx`), stdlib `unittest`.

**Spec:** `docs/specs/2026-08-21-uv-migration-design.md`

---

## File Structure

| Action | Path | Responsibility |
| --- | --- | --- |
| Create | `pyproject.toml` | PEP 621 metadata, build system, entry point, ruff target |
| Create | `.python-version` | pin `3.11` for `uv sync` |
| Create | `lqc/cli.py` | `main(argv=None) -> int`; moved from `bin/lqc` |
| Create | `lqc/test/test_cli.py` | regression test for the version flag |
| Modify | `lqc/__init__.py` | add `__version__ = "0.0.5"` |
| Modify | `.gitignore` | add `.venv/` |
| Modify | `README.md` | uv dev install + Python 3.11 badge/text |
| Modify | `AGENTS.md`, `docs/{architecture,reporting-and-output,development}.md` | replace stale `bin/lqc` / Python-version references |
| Delete | `setup.py` | fully replaced by `pyproject.toml` |
| Delete | `bin/lqc` | moved to `lqc/cli.py` |
| Generate | `uv.lock` | created by `uv sync`, committed |

---

## Task 1: Add `__version__` to `lqc/__init__.py`

**Files:**
- Modify: `lqc/__init__.py:1`

- [ ] **Step 1: Add the version constant**

Edit `lqc/__init__.py`, inserting a first line before the existing `from lqc.cs import CS`:

```python
__version__ = "0.0.5"
```

So the top of the file becomes:

```python
__version__ = "0.0.5"

from lqc.cs import CS
```

- [ ] **Step 2: Verify**

Run: `grep -n "__version__" lqc/__init__.py`
Expected: line 1 shows `__version__ = "0.0.5"`.

- [ ] **Step 3: Commit**

```bash
git add lqc/__init__.py
git commit -m "chore: add __version__ to lqc package"
```

---

## Task 2: Add `pyproject.toml`, `.python-version`, update `.gitignore`, delete `setup.py`

**Files:**
- Create: `pyproject.toml`
- Create: `.python-version`
- Modify: `.gitignore`
- Delete: `setup.py`

- [ ] **Step 1: Create `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "lqc"
description = "The Long-read RNA-seq quality control software."
readme = "README.md"
license = "GPL-3.0-or-later"
requires-python = ">=3.11"
authors = [
    { name = "Zhiheng Liu", email = "wolfsonliu@live.com" },
]
dependencies = [
    "numpy",
    "pandas",
    "matplotlib",
    "pysam",
]
classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dynamic = ["version"]

[project.urls]
Homepage = "https://github.com/gxiaolab/LQC"

[project.scripts]
lqc = "lqc.cli:main"

[tool.hatch.version]
path = "lqc/__init__.py"

[tool.hatch.build.targets.wheel]
packages = ["lqc"]

[tool.ruff]
target-version = "py311"
```

- [ ] **Step 2: Create `.python-version`**

File content (single line, trailing newline):

```text
3.11
```

- [ ] **Step 3: Update `.gitignore`**

Replace the line `__pycache__/` with two lines:

```
__pycache__/
.venv/
```

- [ ] **Step 4: Delete `setup.py`**

```bash
git rm setup.py
```

- [ ] **Step 5: Verify files exist**

Run: `ls -la pyproject.toml .python-version && test ! -e setup.py && echo "setup.py removed"`
Expected: both files listed, `setup.py removed` printed.

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml .python-version .gitignore
git commit -m "build: switch packaging to pyproject.toml (hatchling, uv)"
```

---

## Task 3: Write the failing test for the CLI entry point

**Files:**
- Create: `lqc/test/test_cli.py`

> This test cannot run until dependencies are installed (Task 5); its "failing" state now is
> that `lqc.cli` does not exist yet.

- [ ] **Step 1: Create `lqc/test/test_cli.py`**

```python
import io
import unittest
from contextlib import redirect_stdout

from lqc.cli import main


class TestCli(unittest.TestCase):

    def test_version_flag(self):
        buf = io.StringIO()
        with self.assertRaises(SystemExit) as cm, redirect_stdout(buf):
            main(['--version'])
        self.assertEqual(cm.exception.code, 0)
        self.assertIn('lqc 0.0.5', buf.getvalue())


if __name__ == '__main__':
    unittest.main()
```

- [ ] **Step 2: Confirm it fails for the right reason**

Run: `python3 -m py_compile lqc/test/test_cli.py` (syntax only — no deps needed)
Expected: no output (syntax OK). Importing would fail with `ModuleNotFoundError: lqc.cli`
until Task 4, because `lqc/cli.py` does not exist yet.

- [ ] **Step 3: Commit**

```bash
git add lqc/test/test_cli.py
git commit -m "test: add CLI --version regression test (failing)"
```

---

## Task 4: Move `bin/lqc` into `lqc/cli.py` and delete `bin/lqc`

**Files:**
- Create: `lqc/cli.py` (via `git mv`)
- Delete: `bin/lqc`

- [ ] **Step 1: Move the file and stage the rename**

```bash
git mv bin/lqc lqc/cli.py
```

- [ ] **Step 2: Edit the header (docstring + `import sys`)**

In `lqc/cli.py`, replace:

```
#! /usr/bin/env python3

import os
import argparse
```

with:

```
"""Command-line interface for LQC.

All orchestration lives in main(). See README.md for usage.
"""

import os
import sys
import argparse
```

- [ ] **Step 3: Import the version**

Replace:

```
from lqc import write_readcs
import matplotlib
```

with:

```
from lqc import write_readcs
from lqc import __version__ as VERSION
import matplotlib
```

- [ ] **Step 4: Turn the `__main__` block into `main()`**

Replace:

```
if __name__ == '__main__':
    VMAJOR, VMINOR, VMICRO = 0, 0, 5
    VERSION = '{}.{}.{}'.format(VMAJOR, VMINOR, VMICRO)
```

with:

```
def main(argv = None):
```

- [ ] **Step 5: Make the program name deterministic**

Replace:

```
    parser = argparse.ArgumentParser(
        description='The Long-read RNA-seq quality control software.'
    )
```

with:

```
    parser = argparse.ArgumentParser(
        prog = 'lqc',
        description='The Long-read RNA-seq quality control software.'
    )
```

- [ ] **Step 6: Pass `argv` through to argparse**

Replace:

```
    args = vars(parser.parse_args())
```

with:

```
    args = vars(parser.parse_args(argv))
```

- [ ] **Step 7: Return 0 and add a runnable guard at the end**

Replace:

```
    message = 'All done!'
    logging.info(message)

########################################
```

with:

```
    message = 'All done!'
    logging.info(message)
    return 0


if __name__ == '__main__':
    sys.exit(main())


########################################
```

- [ ] **Step 8: Sanity-check the result**

Run: `python3 -m py_compile lqc/cli.py lqc/test/test_cli.py`
Expected: no output (both compile).

Run: `grep -n "def main\|prog = 'lqc'\|parse_args(argv)\|__version__ as VERSION\|return 0\|sys.exit(main())" lqc/cli.py`
Expected: six matches, one per edit above.

- [ ] **Step 9: Commit**

```bash
git add lqc/cli.py lqc/test/test_cli.py
git commit -m "refactor: move CLI into lqc.cli:main and add console entry point"
```

---

## Task 5: `uv sync`, run tests, verify the CLI, commit the lockfile

**Files:**
- Generate: `uv.lock`

> Needs network to fetch the Python 3.11 interpreter and dependency wheels. All of
> numpy/pandas/matplotlib/pysam ship `cp311` linux x86_64 wheels.

- [ ] **Step 1: Create the environment and lockfile**

Run: `uv sync`
Expected: exit 0; `.venv/` is created; `uv.lock` is written; output ends with the resolved
package list including `lqc`.

- [ ] **Step 2: Run the full test suite**

Run: `uv run python -m unittest lqc.test.test_cs lqc.test.test_cli -v`
Expected: all tests `ok` (exit 0), including `TestCli.test_version_flag`.

- [ ] **Step 3: Verify the installed entry point**

Run: `uv run lqc --version`
Expected: `lqc 0.0.5`.

- [ ] **Step 4: Verify the package imports and version is single-sourced**

Run: `uv run python -c "import lqc; print(lqc.__version__)"`
Expected: `0.0.5`.

- [ ] **Step 5: Commit the lockfile**

```bash
git add uv.lock
git commit -m "build: add uv.lock"
```

---

## Task 6: Update `README.md`

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Fix the Python version badge**

Replace:

```
[![](https://img.shields.io/badge/python-v3.6%2B-brightgreen)](https://www.python.org/)
```

with:

```
[![](https://img.shields.io/badge/python-v3.11%2B-brightgreen)](https://www.python.org/)
```

- [ ] **Step 2: Fix the dependency list Python line**

Replace:

```
* [python3.6+](https://www.python.org/): with os, sys, argparse, re,
```

with:

```
* [python3.11+](https://www.python.org/): with os, sys, argparse, re,
```

- [ ] **Step 3: Replace the "From github" install instructions with uv**

Replace the block from `### From github` through the `python setup.py install` code fence
(i.e. everything between the `### From github` heading and the `### From pip` heading):

```
### From github

Download from github:

```{bash}
git clone https://github.com/gxiaolab/LQC
cd LQC
```

Install the package:

```{bash}
python setup.py install
```
```

with:

```
### From github (development)

Clone and install with [uv](https://docs.astral.sh/uv/):

```{bash}
git clone https://github.com/gxiaolab/LQC
cd LQC
uv sync
```

This creates a `.venv` and installs the `lqc` command. Run the checks and the CLI:

```{bash}
uvx ruff check
uv run python -m unittest lqc.test.test_cs -v
uv run lqc -b <cs-tagged.sorted.indexed.bam> -o out
```
```

- [ ] **Step 4: Verify no stale install command remains**

Run: `grep -n "setup.py install\|v3.6\|python3.6" README.md`
Expected: no output.

- [ ] **Step 5: Commit**

```bash
git add README.md
git commit -m "docs: document uv install flow and Python 3.11 in README"
```

---

## Task 7: Update stale doc references (`bin/lqc` → `lqc/cli.py`, Python version)

**Files:**
- Modify: `AGENTS.md`
- Modify: `docs/architecture.md`
- Modify: `docs/reporting-and-output.md`
- Modify: `docs/development.md`

> Do NOT touch `docs/specs/` — the spec documents the "before" state on purpose.

- [ ] **Step 1: Replace `bin/lqc` with `lqc/cli.py` in the three docs + AGENTS.md**

For each of `AGENTS.md`, `docs/architecture.md`, `docs/reporting-and-output.md`,
replace every occurrence of `bin/lqc` with `lqc/cli.py` (there are none of the static
doc phrase in `docs/development.md` to change yet — see Step 3). In `AGENTS.md` this
changes these two lines:

```
   *(why: `cs` already encodes splice sites, so the genome FASTA can be skipped · when: touching tag detection in `bin/lqc` or `lqc/stat.py` · expire: if precedence is deliberately reversed)*
```
→
```
   *(why: `cs` already encodes splice sites, so the genome FASTA can be skipped · when: touching tag detection in `lqc/cli.py` or `lqc/stat.py` · expire: if precedence is deliberately reversed)*
```

```
   *(why: no display on servers; set in `bin/lqc` and `lqc/report_figure.py` · when: adding/editing any figure code · expire: if a display backend becomes required and safe)*
```
→
```
   *(why: no display on servers; set in `lqc/cli.py` and `lqc/report_figure.py` · when: adding/editing any figure code · expire: if a display backend becomes required and safe)*
```

- [ ] **Step 2: Update the Python version floor in `AGENTS.md`**

Replace the project-overview line and hard constraint #1:

```
LQC is a Python 3.7+ command-line tool that generates quality-control reports — summary
```
→
```
LQC is a Python 3.11+ command-line tool that generates quality-control reports — summary
```

```
1. **Python 3.7+** — do not use newer syntax or stdlib features.
   *(why: `setup.py` sets `python_requires='>=3.7'` · when: writing any code · expire: when that floor is raised)*
```
→
```
1. **Python 3.11+** — do not use newer syntax or stdlib features.
   *(why: `pyproject.toml` sets `requires-python = ">=3.11"` · when: writing any code · expire: when that floor is raised)*
```

- [ ] **Step 3: Fix the version bullet in `docs/development.md`**

Replace the two bullets that still reference `setup.py`/`bin/lqc`:

```
- **Version is duplicated in two files** and must be bumped together:
  `setup.py` (`VMAJOR,VMINOR,VMICRO`) and `bin/lqc` (hardcoded `VERSION`).
  *why: no shared constant exists · when: releasing · expire: when a single version source
  (e.g. `importlib.metadata` or `__version__`) is introduced.*
- Package data is declared in `setup.py` (`lqc.test`, `lqc.template`); adding a new
  `.html`/`.svg`/test-data file means adding it to `package_data`. *why: `zip_safe=False`
  and explicit package_data · expire: if packaging moves fully to `pyproject.toml`.*
- README claims "Python 3.6+" in a badge but `setup.py` sets `>=3.7`; trust `setup.py`.
```

with:

```
- **Version is single-sourced** in `lqc/__init__.py` (`__version__`). Bump only that file;
  hatchling reads it for the build and `lqc/cli.py` reads it for `--version`.
  *why: removes the old `setup.py` + `bin/lqc` duplication · when: releasing · expire: if
  versioning moves to a VCS tag/SCM provider.*
- Packaging lives in `pyproject.toml` (hatchling). Non-Python files under `lqc/` (template
  `.html`/`.svg`, test `.test_data`) are included automatically; no `package_data` config.
  *why: hatchling ships everything under the package by default · when: adding assets ·
  expire: if a build backend with different data-file rules is adopted.*
```

Also remove the now-obsolete "Transition note" paragraph (the blockquote that starts
`> **Transition note (current state):**`) and the trailing standalone sentence
`- README claims "Python 3.6+" ...` if any text remains after the replacements above.

- [ ] **Step 4: Verify no stale references remain**

Run: `grep -rn "bin/lqc\|3\.7+\|3\.6+\|python 3\.7\|setup.py install" AGENTS.md docs/architecture.md docs/reporting-and-output.md docs/development.md`
Expected: no output.

- [ ] **Step 5: Commit**

```bash
git add AGENTS.md docs/architecture.md docs/reporting-and-output.md docs/development.md
git commit -m "docs: update references for lqc.cli entry point and Python 3.11"
```

---

## Task 8: Final verification (`uv build`, `uvx ruff check`)

- [ ] **Step 1: Build distributions**

Run: `uv build`
Expected: exit 0; `dist/` contains `lqc-0.0.5.tar.gz` and `lqc-0.0.5-py3-none-any.whl`.

- [ ] **Step 2: Confirm the wheel carries the entry point**

Run: `python3 -c "import zipfile,glob; w=glob.glob('dist/*.whl')[0]; print(sorted(n for n in zipfile.ZipFile(w).namelist() if n.endswith('entry_points.txt') or 'template/' in n))"`
Expected: lists `lqc-0.0.5.dist-info/entry_points.txt` and `lqc/template/...` files.

- [ ] **Step 3: Lint the changed file**

Run: `uvx ruff check lqc/cli.py`
Expected: exit 0 (or only pre-existing-style findings, which are out of scope — do NOT
run `--fix` on the whole tree).

- [ ] **Step 4: Full test suite once more**

Run: `uv run python -m unittest lqc.test.test_cs lqc.test.test_cli -v`
Expected: all tests pass.

---

## Self-Review Notes

- **Spec coverage:** every design item maps to a task — pyproject (T2), entry point
  (T1+T4), version single-source (T1+T4), `.python-version`/`uv.lock` (T2+T5), README
  (T6), doc references (T7), out-of-scope items are explicitly not touched.
- **Consistency:** `__version__` is defined in T1 and consumed by T4 and T5; the entry
  point string `lqc = "lqc.cli:main"` matches the `def main(argv = None)` created in T4;
  the test in T3 targets the exact `lqc 0.0.5` output produced by `prog = 'lqc'` (T4).
- **Out of scope:** no pytest, no ruff rules/config expansion, no version bump, no
  `tests/` restructuring, no commit of unrelated uncommitted files (`AGENTS.md`/`docs`
  are committed only via the Task 7 doc refactoring and are otherwise left as-is;
  the pre-existing `.gitignore` `tmp/*` edit is preserved, not committed).