# Development Workflow

> **When to read:** setting up a fresh checkout, running tests/lint, or releasing a
> version. **Expire:** update as tooling changes.

## Environment

Workflow targets [uv](https://docs.astral.sh/uv/). There is no `Makefile` and no CI yet —
these commands are the canonical equivalents:

```bash
uv sync                                    # create .venv, install deps + package + lqc CLI
uv run ruff check                            # lint
uv run pytest
uv run lqc -b tmp/data/ENCFF417VHJ.chr22.sorted.bam -o out   # smoke run
```

Dependencies are all runtime: `numpy`, `pandas`, `matplotlib`, `pysam` (from `pyproject.toml`).
`ruff` is intended as a dev-only tool. *why: numpy/pandas/matplotlib/pysam are imported by
the package · expire: if the dependency set changes.*

## Running & testing without a full BAM

- The unit tests need only `tests/data/cs_test.test_data` (committed) — no reference
  FASTA, no BAM, no network. They assert the cs path and the CIGAR+MD path produce
  identical results (see `tests/test_cs.py`).
- `tmp/data/` (gitignored) holds one real indexed BAM (≈200 MB) for manual smoke runs. It
  is **not** part of automated testing. *why: too large and slow for CI · expire: if CI
  re-uses it via an artifact cache.*
- Regenerate the CS fixture with `scripts/generate_cs_test_data.py <genome.fa>
  <cs.bam> <md.bam> <out_prefix>`; the committed `.test_data` is generated, not handwritten.

## Packaging & versioning

- **Version is single-sourced** in `src/lqc/__init__.py` (`__version__`). Bump only that file;
  hatchling reads it for the build (`[tool.hatch.version] path = "src/lqc/__init__.py"`) and
  `src/lqc/cli.py` reads it for `--version`.
  *why: removes the old `setup.py` + `bin/lqc` duplication · when: releasing · expire: if
  versioning moves to a VCS tag/SCM provider.*
- Packaging lives in `pyproject.toml` (hatchling). Non-Python assets under `src/lqc/`
  (`template/*.html`, `template/*.svg`) ship automatically; no `package_data` config is
  needed. *why: hatchling includes everything under the package by default · when: adding
  assets · expire: if a build backend with different data-file rules is adopted.*

## Style

- No linter is currently enforced in-repo; adopt `uv run ruff check` as the gate (see
  AGENTS.md hard constraint #11). Keep existing PEP-8-ish spacing, single quotes, and
  `name = value` argument style. *why: consistency with the existing ~3k-line codebase ·
  when: writing new code · expire: once `restructuredtext`/`ruff format` is agreed and
  applied tree-wide.*