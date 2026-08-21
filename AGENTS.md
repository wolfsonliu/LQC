# AGENTS.md

## Project Overview

LQC is a Python 3.11+ command-line tool that generates quality-control reports — summary
tables, figures, and an HTML page — from SAM/BAM alignments of long-read RNA-seq data
(PacBio, Oxford Nanopore). It parses the `cs` or `MD` tag of each read to report read
statistics, insertions, deletions, mismatches, and splice sites.

## Quick Start

Development uses [uv](https://docs.astral.sh/uv/). There is no `Makefile` or CI yet.

- Environment + deps + CLI entry point: `uv sync`
- Lint: `uvx ruff check`  (append `--fix` to auto-fix what it can)
- Tests: `uv run pytest`
- Smoke-run the CLI: `uv run lqc -b <cs-tagged.sorted.indexed.bam> -o out`

## Hard Constraints

1. **Python 3.11+** — do not use newer syntax or stdlib features.
   *(why: `pyproject.toml` sets `requires-python = ">=3.11"` · when: writing any code · expire: when that floor is raised)*

2. **Input alignments must carry `cs` tags or `MD` tags**; MD-only input also requires `--genome-fasta`.
   *(why: `check_bam_with_cs_or_md` rejects otherwise, `stat.py` needs a tag path · when: accepting/running any input · expire: when a new tag path is supported)*

3. **`cs` wins over `MD` when both are present** (the `'both'` case is downgraded to `'cs'`).
   *(why: `cs` already encodes splice sites, so the genome FASTA can be skipped · when: touching tag detection in `lqc/cli.py` or `lqc/stat.py` · expire: if precedence is deliberately reversed)*

4. **Input BAM must be sorted and indexed** (`.bai`/`.csi` present).
   *(why: `pysam.AlignmentFile.fetch(contig)` requires an index · when: running against real BAMs, writing fixtures · expire: if `fetch` is replaced by an index-free scan)*

5. **matplotlib backend must remain `Agg`** (headless).
   *(why: no display on servers; set in `lqc/cli.py` and `lqc/report_figure.py` · when: adding/editing any figure code · expire: if a display backend becomes required and safe)*

6. **`stat_element_from_bam_by_contig(...)` accepts `method` only in `['cs', 'MD', 'both']`.**
   *(why: `assert` guard in `lqc/stat.py` · when: calling stat functions directly (library use) · expire: when the API changes)*

7. **`label` on `ReadStat`/`Indel`/`Mismatch`/`Splice` must be a `str`; `'Total'` is reserved for the aggregate row.**
   *(why: constructors raise `TypeError` otherwise; `'Total'` is the sentinel across tables, figures, and HTML · when: creating/merging stats · expire: if label semantics change)*

8. **Summary-table column names are a contract** — downstream code reads specific columns by name (e.g. `mean_mismatch_per_read_per_kb`, `mean_intron_per_read`).
   *(why: `report_html.py` does `.loc[...]` lookups on those names · when: changing `lqc/report_table.py` columns or wiring · expire: when the name-based coupling is removed)*

9. **Parsed coordinates are 0-based half-open (`[low, high)`) on the reference**; any 1-based display conversion happens only at the boundary.
   *(why: defined in `lqc/cs.py`/`lqc/stat.py` · when: touching positions in `cs.py`, `stat.py`, `utils.py` · expire: only if the whole pipeline switches conventions, in one coordinated change)*

10. **License is GPLv3+** — any vendored or derived code must remain GPL-compatible.
    *(why: `LICENSE` and `pyproject.toml` classifier · when: adding third-party or copied code · expire: if the project is relicensed)*

11. **Code must pass `uvx ruff check`.**
    *(why: lint is the agreed quality gate for the upcoming improvement work · when: before finishing any change · expire: replaced by another agreed lint command)*

## Topic Docs

- **Architecture & data flow** — `docs/architecture.md`. Read when you need the big picture, are adding a feature, or tracing a cross-module bug.
- **Input & tag conventions** — `docs/input-and-tags.md`. Read before touching `lqc/cs.py`, `lqc/stat.py`, or `lqc/utils.py`, or when diagnosing parsing.
- **Reporting & output format** — `docs/reporting-and-output.md`. Read before touching `lqc/report_table.py`, `lqc/report_figure.py`, `lqc/report_html.py`, or `lqc/template/`.
- **Development workflow** — `docs/development.md`. Read when setting up a fresh checkout, running tests/lint, or releasing a version.