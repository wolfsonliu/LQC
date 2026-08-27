# Architecture Cleanup (Maintainability) Design

## Problem Statement

LQC went through a performance optimization pass (M0–M3) that moved the pipeline
from a materializing prefetch to a worker-side fetch (`plan_tasks` + `stat_region`),
and from boxed Python lists to numpy columnar storage. That work left a residue of
structural debt around it:

- `cli.py` (852 lines) is a "composition root" whose output-artifact wiring is
  largely declarative but written as imperative boilerplate, with a 40+ entry
  `o_files` dict whose string keys must match f-strings elsewhere by accident
  (dict keys even contain spaces, e.g. `'f_readstat_bar_Read count'`).
- The four accumulator classes (`ReadStat`, `Indel`, `Mismatch`, `Splice`) each
  duplicate the same `label` guard and `__radd__`, with no shared base.
- The `'Total'` aggregate is built in the composition root via a `copy.copy` +
  relabel hack that exists only because `sum([x]) is x`.
- The 12-element mismatch-type list is hardcoded in three modules, and the summary
  table column names form a cross-module contract (AGENTS rule 8) with no single
  source of truth.
- `report_html.py` has near-verbatim copy-paste (`html_add_insertion_table` vs
  `html_add_deletion_table`) plus a repetitive `re.sub` placeholder chain.
- BAM open/close boilerplate is copied six times; `prefetch_records` and
  `write_readcs` are half-dead relative to the M0 worker path; and
  `docs/architecture.md` still describes the pre-M0 pipeline.

This work is **maintainability only**: no column names, output filenames,
numbers, merge formulas, or runtime hot paths change. Every refactor must be
behavior-preserving, verified by the existing byte-parity / regression suite.

## Proposed Solution

A sequence of small, independently-green refactors that (1) collapse duplicated
boilerplate into named single sources of truth, (2) rewrite the composition
root's figure/table emission as declarative specs, and (3) reconcile the "five
stat objects" tuple ordering and the column-name contract. Development follows
TDD (test-first for every new abstraction) and is executed by subagents against
a detailed implementation plan.

## Detailed Design

### ① Single source of truth — new `src/lqc/constants.py`

New module holding the cross-module contract values:

```python
"""Shared constants used across the parsing, stat, and reporting layers."""
TOTAL_LABEL = 'Total'
MISMATCH_TYPES = (
    'ac', 'ag', 'at', 'ca', 'cg', 'ct',
    'ga', 'gc', 'gt', 'ta', 'tc', 'tg',
)
```

- `report_table.MISMATCH_TYPES_IN_ORDER` → imports `MISMATCH_TYPES` as its value
  (kept as a module attribute so existing readers are unaffected).
- `report_html.mis_types` and the `report_figure` local list are replaced by the
  shared `MISMATCH_TYPES`.
- The `'Total'` sentinel strings in `cli.py` / `report_html.py` (and the plan's
  tests) use `TOTAL_LABEL` (preserves AGENTS rule 7's reserved word).

### ② Shared accumulator base — new `src/lqc/_base.py`

```python
"""Internal base for the four labelled accumulator classes."""
class _LabelledStat:
    def __init__(self, label=''):
        if not isinstance(label, str):
            raise TypeError('label should be string.')
        self.label = label

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)

    def _require_same_type(self, other):
        if not isinstance(other, type(self)):
            raise TypeError('wrong object to add')
        return other
```

`ReadStat`/`Indel`/`Mismatch`/`Splice` call `super().__init__(label)`, delete
their own guard and `__radd__`. `__add__` uses `other = self._require_same_type(other)`
instead of the inline `assert` (same failure, cleaner message).

### ③ Unify the "five stats" tuple — `src/lqc/stat.py` + `cli.py`

Three places currently encode the same 5-element ordering implicitly: the
`StatBlock` NamedTuple fields, the bare tuple returned by
`reduce_blocks_to_contigs` / `stat_element_from_bam_by_contig`, and
`cli.get_stat_list`'s `0..4` int map. Add:

```python
class ContigStats(NamedTuple):
    readstat: ReadStat
    insertion: Indel
    deletion: Indel
    mismatch: Mismatch
    splice: Splice
```

- Both `reduce_blocks_to_contigs` and `stat_element_from_bam_by_contig` return
  `ContigStats(...)`.
- `cli.py` uses attribute access (`[cs.readstat for cs in result]`) and deletes
  `get_stat_list` entirely.
- `ContigStats` is a `NamedTuple` (still a tuple), so `test_stat.py`'s
  `reduced[0][0]`-style indexing keeps working unchanged.

### ④ Remove the `relabel()` patch — `concat_stats` in `src/lqc/_base.py`

```python
from copy import copy

def concat_stats(iterable, label):
    """Fold accumulator objects (via ``__add__``) into one relabelled object.

    Always returns a new object so relabelling never aliases a caller's
    per-contig instance (the exact bug the old ``sum([x]) is x`` path guarded
    against). Elements are shallow-copied only for the first element to avoid
    duplicating read-only numpy arrays.
    """
    it = iter(iterable)
    try:
        acc = copy(next(it))
    except StopIteration:
        raise ValueError('concat_stats requires at least one object') from None
    for obj in it:
        acc = acc + obj
    acc.label = label
    return acc
```

`cli.py` builds the five `'Total'` objects with `concat_stats(l_readstat, TOTAL_LABEL)`
etc., removing `relabel` and `sum` from the composition root. `reduce_blocks_to_contigs`
keeps its `sum(...)` fold (chunk objects are throwaway); `test_regressions.py`
still verifies `sum([x]) is x` for library semantics.

### ⑤ BAM open context manager — `src/lqc/utils.py`

```python
from contextlib import contextmanager

@contextmanager
def open_alignment_file(path):
    """Yield an opened pysam.AlignmentFile in the mode matching the file type."""
    fh = pysam.AlignmentFile(
        path, 'rb' if bam_or_sam(path) == 'BAM' else 'r'
    )
    try:
        yield fh
    finally:
        fh.close()
```

Replaces the six copied open/close blocks in `utils.py` (3×) and `stat.py` (3×).

### ⑥ Data-driven figure/table emission — `src/lqc/cli.py`

- Readstat 12 bar features: a module-level spec list of `(feature, stem)` pairs
  (the stems reproduce every existing key verbatim; note six "Insertions /
  Deletions / Mismatches" stems are lowercase while the feature strings are
  capitalized, so the stem is spelled out explicitly rather than derived),
  then `plot_readstat_bar(l_readstat, feature)` per pair.
- `generate_multiple_figs` plots: a module-level spec list
  `(stem, plot_func, data_key)` plus a per-run `data_sets` dict (key →
  `(list, sum)`) built in `main`, and one loop.
- Element bar counts: a spec list `(stem, data_key, kind)` and one loop.
- Tables: a `(o_files_key, dataframe)` list driven by one `to_csv` loop.

`o_files` shrinks to the non-bulk artifacts (pickle, read.cs, html). No filename,
plot, or number changes; the ~200 lines removed are pure boilerplate.

### ⑦ `report_html.py` de-duplication

- `_html_add_indel_table(html_string, table, kind)` shared by the two near-identical
  `html_add_insertion_table` / `html_add_deletion_table` functions (kept as thin
  public wrappers for the unchanged `__all__` API).
- `_replace_tokens(html_string, mapping)` wraps the scalar `{%token%}` `re.sub`
  chains; callable-based replacements (`html_add_bootstrap`, `html_add_data`,
  `inline_figures`) stay as-is.

### ⑧ Column-name contract + lazy imports

- `src/lqc/report_table.py` gains module-level column-name constants for every
  column `report_html.py` or `cli.py` reads across modules, plus
  `def total_row(table, column): return table.loc[table[TOTAL_LABEL_COL] == TOTAL_LABEL, column].iloc[0]`.
  `cli.py`'s four `.loc[label=='Total', col].values[0]` call sites use it.
- `src/lqc/__init__.py` switches to PEP 562 lazy imports (`__getattr__`) so
  `import lqc` no longer eagerly imports pandas/matplotlib/pysam; `__all__` and
  `from lqc import X` semantics are preserved.

### ⑨ Legacy code + doc drift

- `prefetch_records` (`stat.py`) — keep (used by `tests/` and `tmp/profile_lqc.py`)
  but fix its docstring to say "serial/test helper, not the M0+ CLI path".
- `write_readcs` (`utils.py`) — fix docstring; optionally extract the shared
  per-read cs-line writer used by `stat_records`.
- `docs/architecture.md` — update the pipeline section (lines 13–18) to
  `plan_tasks` + `stat_region`, and fix the module map's stale `stat.py` symbols
  and `tmp/data/` entry.

## Commit Sequence

Each step lands green and independently:

1. ① + ② — constants + base class.
2. ③ + ④ — `ContigStats` + `concat_stats`.
3. ⑤ + ⑦ + ⑧ — BAM context manager, HTML de-dup, column contract.
4. ⑥ — composition-root data-driven emission (largest; byte-parity gate).
5. ⑧-b lazy import + ⑨ doc/legacy.

## Success Criteria

- `uv run ruff check` clean (AGENTS rule 11).
- `uv run pytest` green, in particular `test_real_bam_smoke.py`,
  `test_regressions.py`, `test_stat.py`, `test_report_html.py`,
  `test_report_table.py`.
- With a real BAM under `tmp/`, a CLI run is byte-identical (tables / figures /
  HTML) to the pre-refactor baseline, per the M1–M3 `tmp/verify` convention.
- `docs/architecture.md` matches the implemented module boundaries.

## Out of Scope

- Renaming any summary-table column or output filename (rule 8 contract).
- Touching M0–M3 numpy layout / merge formulas / any perf hot path.
- Splitting `report_figure.py` / `cs.py`; type-hint overhaul; CI/Makefile/packaging.

## Open Questions

- Keep `formatting.py` as-is (with `FLOAT_FORMAT`) rather than folding it into
  `constants.py` — **recommended: keep separate**.
- Keep `write_readcs` exported vs. fold its cs-line writer with `stat_records` —
  recommended: keep exported, share the writer only if trivial.