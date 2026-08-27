# Architecture Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the structural debt documented in `docs/specs/2026-08-27-architecture-cleanup-design.md` — collapse duplicated boilerplate, unify the five-stat tuple and column-name contract, and drive the CLI output with declarative specs — without changing any output, column name, filename, or runtime behavior.

**Architecture:** Small, independently-green refactors. New shared code goes in one of two new modules (`src/lqc/constants.py`, `src/lqc/_base.py`). Each accumulator class keeps its own `__add__` but inherits the label guard and `__radd__` from `_LabelledStat`. `stat.py` gains a `ContigStats` NamedTuple that becomes the single representation of "the five stats". The composition root `cli.py` swaps its imperative output section for declarative spec lists.

**Tech Stack:** Python 3.11+, pytest, ruff (0.16.4), pandas, numpy, matplotlib (Agg), pysam.

---

## File Structure

- Create `src/lqc/constants.py` — shared constants (`TOTAL_LABEL`, `MISMATCH_TYPES`).
- Create `src/lqc/_base.py` — `_LabelledStat` base + `concat_stats()` helper.
- Modify `src/lqc/readstat.py`, `indel.py`, `mismatch.py`, `splice.py` — inherit the base.
- Modify `src/lqc/stat.py` — `ContigStats` NamedTuple, return it, use `open_alignment_file`.
- Modify `src/lqc/utils.py` — `open_alignment_file` context manager; legacy docstring.
- Modify `src/lqc/report_table.py` — column constants + `total_row()`; import `MISMATCH_TYPES`.
- Modify `src/lqc/report_html.py` — `_html_add_indel_table`, `_replace_tokens`; use constants.
- Modify `src/lqc/report_figure.py` — import `MISMATCH_TYPES`.
- Modify `src/lqc/cli.py` — declarative figure/table emission; use `ContigStats`/`concat_stats`/`total_row`.
- Modify `src/lqc/__init__.py` — lazy imports via `__getattr__`.
- Modify `docs/architecture.md` — sync pipeline to M0+ reality.
- Create `tests/test_constants.py`, `tests/test_base.py`; extend `tests/test_stat.py`, `tests/test_report_html.py`, `tests/test_report_table.py`, `tests/test_cli.py`, `tests/test_regressions.py`, `tests/test_utils.py`.

---

## Task 1: `src/lqc/constants.py`

**Files:**
- Create: `src/lqc/constants.py`
- Create: `tests/test_constants.py`
- Modify: `src/lqc/report_table.py` (import `MISMATCH_TYPES`)
- Modify: `src/lqc/report_figure.py` (import `MISMATCH_TYPES`)

- [ ] **Step 1: Write the failing test**

Create `tests/test_constants.py`:

```python
from lqc import constants
from lqc.report_table import MISMATCH_TYPES_IN_ORDER


def test_total_label_sentinel():
    assert constants.TOTAL_LABEL == 'Total'


def test_mismatch_types_order():
    assert len(constants.MISMATCH_TYPES) == 12
    assert constants.MISMATCH_TYPES[0] == 'ac'
    assert constants.MISMATCH_TYPES[-1] == 'tg'


def test_report_table_uses_shared_mismatch_types():
    assert tuple(MISMATCH_TYPES_IN_ORDER) == constants.MISMATCH_TYPES
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_constants.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'lqc.constants'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/lqc/constants.py`:

```python
"""Shared constants used across the parsing, stat, and reporting layers."""

TOTAL_LABEL = 'Total'

MISMATCH_TYPES = (
    'ac', 'ag', 'at', 'ca', 'cg', 'ct',
    'ga', 'gc', 'gt', 'ta', 'tc', 'tg',
)
```

In `src/lqc/report_table.py`, replace the module-level definition:

```python
MISMATCH_TYPES_IN_ORDER = ['ac', 'ag', 'at', 'ca', 'cg', 'ct',
                           'ga', 'gc', 'gt', 'ta', 'tc', 'tg']
```

with:

```python
from lqc.constants import MISMATCH_TYPES

MISMATCH_TYPES_IN_ORDER = list(MISMATCH_TYPES)
```

In `src/lqc/report_figure.py`, find the local list in `plot_mismatch_type_count` (around line 796):

```python
    mistypes = ["ac", "ag", "at",
                "ca", "cg", "ct",
                "ga", "gc", "gt",
                "ta", "tc", "tg"]
```

and replace with:

```python
    mistypes = [t for t in MISMATCH_TYPES]
```

Add the import at the top of `src/lqc/report_figure.py`:

```python
from lqc.constants import MISMATCH_TYPES
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_constants.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/constants.py tests/test_constants.py src/lqc/report_table.py src/lqc/report_figure.py
git commit -m "refactor: add shared constants module"
```

---

## Task 2: `_LabelledStat` base class

**Files:**
- Create: `src/lqc/_base.py`
- Create: `tests/test_base.py`
- Modify: `src/lqc/readstat.py`, `src/lqc/indel.py`, `src/lqc/mismatch.py`, `src/lqc/splice.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_base.py`:

```python
import pytest

from lqc._base import _LabelledStat
from lqc.mismatch import Mismatch
from lqc.readstat import ReadStat


def test_labelled_stat_rejects_nonstring_label():
    with pytest.raises(TypeError):
        _LabelledStat(123)


def test_labelled_stat_radd_zero_returns_self():
    stat = _LabelledStat('x')
    assert (0 + stat) is stat


def test_readstat_still_rejects_nonstring_label():
    with pytest.raises(TypeError):
        ReadStat(123)


def test_readstat_add_wrong_type_raises_typeerror():
    with pytest.raises(TypeError):
        ReadStat('a') + Mismatch('b')
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_base.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'lqc._base'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/lqc/_base.py`:

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

In each of the four classes, replace the `__init__` body:

```python
    def __init__(self, label = ''):
        super().__init__()
        if isinstance(label, str):
            pass
        else:
            raise TypeError(
                "label should be string."
            )
        self.label = label
```

with:

```python
    def __init__(self, label = ''):
        super().__init__(label)
```

Add the import to each class file:

```python
from lqc._base import _LabelledStat
```

and change the class declaration to inherit, e.g. `class ReadStat(_LabelledStat):` (likewise `Indel`, `Mismatch`, `Splice`).

In each class, delete the `__radd__` method entirely, and in `__add__` replace:

```python
        assert isinstance(other, type(self)),\
            'wrong object to add'
```

with:

```python
        other = self._require_same_type(other)
```

For `Mismatch.__add__` the inline form `assert isinstance(other, type(self)), 'wrong object to add'` becomes the same `other = self._require_same_type(other)` line.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_base.py tests/test_readstat.py tests/test_indel.py tests/test_mismatch.py tests/test_splice.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/_base.py tests/test_base.py src/lqc/readstat.py src/lqc/indel.py src/lqc/mismatch.py src/lqc/splice.py
git commit -m "refactor: shared labelled-stat base class"
```

---

## Task 3: `concat_stats` + drop the `relabel()` patch

**Files:**
- Modify: `src/lqc/_base.py` (add `concat_stats`)
- Modify: `tests/test_base.py` (add tests)
- Modify: `src/lqc/cli.py` (use `concat_stats` for the five Total objects)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_base.py`:

```python
from lqc._base import concat_stats


def test_concat_stats_sums_and_relabels():
    a = ReadStat('chr1')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    b = ReadStat('chr2')
    b.add_read(200, insertion=0, deletion=1, mismatch=0, intron=1)
    total = concat_stats([a, b], label='Total')
    assert total.get_read_count() == 2
    assert total.label == 'Total'
    assert a.label == 'chr1'
    assert b.label == 'chr2'


def test_concat_stats_single_element_does_not_alias():
    a = ReadStat('chr1')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    total = concat_stats([a], label='Total')
    assert total is not a
    assert total.label == 'Total'
    assert a.label == 'chr1'
    assert total.get_read_count() == 1


def test_concat_stats_requires_nonempty():
    with pytest.raises(ValueError):
        concat_stats([], label='Total')
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_base.py -k concat_stats -v`
Expected: FAIL with `ImportError: cannot import name 'concat_stats'`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/lqc/_base.py`:

```python
from copy import copy


def concat_stats(iterable, label):
    """Fold accumulator objects (via ``__add__``) into one relabelled object.

    Always returns a new object so relabelling never aliases a caller's
    per-contig instance (the exact bug the old ``sum([x]) is x`` path guarded
    against). Only the first element is shallow-copied, sharing the read-only
    numpy arrays rather than duplicating them.
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

Add to the top-level import block of `src/lqc/cli.py`:

```python
from lqc._base import concat_stats
from lqc.constants import TOTAL_LABEL
```

Remove the now-unused `from copy import copy` import (line 12).

In `src/lqc/cli.py`, replace the `sum(...)` + `relabel` block (the 12 lines beginning with the `relabel` definition comment and ending with `ssplice = relabel(sum(l_splice))`) with:

```python
    message = 'Sum of statistics from each contig.'
    logger.debug(message)
    sreadstat = concat_stats(l_readstat, TOTAL_LABEL)
    sinsertion = concat_stats(l_insertion, TOTAL_LABEL)
    sdeletion = concat_stats(l_deletion, TOTAL_LABEL)
    smismatch = concat_stats(l_mismatch, TOTAL_LABEL)
    ssplice = concat_stats(l_splice, TOTAL_LABEL)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_base.py -k concat_stats -v && uv run pytest tests/test_cli.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/_base.py tests/test_base.py src/lqc/cli.py
git commit -m "refactor: concat_stats helper replaces relabel patch"
```

---

## Task 4: `ContigStats` NamedTuple + delete `get_stat_list`

**Files:**
- Modify: `src/lqc/stat.py` (add `ContigStats`, return it)
- Modify: `src/lqc/cli.py` (attribute access, delete `get_stat_list`)
- Modify: `tests/test_stat.py` (add a named-field test)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_stat.py`:

```python
def test_stat_element_returns_contigstats_named_fields(cs_bam):
    cs = stat_element_from_bam_by_contig(cs_bam, None, 'chr1', 'cs')
    assert cs._fields == (
        'readstat', 'insertion', 'deletion', 'mismatch', 'splice'
    )
    assert cs.readstat.get_read_count() == 2
    assert cs.insertion.label == 'chr1'
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_stat.py::test_stat_element_returns_contigstats_named_fields -v`
Expected: FAIL with `AttributeError` (plain tuple has no `_fields`).

- [ ] **Step 3: Write minimal implementation**

In `src/lqc/stat.py`, after the `StatBlock` definition add:

```python
class ContigStats(NamedTuple):
    """The five per-contig accumulators, in canonical order."""
    readstat: ReadStat
    insertion: Indel
    deletion: Indel
    mismatch: Mismatch
    splice: Splice
```

Change `stat_element_from_bam_by_contig`'s `return readstat, insertion, deletion, mismatch, splice` to:

```python
    return ContigStats(readstat, insertion, deletion, mismatch, splice)
```

Change `reduce_blocks_to_contigs`'s `result.append((readstat, insertion, deletion, mismatch, splice))` to:

```python
        result.append(ContigStats(
            readstat, insertion, deletion, mismatch, splice
        ))
```

In `src/lqc/cli.py`, delete the entire `get_stat_list` function (lines 63–79) and replace the five calls:

```python
    l_readstat = get_stat_list(result, 'readstat')
    l_insertion = get_stat_list(result, 'insertion')
    l_deletion = get_stat_list(result, 'deletion')
    l_mismatch = get_stat_list(result, 'mismatch')
    l_splice = get_stat_list(result, 'splice')
```

with:

```python
    l_readstat = [cs.readstat for cs in result]
    l_insertion = [cs.insertion for cs in result]
    l_deletion = [cs.deletion for cs in result]
    l_mismatch = [cs.mismatch for cs in result]
    l_splice = [cs.splice for cs in result]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_stat.py tests/test_cli.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/stat.py src/lqc/cli.py tests/test_stat.py
git commit -m "refactor: ContigStats named tuple for per-contig stats"
```

---

## Task 5: `open_alignment_file` context manager

**Files:**
- Modify: `src/lqc/utils.py` (add context manager, refactor its 3 open sites)
- Modify: `src/lqc/stat.py` (refactor its 3 open sites)
- Modify: `tests/test_utils.py` (add tests)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_utils.py`:

```python
from lqc.utils import open_alignment_file


def test_open_alignment_file_yields_and_closes_bam(cs_bam):
    with open_alignment_file(cs_bam) as bam:
        assert bam is not None
        assert list(bam.references) == ['chr1']
    # Exiting the context manager must close the file handle.
    assert bam.closed
```

(The `cs_bam` fixture is injected by name by pytest; no import needed. `test_utils.py` already imports `pytest`.)

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_utils.py::test_open_alignment_file_yields_and_closes_bam -v`
Expected: FAIL with `ImportError: cannot import name 'open_alignment_file'`.

- [ ] **Step 3: Write minimal implementation**

In `src/lqc/utils.py`, add:

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

Refactor the three call sites in `utils.py`:

`list_bam_contigs` → replace:

```python
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    contigs = list(bam.references)
    bam.close()
    return contigs
```

with:

```python
    with open_alignment_file(bam_file) as bam:
        return list(bam.references)
```

`check_bam_with_cs_or_md` → replace its open (`file_type = ... file_read = ... bam = pysam.AlignmentFile(...)`) with `with open_alignment_file(bam_file) as bam:` and change its final `bam.close()` to remove it (the `with` handles it); keep the `return bam_type` statements and the enumeration body inside the `with`.

`write_readcs` → replace its open/close pair with `with open_alignment_file(bam_file) as bam:` (indent the loop body), removing `bam.close()`.

In `src/lqc/stat.py`, import `open_alignment_file` (extend the `from lqc.utils import ...` line) and refactor the three sites (`prefetch_records`, `plan_tasks`, `stat_region`) the same way — `plan_tasks` and `stat_region` currently wrap the loop in `try/finally: bam.close()`, which the `with` replaces.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_utils.py tests/test_stat.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/utils.py src/lqc/stat.py tests/test_utils.py
git commit -m "refactor: open_alignment_file context manager"
```

---

## Task 6: `report_html.py` de-duplication

**Files:**
- Modify: `src/lqc/report_html.py`
- Modify: `tests/test_report_html.py` (add tests)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_html.py`:

```python
from lqc.report_html import (
    _replace_tokens,
    html_add_deletion_table,
    html_add_insertion_table,
)


def test_replace_tokens_substitutes_all_in_order():
    html = _replace_tokens(
        '<p>{%a%}-{%b%}</p>', {'{%a%}': '1', '{%b%}': '2'}
    )
    assert html == '<p>1-2</p>'


def test_indel_tables_share_implementation():
    import pandas as pd
    table = pd.DataFrame([{
        'label': 'Total', 'total_count': 1, 'total_length': 2,
        'mean_length': 2.0, 'median_length': 2.0,
    }])
    ins = html_add_insertion_table('<h>{%insertion_table%}</h>', table, 0.0)
    dele = html_add_deletion_table('<h>{%deletion_table%}</h>', table, 0.0)
    assert '{%insertion_table%}' not in ins
    assert '{%deletion_table%}' not in dele
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_report_html.py::test_replace_tokens_substitutes_all_in_order -v`
Expected: FAIL with `ImportError: cannot import name '_replace_tokens'`.

- [ ] **Step 3: Write minimal implementation**

In `src/lqc/report_html.py`, add after the imports:

```python
def _replace_tokens(html_string, tokens):
    """Apply each ``{token: replacement}`` via ``re.sub`` in insertion order."""
    result = html_string
    for token, replacement in tokens.items():
        result = re.sub(token, replacement, result)
    return result
```

Add the shared indel implementation and thin wrappers. Replace the two existing
`html_add_insertion_table` / `html_add_deletion_table` definitions with:

```python
def _html_add_indel_table(html_string, indel_table, per_kb, kind):
    token = kind  # 'insertion' or 'deletion'
    rowstring_list = []
    total_count = total_length = 0
    mean_length = median_length = 0.0
    for _, row in indel_table.iterrows():
        tmprow_list = [
            '<th scope="row">{}</th>'.format(row['label']),
            '<td>{}</td>'.format(row['total_count']),
            '<td>{}</td>'.format(row['total_length']),
            '<td>{:.4g}</td>'.format(float(row['mean_length'])),
            '<td>{:.4g}</td>'.format(float(row['median_length']))
        ]
        if row['label'] == TOTAL_LABEL:
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + "</tr>"
            )
            total_count = row['total_count']
            total_length = row['total_length']
            mean_length = float(row['mean_length'])
            median_length = float(row['median_length'])
        else:
            rowstring_list.append(
                "<tr>" + '\n'.join(tmprow_list) + "</tr>"
            )

    return _replace_tokens(html_string, {
        r"\{%" + token + r"_total_" + token + r"_number%\}":
            f'{total_count}',
        r"\{%" + token + r"_total_" + token + r"_length%\}":
            f'{total_length}',
        r"\{%" + token + r"_mean_length%\}": f'{mean_length:.4g}',
        r"\{%" + token + r"_median_length%\}": f'{median_length:.4g}',
        r"\{%" + token + r"_mean_" + token + r"_per_read_per_kb%\}":
            f'{per_kb:.4g}',
        r"\{%" + token + r"_table%\}": '\n'.join(rowstring_list),
    })


def html_add_insertion_table(html_string, insertion_table,
                             mean_insertion_per_read_per_kb):
    return _html_add_indel_table(
        html_string, insertion_table,
        mean_insertion_per_read_per_kb, 'insertion'
    )


def html_add_deletion_table(html_string, deletion_table,
                            mean_deletion_per_read_per_kb):
    return _html_add_indel_table(
        html_string, deletion_table,
        mean_deletion_per_read_per_kb, 'deletion'
    )
```

Add `from lqc.constants import TOTAL_LABEL` to `report_html.py`. Use the same
`_replace_tokens` helper to collapse the scalar-token chains in
`html_add_readstat_table`, `html_add_mismatch_table`, `html_add_splice_table`,
and `html_add_mapping` (leave `html_add_bootstrap`, `html_add_data`, and
`inline_figures` using their callable `re.sub` as-is). Finally replace
`html_add_mismatch_table`'s local `mis_types = ['ac', ...]` with
`mis_types = list(MISMATCH_TYPES)` (import `MISMATCH_TYPES`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_html.py tests/test_report_html_assets.py tests/test_regressions.py -v`
Expected: PASS (the byte-level asset/parity tests protect the exact HTML output).

- [ ] **Step 5: Commit**

```bash
git add src/lqc/report_html.py tests/test_report_html.py
git commit -m "refactor: deduplicate HTML table renderers"
```

---

## Task 7: Column-name contract + `total_row`

**Files:**
- Modify: `src/lqc/report_table.py` (column constants + `total_row`)
- Modify: `src/lqc/cli.py` (use `total_row`)
- Modify: `tests/test_report_table.py` (add tests)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_table.py`:

```python
from lqc.report_table import total_row, COL_LABEL, COL_READ_COUNT


def test_total_row_pulls_total_cell():
    import pandas as pd
    df = pd.DataFrame([
        {'label': 'chr1', 'read_count': 5},
        {'label': 'Total', 'read_count': 7},
    ])
    assert total_row(df, 'read_count') == 7
    assert COL_LABEL == 'label'
    assert COL_READ_COUNT == 'read_count'
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_report_table.py::test_total_row_pulls_total_cell -v`
Expected: FAIL with `ImportError: cannot import name 'total_row'`.

- [ ] **Step 3: Write minimal implementation**

In `src/lqc/report_table.py`, add at the top (after imports):

```python
from lqc.constants import TOTAL_LABEL

COL_LABEL = 'label'
COL_READ_COUNT = 'read_count'


def total_row(table, column):
    """Return the ``column`` value from the ``Total`` row of a summary table."""
    return table.loc[table[COL_LABEL] == TOTAL_LABEL, column].iloc[0]
```

In `src/lqc/cli.py`, replace the four `.loc[...].values[0]` call sites with
`total_row(t_readstat, column)` using the column names currently in those
expressions (`'mean_mismatch_per_read_per_kb'`, `'mean_insertion_per_read_per_kb'`,
`'mean_deletion_per_read_per_kb'`, `'mean_intron_per_read'`), and add
`total_row` to the `lqc.report_table` import list. New code should read, e.g.:

```python
    new_html_string = html_add_mismatch_table(
        new_html_string, t_mismatch,
        smismatch.get_total_count(),
        total_row(t_readstat, 'mean_mismatch_per_read_per_kb'),
        mismatch_type_counter
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_table.py tests/test_cli.py tests/test_report_html.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/report_table.py src/lqc/cli.py tests/test_report_table.py
git commit -m "refactor: column-name contract and total_row helper"
```

---

## Task 8: Data-driven figure/table emission in `cli.py`

**Files:**
- Modify: `src/lqc/cli.py`
- Modify: `tests/test_cli.py` (add a fixture-level smoke assertion if useful)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_cli.py`:

```python
def test_main_outputs_expected_figure_files(cs_bam, tmp_path):
    out = tmp_path / 'out'
    rc = main(['-b', cs_bam, '-o', str(out), '-t', '1',
               '-c', 'chr1', '--log-level', 'DEBUG'])
    assert rc == 0
    produced = {p.name for p in (out / 'fig').glob('*')}
    for expected in ['readstat_bar_Read_count.png',
                     'readstat_bar_mean_element_per_read.Total.png',
                     'insertion_bar_count.png',
                     'splice_type.png']:
        assert expected in produced
```

(`tests/test_cli.py` already imports `main` from `lqc.cli`.)

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_cli.py::test_main_outputs_expected_figure_files -v`
Expected: PASS on current code (this is a *characterization* test pinning the existing filenames before the rewrite). Then proceed; after the rewrite it must still pass. If it does not pass on current code, stop and fix the expected filenames in the test first.

- [ ] **Step 3: Rewrite the figure/table emission**

Delete every `f_*` key from `o_files` (lines 268–360). `o_files` keeps only `pickle`, `cs`, the seven `t_*` keys, and `html`.

Replace the entire `# plot figures` section (from the `message = 'Output figures.'` line through the end of the `plot_mapping_*` block, ending just before `message = 'Output figures finished.'`) and the table-write block with declarative loops. Use module-level spec lists (place them near `savefig`, above `main`):

```python
READSTAT_BAR_FEATURES = [
    'Read count', 'Median read length', 'Mean read length',
    'Insertions per read', 'Insertions per read per kb',
    'Deletions per read', 'Deletions per read per kb',
    'Mismatches per read', 'Mismatches per read per kb',
    'Mean intron number', 'N50', 'L50',
]

MULTI_FIG_SPECS = [
    ('readstat_bar_mean_element_per_read', plot_readstat_bar_mean_element_per_read, 'readstat'),
    ('readstat_bar_mean_element_per_read_per_kb', plot_readstat_bar_mean_element_per_read_per_kb, 'readstat'),
    ('readstat_line_cumulative_length', plot_readstat_cumulative_length, 'readstat'),
    ('readstat_bar_ratio_with_element', plot_readstat_bar_ratio_with_element, 'readstat'),
    ('readstat_hist_length', plot_readstat_length_hist, 'readstat'),
    ('insertion_hist_length', plot_indel_hist_length, 'insertion'),
    ('deletion_hist_length', plot_indel_hist_length, 'deletion'),
    ('insertion_hist_location', plot_indel_hist_location, 'insertion'),
    ('deletion_hist_location', plot_indel_hist_location, 'deletion'),
    ('mismatch_type', plot_mismatch_type_count, 'mismatch'),
    ('mismatch_hist_location', plot_mismatch_hist_location, 'mismatch'),
    ('splice_type', plot_splice_type_count, 'splice'),
    ('mapping_hist_mapq', plot_mapping_mapq_hist, 'readstat'),
    ('mapping_hist_aligned_fraction', plot_mapping_aligned_fraction_hist, 'readstat'),
    ('mapping_scatter_aligned_vs_query', plot_mapping_aligned_vs_query, 'readstat'),
]

ELEMENT_BAR_SPECS = [
    ('insertion_bar_count', 'insertion', 'Insertion'),
    ('deletion_bar_count', 'deletion', 'Deletion'),
    ('mismatch_bar_count', 'mismatch', 'Mismatch'),
    ('intron_bar_count', 'splice', 'Intron'),
]
```

Inside `main`, after building the `l_*` and `s*` objects, add:

```python
    data_sets = {
        'readstat': (l_readstat, sreadstat),
        'insertion': (l_insertion, sinsertion),
        'deletion': (l_deletion, sdeletion),
        'mismatch': (l_mismatch, smismatch),
        'splice': (l_splice, ssplice),
    }
```

and replace the figure/table blocks with:

```python
    ####################
    # write output tables
    for key, table in [
            ('t_readstat', t_readstat),
            ('t_insertion', t_insertion),
            ('t_deletion', t_deletion),
            ('t_mismatch', t_mismatch),
            ('t_splice', t_splice),
            ('t_mapping', t_mapping),
            ('t_splice_all', t_splice_all),
    ]:
        table.to_csv(
            o_files[key], sep = '\t', index = False,
            float_format = FLOAT_FORMAT
        )

    ####################
    # plot figures
    figdir = o_dirs['fig']
    for feature in READSTAT_BAR_FEATURES:
        fig = plot_readstat_bar(l_readstat, feature)
        savefig(fig, os.path.join(
            figdir, 'readstat_bar_' + feature.replace(' ', '_')
        ))
        plt.close('all')

    for stem, plot_func, data_key in MULTI_FIG_SPECS:
        data_list, data_sum = data_sets[data_key]
        generate_multiple_figs(
            plot_func, data_list, data_sum,
            os.path.join(figdir, stem), width = 5, height = 4
        )

    for stem, data_key, kind in ELEMENT_BAR_SPECS:
        fig = plot_element_total_count(data_sets[data_key][0], kind)
        savefig(fig, os.path.join(figdir, stem))
        plt.close('all')
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_cli.py -v && uv run ruff check src/lqc/cli.py`
Expected: PASS / clean.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/cli.py tests/test_cli.py
git commit -m "refactor: declarative figure and table emission in CLI"
```

---

## Task 9: Lazy imports in `__init__.py`

**Files:**
- Modify: `src/lqc/__init__.py`
- Create: `tests/test_lazy_imports.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_lazy_imports.py`:

```python
import lqc


def test_all_names_resolve_lazily():
    for name in lqc.__all__:
        assert getattr(lqc, name) is not None


def test_unknown_attribute_raises():
    import pytest
    with pytest.raises(AttributeError):
        lqc.__getattr__('definitely_not_real')
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_lazy_imports.py -v`
Expected: the `test_unknown_attribute_raises` test FAILs (currently `__getattr__` does not exist, but Python falls back to module default → `lqc.__getattr__` raises `AttributeError` already; the missing part is `__getattr__`). If it passes on current code, it is a characterization/regression test; proceed.

- [ ] **Step 3: Write minimal implementation**

Replace the eager `from lqc.* import ...` bodies in `src/lqc/__init__.py` with a lazy loader, keeping `__version__` and `__all__`. Use this skeleton:

```python
import importlib

__version__ = "0.0.8"

__all__ = [
    'CS', 'Indel', 'Mismatch', 'Splice',
    'bam_or_sam', 'check_bam_with_cs_or_md', 'copy_logo',
    'create_indel_summary_table', 'create_mapping_table',
    'create_mismatch_normalized_read_location_table',
    'create_readstat_table', 'create_splice_all_table',
    'create_splice_table',
    'get_html_template', 'html_add_bootstrap', 'html_add_data',
    'html_add_deletion_table', 'html_add_insertion_table',
    'html_add_mapping', 'html_add_mismatch_table',
    'html_add_readstat_table', 'html_add_splice_table',
    'inline_figures', 'list_bam_contigs',
    'plot_element_total_count', 'plot_indel_hist_length',
    'plot_indel_hist_location', 'plot_mapping_aligned_fraction_hist',
    'plot_mapping_aligned_vs_query', 'plot_mapping_mapq_hist',
    'plot_mismatch_hist_location', 'plot_mismatch_type_count',
    'plot_readstat_bar', 'plot_readstat_bar_mean_element_per_read',
    'plot_readstat_bar_mean_element_per_read_per_kb',
    'plot_readstat_bar_ratio_with_element',
    'plot_readstat_cumulative_length', 'plot_readstat_length_hist',
    'plot_splice_type_count', 'stat_element_from_bam_by_contig',
    'write_readcs',
]

_LAZY = {
    'CS': 'lqc.cs',
    'Indel': 'lqc.indel',
    'Mismatch': 'lqc.mismatch',
    'Splice': 'lqc.splice',
    'bam_or_sam': 'lqc.utils',
    'check_bam_with_cs_or_md': 'lqc.utils',
    'list_bam_contigs': 'lqc.utils',
    'write_readcs': 'lqc.utils',
    'stat_element_from_bam_by_contig': 'lqc.stat',
    'copy_logo': 'lqc.template',
    'get_html_template': 'lqc.template',
    'create_indel_summary_table': 'lqc.report_table',
    'create_mapping_table': 'lqc.report_table',
    'create_mismatch_normalized_read_location_table': 'lqc.report_table',
    'create_readstat_table': 'lqc.report_table',
    'create_splice_all_table': 'lqc.report_table',
    'create_splice_table': 'lqc.report_table',
    'html_add_bootstrap': 'lqc.report_html',
    'html_add_data': 'lqc.report_html',
    'html_add_deletion_table': 'lqc.report_html',
    'html_add_insertion_table': 'lqc.report_html',
    'html_add_mapping': 'lqc.report_html',
    'html_add_mismatch_table': 'lqc.report_html',
    'html_add_readstat_table': 'lqc.report_html',
    'html_add_splice_table': 'lqc.report_html',
    'inline_figures': 'lqc.report_html',
    'plot_element_total_count': 'lqc.report_figure',
    'plot_indel_hist_length': 'lqc.report_figure',
    'plot_indel_hist_location': 'lqc.report_figure',
    'plot_mapping_aligned_fraction_hist': 'lqc.report_figure',
    'plot_mapping_aligned_vs_query': 'lqc.report_figure',
    'plot_mapping_mapq_hist': 'lqc.report_figure',
    'plot_mismatch_hist_location': 'lqc.report_figure',
    'plot_mismatch_type_count': 'lqc.report_figure',
    'plot_readstat_bar': 'lqc.report_figure',
    'plot_readstat_bar_mean_element_per_read': 'lqc.report_figure',
    'plot_readstat_bar_mean_element_per_read_per_kb': 'lqc.report_figure',
    'plot_readstat_bar_ratio_with_element': 'lqc.report_figure',
    'plot_readstat_cumulative_length': 'lqc.report_figure',
    'plot_readstat_length_hist': 'lqc.report_figure',
    'plot_splice_type_count': 'lqc.report_figure',
}


def __getattr__(name):
    module_name = _LAZY.get(name)
    if module_name is not None:
        value = getattr(importlib.import_module(module_name), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lqc' has no attribute {name!r}")
```

Every `__all__` name maps to the module whose attribute has the identical name;
`__all__` remains the source of truth and the test iterates it.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_lazy_imports.py tests/test_cli.py tests/test_report_html_assets.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/__init__.py tests/test_lazy_imports.py
git commit -m "refactor: lazy imports for the lqc package surface"
```

---

## Task 10: Legacy code docstrings + `docs/architecture.md` sync

**Files:**
- Modify: `src/lqc/stat.py` (`prefetch_records` docstring)
- Modify: `src/lqc/utils.py` (`write_readcs` docstring)
- Modify: `docs/architecture.md`

- [ ] **Step 1: Fix `prefetch_records` docstring**

Change its first paragraph:

```python
    """Return ``[(contig, [ReadRecord, ...]), ...]`` in BAM traversal order.

    Materializes every read's ``ReadRecord`` in memory before pooling. For the
    cs path this is small (a cs string plus a few ints); for the MD path each
    record additionally carries the full ``query_sequence``.
    """
```

to:

```python
    """Return ``[(contig, [ReadRecord, ...]), ...]`` in BAM traversal order.

    Serial/test helper: the M0+ CLI path uses ``plan_tasks`` + ``stat_region``
    (worker-side fetch) instead of materializing every read up front. Keep for
    ``tests/`` and ``tmp/profile_lqc.py``.
    """
```

- [ ] **Step 2: Fix `write_readcs` docstring**

Replace its current body lead sentence (which says "The CLI no longer calls this") with a note that it remains the documented library API for a full-BAM dump, mirroring the worker-side `read.cs` emission in `stat_records`.

- [ ] **Step 3: Sync `docs/architecture.md`**

In the "Pipeline" section, replace steps 3–4 (the `prefetch_records()` + `_chunk` + `mp.Pool` description) with the M0+ flow: `plan_tasks()` splits each contig into coordinate windows, `mp.Pool.map(stat_region)` fetches/accumulates each window, and `reduce_blocks_to_contigs()` folds `StatBlock`s back into per-contig `ContigStats`. Also update the module map's `stat.py` row (key symbols: `ReadRecord`, `StatBlock`, `ContigStats`, `plan_tasks`, `stat_region`, `reduce_blocks_to_contigs`, `stat_element_from_bam_by_contig`) and the `src/lqc/cli.py` note if it references `_chunk`.

- [ ] **Step 4: Run checks**

Run: `uv run ruff check && uv run pytest -q`
Expected: clean / all pass.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/stat.py src/lqc/utils.py docs/architecture.md
git commit -m "docs: sync architecture and legacy docstrings"
```

---

## Final Verification

- [ ] Run `uv run ruff check` — clean.
- [ ] Run `uv run pytest -q` — all green.
- [ ] If a real BAM exists under `tmp/`, run the CLI and diff tables/figures/HTML
      against the pre-refactor baseline under `tmp/verify/` to confirm byte parity.