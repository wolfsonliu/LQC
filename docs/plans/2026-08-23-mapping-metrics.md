# Mapping Metrics & Output Polish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-read mapping-quality and aligned-fraction statistics (new `table/mapping.txt`, three figures, an HTML "Mapping" section), and in the same pass fix the output-quality issues from the review: round TSV floats, share a categorical palette, collapse `splice.txt` to four categories (full matrix → `splice_all.txt`), fix `mismatch.txt` type ordering + add `bin_total`, and resolve the per-kb normalization ambiguity (`X_per_kb` → `X_per_query_kb`; add `X_per_aligned_kb`).

**Architecture:** Extend the existing `ReadStat` to carry `mapping_quality` and `aligned_length` per read (Approach A). `stat.py` captures the two fields in `ReadRecord`. New `create_mapping_table` and three `plot_mapping_*` functions follow the existing table/figure patterns. A `FLOAT_FORMAT` constant is applied to all `to_csv` calls. All changes are additive except the three `*_per_kb` renames and the reshaped `splice.txt` columns, which are documented as external-contract breaks.

**Tech Stack:** Python 3.11+, pysam, pandas, matplotlib (Agg only), numpy, uv + pytest + ruff.

**Design spec:** `docs/specs/2026-08-23-mapping-metrics-design.md`

---

## File Structure

- **Create `src/lqc/formatting.py`** — `FLOAT_FORMAT` TSV float format constant.
- **Modify `src/lqc/readstat.py`** — store `mapping_quality`/`aligned_length` per read; new getters; sum aligned base in `__add__`; extend `__repr__`.
- **Modify `src/lqc/stat.py`** — add two `ReadRecord` fields; capture them in `record_from_read`; pass them in `process_record`.
- **Modify `src/lqc/report_table.py`** — new `create_mapping_table`, `create_splice_all_table`; rework `create_readstat_table`, `create_splice_table`, `create_mismatch_normalized_read_location_table`.
- **Modify `src/lqc/report_figure.py`** — shared palette in the two bar plotters; readable `plot_indel_hist_length`; three new `plot_mapping_*`.
- **Modify `src/lqc/report_html.py`** — new `html_add_mapping`; simplify `html_add_splice_table`.
- **Modify `src/lqc/template/template.html`** — "Mapping" nav + accordion section.
- **Modify `src/lqc/cli.py`** — `o_files`, build/write `mapping.txt` + `splice_all.txt`, `FLOAT_FORMAT` on TSV, generate mapping figures, `html_add_mapping`, `lqc-data` mapping entry.
- **Modify `src/lqc/__init__.py`** — export new symbols; keep `__all__` in sync.
- **Modify `docs/reporting-and-output.md`** — document new files/columns/tokens.
- **Tests** — `tests/test_readstat.py`, `tests/test_stat.py`, `tests/test_report_table.py`, `tests/test_report_figure.py`, `tests/test_report_html.py`, `tests/test_regressions.py`.

---

## Task 1: Extend `ReadStat` with mapping fields and getters

**Files:**
- Modify: `src/lqc/readstat.py`
- Test: `tests/test_readstat.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_readstat.py`:

```python
def test_mapping_fields_default_to_zero(rs):
    assert rs.get_mapping_qualities() == [0, 0]
    assert rs.get_aligned_lengths() == [0, 0]
    assert rs.get_aligned_fractions() == [0.0, 0.0]
    assert rs.get_total_aligned_base() == 0


def test_mapping_getters():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
               mapping_quality=60, aligned_length=180)
    r.add_read(100, insertion=0, deletion=0, mismatch=0, intron=0,
               mapping_quality=30, aligned_length=100)
    assert r.get_mapping_qualities() == [60, 30]
    assert r.get_aligned_lengths() == [180, 100]
    assert r.get_aligned_fractions() == [0.9, 1.0]
    assert r.get_total_aligned_base() == 280
    assert r.get_mean_mapping_quality() == 45
    assert r.get_median_mapping_quality() == 45
    assert r.get_mean_aligned_fraction() == pytest.approx(0.95)
    assert r.get_median_aligned_fraction() == pytest.approx(0.95)
    assert r.get_read_count_with_aligned_fraction_below(0.9) == 0
    assert r.get_read_count_with_aligned_fraction_below(0.95) == 1
    assert r.get_read_count_fully_aligned() == 1


def test_aligned_base_rates():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=0,
               mapping_quality=60, aligned_length=100)
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=0,
               mapping_quality=60, aligned_length=100)
    assert r.get_total_aligned_base() == 200
    assert r.insertions_per_aligned_base() == pytest.approx(4 / 200)
    assert r.deletions_per_aligned_base() == pytest.approx(8 / 200)
    assert r.mismatches_per_aligned_base() == pytest.approx(12 / 200)


def test_aligned_base_rate_zero_when_no_aligned():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=0)
    assert r.insertions_per_aligned_base() == 0
    assert r.deletions_per_aligned_base() == 0
    assert r.mismatches_per_aligned_base() == 0


def test_add_two_readstats_sums_aligned_base():
    a = ReadStat('a')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0,
               mapping_quality=60, aligned_length=90)
    b = ReadStat('b')
    b.add_read(200, insertion=0, deletion=1, mismatch=0, intron=1,
               mapping_quality=30, aligned_length=200)
    c = a + b
    assert c.get_total_aligned_base() == 290
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_readstat.py -v`
Expected: FAIL with `AttributeError: 'ReadStat' object has no attribute 'get_mapping_qualities'`

- [ ] **Step 3: Add `_total_aligned_base` and extend `add_read`**

In `src/lqc/readstat.py`, replace the `__init__` body and `add_read`:

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
        self._read_count = 0
        self._total_base = 0
        self._total_aligned_base = 0
        self._reads = []

    def add_read(self,
                 length,
                 insertion,
                 deletion,
                 mismatch,
                 intron,
                 mapping_quality = 0,
                 aligned_length = 0):
        self._read_count += 1
        self._reads.append([
            length,
            insertion, deletion,
            mismatch, intron,
            mapping_quality, aligned_length
        ])
        self._total_base += length
        self._total_aligned_base += aligned_length
```

- [ ] **Step 4: Add the new getters**

In `src/lqc/readstat.py`, insert the following block immediately before `_get_median` (the method that starts `def _get_median(self, item_list):`):

```python
    def get_mapping_qualities(self):
        return [a[5] for a in self._reads]

    def get_aligned_lengths(self):
        return [a[6] for a in self._reads]

    def get_aligned_fractions(self):
        return [
            a[6] / a[0] if a[0] > 0 else 0.0
            for a in self._reads
        ]

    def get_total_aligned_base(self):
        return self._total_aligned_base

    def get_mean_mapping_quality(self):
        qualities = self.get_mapping_qualities()
        if not qualities:
            return 0.0
        return sum(qualities) / len(qualities)

    def get_median_mapping_quality(self):
        return self._get_median(self.get_mapping_qualities())

    def get_mean_aligned_fraction(self):
        fractions = self.get_aligned_fractions()
        if not fractions:
            return 0.0
        return sum(fractions) / len(fractions)

    def get_median_aligned_fraction(self):
        return self._get_median(self.get_aligned_fractions())

    def get_read_count_with_aligned_fraction_below(self, threshold = 0.9):
        return sum(
            1 for a in self._reads
            if (a[6] / a[0] if a[0] > 0 else 0.0) < threshold
        )

    def get_read_count_fully_aligned(self):
        return sum(
            1 for a in self._reads
            if a[6] == a[0]
        )

    def insertions_per_aligned_base(self):
        if self._total_aligned_base == 0:
            return 0.0
        return sum(self.get_insertions()) / self._total_aligned_base

    def deletions_per_aligned_base(self):
        if self._total_aligned_base == 0:
            return 0.0
        return sum(self.get_deletions()) / self._total_aligned_base

    def mismatches_per_aligned_base(self):
        if self._total_aligned_base == 0:
            return 0.0
        return sum(self.get_mismatches()) / self._total_aligned_base
```

- [ ] **Step 5: Sum aligned base in `__add__`**

In `src/lqc/readstat.py`, replace the `__add__` method (lines 247-262) with:

```python
    def __add__(self, other):
        assert isinstance(other, type(self)),\
            'wrong object to add'
        sumReadStat = type(self)(
            f'{self.label} {other.label}'
        )
        sumReadStat._reads = deepcopy(
            self._reads + other._reads
        )
        sumReadStat._read_count = deepcopy(
            self._read_count + other._read_count
        )
        sumReadStat._total_base = deepcopy(
            self._total_base + other._total_base
        )
        sumReadStat._total_aligned_base = deepcopy(
            self._total_aligned_base + other._total_aligned_base
        )
        return sumReadStat
```

- [ ] **Step 6: Extend `__repr__`**

In `src/lqc/readstat.py`, in the `__repr__` list (lines 227-240) add two lines after the `mean of introns per read` entry:

```python
            f"  mean of aligned fraction: {float(self.get_mean_aligned_fraction()):.4}",
            f"  median of mapping quality: {float(self.get_median_mapping_quality()):.4}"
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `uv run pytest tests/test_readstat.py -v`
Expected: PASS (new tests and all pre-existing readstat tests)

- [ ] **Step 8: Commit**

```bash
git add src/lqc/readstat.py tests/test_readstat.py
git commit -m "feat: track mapping quality and aligned length in ReadStat"
```

---

## Task 2: Capture mapping fields in `ReadRecord`

**Files:**
- Modify: `src/lqc/stat.py`
- Test: `tests/test_stat.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_stat.py`:

```python
def test_stat_records_captures_mapping_fields(cs_bam):
    contig, records = prefetch_records(cs_bam, ['chr1'], 'cs')[0]
    read0 = records[0]
    assert read0.mapping_quality == 60
    assert read0.aligned_length == 10
    block = stat_records((0, contig, records), None, 'cs')
    assert block.readstat.get_mapping_qualities() == [60, 60]
    assert block.readstat.get_aligned_lengths() == [10, 10]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_stat.py::test_stat_records_captures_mapping_fields -v`
Expected: FAIL with `AttributeError: 'ReadRecord' object has no attribute 'mapping_quality'`

- [ ] **Step 3: Add the `ReadRecord` fields**

In `src/lqc/stat.py`, replace the `ReadRecord` field list (lines 16-32) so the two new required fields sit between `query_length` and `query_name`, before the defaults:

```python
class ReadRecord(NamedTuple):
    """Lightweight per-read fields needed to build a CS and accumulate stats.

    ``method`` is ``'cs'`` or ``'md'``. ``cs_string`` is set for the cs path;
    the remaining optional fields are set for the MD+CIGAR path.
    """
    method: str
    contig: str
    start_pos: int
    strand: str
    query_length: int
    mapping_quality: int
    aligned_length: int
    query_name: str
    cs_string: Optional[str] = None
    cigar: Optional[str] = None
    md_string: Optional[str] = None
    query_sequence: Optional[str] = None
    reference_end: Optional[int] = None
```

- [ ] **Step 4: Populate the fields in `record_from_read`**

In `src/lqc/stat.py`, replace the `record_from_read` function (lines 46-70) with:

```python
def record_from_read(read, contig, method):
    """Extract the minimal fields from a pysam read for stat accumulation."""
    strand = '-' if read.is_reverse else '+'
    if method == 'cs':
        return ReadRecord(
            method='cs',
            contig=contig,
            start_pos=read.reference_start,
            strand=strand,
            query_length=len(read.query_sequence),
            mapping_quality=read.mapping_quality,
            aligned_length=read.query_alignment_length,
            query_name=read.query_name,
            cs_string=read.get_tag('cs'),
        )
    return ReadRecord(
        method='md',
        contig=contig,
        start_pos=read.reference_start,
        strand=strand,
        query_length=len(read.query_sequence),
        mapping_quality=read.mapping_quality,
        aligned_length=read.query_alignment_length,
        query_name=read.query_name,
        cigar=read.cigarstring,
        md_string=read.get_tag('MD'),
        query_sequence=read.query_sequence,
        reference_end=read.reference_end,
    )
```

- [ ] **Step 5: Pass the fields in `process_record`**

In `src/lqc/stat.py`, replace the `readstat.add_read(...)` call inside `process_record` (lines 100-106) with:

```python
    readstat.add_read(
        length=record.query_length,
        insertion=cs.get_insertion_count(),
        deletion=cs.get_deletion_count(),
        mismatch=cs.get_mismatch_count(),
        intron=cs.get_intron_count(),
        mapping_quality=record.mapping_quality,
        aligned_length=record.aligned_length,
    )
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `uv run pytest tests/test_stat.py -v`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/lqc/stat.py tests/test_stat.py
git commit -m "feat: capture mapping quality and aligned length per read"
```

---

## Task 3: Add the TSV float format and apply it

**Files:**
- Create: `src/lqc/formatting.py`
- Modify: `src/lqc/cli.py`

- [ ] **Step 1: Create the formatting constant**

Create `src/lqc/formatting.py`:

```python
"""Shared output-formatting constants."""


# Significant digits for float columns in the TSV summary tables. Matches the
# ``{:.4}`` significant-digit formatting already used throughout report_html.py.
FLOAT_FORMAT = '%.4g'
```

- [ ] **Step 2: Apply `FLOAT_FORMAT` to every TSV write**

In `src/lqc/cli.py`, add the import near the other `from lqc import (...)` block:

```python
from lqc.formatting import FLOAT_FORMAT
```

Then update the five existing `to_csv` calls (lines 505-524) to pass `float_format=FLOAT_FORMAT`:

```python
    t_readstat.to_csv(
        o_files['t_readstat'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
    t_insertion.to_csv(
        o_files['t_insertion'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
    t_deletion.to_csv(
        o_files['t_deletion'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
    t_mismatch.to_csv(
        o_files['t_mismatch'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
    t_splice.to_csv(
        o_files['t_splice'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
```

- [ ] **Step 3: Run the fast suite to verify nothing broke**

Run: `uv run pytest tests/test_cli.py tests/test_report_table.py -v`
Expected: PASS (note: `to_csv(float_format=...)` does not affect int or string columns)

- [ ] **Step 4: Commit**

```bash
git add src/lqc/formatting.py src/lqc/cli.py
git commit -m "feat: round TSV float columns to 4 significant digits"
```

---

## Task 4: Add `create_mapping_table`

**Files:**
- Modify: `src/lqc/report_table.py`
- Modify: `src/lqc/__init__.py`
- Test: `tests/test_report_table.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_table.py`:

```python
def test_mapping_table_columns_and_values():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
                mapping_quality=60, aligned_length=180)
    r2 = ReadStat('chr2')
    r2.add_read(100, insertion=0, deletion=0, mismatch=0, intron=0,
                mapping_quality=30, aligned_length=100)
    total = sum([r1, r2])
    total.label = 'Total'
    df = create_mapping_table([r1, r2], total)

    assert list(df.columns) == [
        'label', 'read_count', 'query_base', 'aligned_base',
        'aligned_fraction_mean', 'aligned_fraction_median',
        'mapq_mean', 'mapq_median',
        'reads_aligned_fraction_lt_0.9', 'reads_fully_aligned'
    ]
    assert list(df['label']) == ['chr1', 'chr2', 'Total']
    t = df[df['label'] == 'Total'].iloc[0]
    assert t['read_count'] == 2
    assert t['query_base'] == 300
    assert t['aligned_base'] == 280
    assert t['aligned_fraction_mean'] == pytest.approx(0.95)
    assert t['aligned_fraction_median'] == pytest.approx(0.95)
    assert t['mapq_mean'] == pytest.approx(45)
    assert t['reads_aligned_fraction_lt_0.9'] == 0
    assert t['reads_fully_aligned'] == 1
```

Fix the import block at the top of `tests/test_report_table.py` (it currently imports `create_indel_summary_table, create_readstat_table, create_splice_table`) to also include `create_mapping_table`:

```python
from lqc.report_table import (
    create_indel_summary_table,
    create_mapping_table,
    create_readstat_table,
    create_splice_table,
)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_report_table.py::test_mapping_table_columns_and_values -v`
Expected: FAIL with `ImportError: cannot import name 'create_mapping_table'`

- [ ] **Step 3: Implement `create_mapping_table`**

Append to `src/lqc/report_table.py` (after `create_readstat_table`, before `create_mismatch_normalized_read_location_table`):

```python
def create_mapping_table(readstat_list, readstat_sum):
    colnames = [
        'label', 'read_count', 'query_base', 'aligned_base',
        'aligned_fraction_mean', 'aligned_fraction_median',
        'mapq_mean', 'mapq_median',
        'reads_aligned_fraction_lt_0.9', 'reads_fully_aligned'
    ]

    def _row(a):
        return [
            a.label,
            a.get_read_count(),
            a.get_total_base(),
            a.get_total_aligned_base(),
            float(a.get_mean_aligned_fraction()),
            float(a.get_median_aligned_fraction()),
            float(a.get_mean_mapping_quality()),
            float(a.get_median_mapping_quality()),
            a.get_read_count_with_aligned_fraction_below(0.9),
            a.get_read_count_fully_aligned()
        ]

    rows = [_row(a) for a in readstat_list] + [_row(readstat_sum)]
    return pd.DataFrame(rows, columns = colnames)
```

- [ ] **Step 4: Export `create_mapping_table`**

In `src/lqc/__init__.py`, add `create_mapping_table` to the `from lqc.report_table import (...)` block and to `__all__` (insert alphabetically next to the other `create_*` names).

- [ ] **Step 5: Run test to verify it passes**

Run: `uv run pytest tests/test_report_table.py::test_mapping_table_columns_and_values -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/lqc/report_table.py src/lqc/__init__.py tests/test_report_table.py
git commit -m "feat: add mapping summary table"
```

---

## Task 5: Rework `create_readstat_table` (normalization columns)

**Files:**
- Modify: `src/lqc/report_table.py`
- Test: `tests/test_report_table.py`

- [ ] **Step 1: Update the column-contract test**

In `tests/test_report_table.py`, replace the `READSTAT_COLUMNS` list (lines 11-22) with:

```python
READSTAT_COLUMNS = [
    'label', 'read_count', 'total_base', 'aligned_base',
    'read_length_mean', 'read_length_median',
    'read_length_N50', 'read_length_L50',
    'mean_insertion_per_read', 'mean_insertion_per_read_per_kb',
    'insertion_per_query_kb', 'insertion_per_aligned_kb',
    'mean_deletion_per_read', 'mean_deletion_per_read_per_kb',
    'deletion_per_query_kb', 'deletion_per_aligned_kb',
    'mean_mismatch_per_read', 'mean_mismatch_per_read_per_kb',
    'mismatch_per_query_kb', 'mismatch_per_aligned_kb',
    'mean_intron_per_read', 'mean_intron_per_read_per_kb',
]
```

Also update the `_two_readstat` helper (lines 25-32) to pass aligned fields, and add a value assertion for the new aligned columns:

```python
def _two_readstat():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
                mapping_quality=60, aligned_length=190)
    r2 = ReadStat('chr2')
    r2.add_read(200, insertion=4, deletion=8, mismatch=12, intron=16,
                mapping_quality=60, aligned_length=190)
    total = sum([r1, r2])
    total.label = 'Total'
    return r1, r2, total


def test_readstat_table_aligned_columns():
    r1, r2, total = _two_readstat()
    df = create_readstat_table([r1, r2], total)
    total_row = df[df['label'] == 'Total'].iloc[0]
    # 2 reads x aligned 190 = 380 aligned bases
    assert total_row['aligned_base'] == 380
    assert total_row['insertion_per_query_kb'] == pytest.approx(6 / 400 * 1000)
    assert total_row['insertion_per_aligned_kb'] == pytest.approx(6 / 380 * 1000)
    assert total_row['mismatch_per_aligned_kb'] == pytest.approx(18 / 380 * 1000)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_report_table.py::test_readstat_table_column_contract tests/test_report_table.py::test_readstat_table_aligned_columns -v`
Expected: FAIL (`create_readstat_table` still returns the old columns)

- [ ] **Step 3: Implement the reworked `create_readstat_table`**

In `src/lqc/report_table.py`, replace `create_readstat_table` (lines 4-71) with:

```python
def create_readstat_table(readstat_list, readstat_sum):
    colnames = [
        'label', 'read_count', 'total_base', 'aligned_base',
        'read_length_mean', 'read_length_median',
        'read_length_N50', 'read_length_L50',
        'mean_insertion_per_read', 'mean_insertion_per_read_per_kb',
        'insertion_per_query_kb', 'insertion_per_aligned_kb',
        'mean_deletion_per_read', 'mean_deletion_per_read_per_kb',
        'deletion_per_query_kb', 'deletion_per_aligned_kb',
        'mean_mismatch_per_read', 'mean_mismatch_per_read_per_kb',
        'mismatch_per_query_kb', 'mismatch_per_aligned_kb',
        'mean_intron_per_read', 'mean_intron_per_read_per_kb'
    ]

    def _row(a):
        N50, L50 = a.get_length_NL(50)
        return [
            a.label,
            a.get_read_count(),
            a.get_total_base(),
            a.get_total_aligned_base(),
            a.get_mean_length(),
            a.get_median_length(),
            N50, L50,
            a.get_mean_insertions(),
            a.get_mean_length_normalized_insertions() * 1000,
            a.insertions_per_base() * 1000,
            a.insertions_per_aligned_base() * 1000,
            a.get_mean_deletions(),
            a.get_mean_length_normalized_deletions() * 1000,
            a.deletions_per_base() * 1000,
            a.deletions_per_aligned_base() * 1000,
            a.get_mean_mismatches(),
            a.get_mean_length_normalized_mismatches() * 1000,
            a.mismatches_per_base() * 1000,
            a.mismatches_per_aligned_base() * 1000,
            a.get_mean_introns(),
            a.get_mean_length_normalized_introns() * 1000
        ]

    rows = [_row(a) for a in readstat_list] + [_row(readstat_sum)]
    return pd.DataFrame(rows, columns = colnames)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_table.py -v`
Expected: PASS (column-contract, value, aligned-column, and existing table tests)

- [ ] **Step 5: Commit**

```bash
git add src/lqc/report_table.py tests/test_report_table.py
git commit -m "feat: clarify read_stat normalization; add aligned-base rates"
```

---

## Task 6: Collapse `splice.txt` to four categories; add `splice_all.txt`

**Files:**
- Modify: `src/lqc/report_table.py` (`create_splice_table`, new `create_splice_all_table`)
- Modify: `src/lqc/report_html.py` (`html_add_splice_table`)
- Modify: `src/lqc/__init__.py`
- Test: `tests/test_report_table.py`, `tests/test_regressions.py`

- [ ] **Step 1: Write the failing tests**

In `tests/test_report_table.py`, replace `test_splice_table` (lines 68-77) with:

```python
def test_splice_table():
    s1 = Splice('chr1')
    s1.add_splice_pair('gt-ag')
    s2 = Splice('chr2')
    s2.add_splice_pair('gt-ag')
    total = s1 + s2
    total.label = 'Total'

    df = create_splice_table([s1, s2], total)
    assert list(df.columns) == [
        'label', 'gt-ag', 'gt-ag_pct', 'gc-ag', 'gc-ag_pct',
        'at-ac', 'at-ac_pct', 'other', 'other_pct',
    ]
    total_row = df[df['label'] == 'Total'].iloc[0]
    assert total_row['gt-ag'] == 2
    assert total_row['gt-ag_pct'] == 100.0
    assert total_row['other'] == 0
    assert total_row['other_pct'] == 0.0
```

Add a test for the full matrix alongside it:

```python
def test_splice_all_table_preserves_full_matrix():
    s = Splice('chr1')
    s.add_splice_pair_count_dict({'gt-ag': 3, 'gc-ag': 1, 'tt-tt': 4})
    total = Splice('Total')
    total.add_splice_pair_count_dict({'gt-ag': 3, 'gc-ag': 1, 'tt-tt': 4})

    df = create_splice_all_table([s], total)
    assert list(df.columns) == ['label', 'gt-ag', 'gc-ag', 'tt-tt']
    total_row = df[df['label'] == 'Total'].iloc[0]
    assert total_row['gt-ag'] == 3
    assert total_row['tt-tt'] == 4
```

Update the import block at the top of `tests/test_report_table.py` so `from lqc.report_table import (...)` reads:

```python
from lqc.report_table import (
    create_indel_summary_table,
    create_mapping_table,
    create_readstat_table,
    create_splice_all_table,
    create_splice_table,
)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_report_table.py::test_splice_table tests/test_report_table.py::test_splice_all_table_preserves_full_matrix -v`
Expected: FAIL (`create_splice_table` returns `['label', 'gt-ag']`, and `create_splice_all_table` does not exist)

- [ ] **Step 3: Rework `create_splice_table` and add `create_splice_all_table`**

In `src/lqc/report_table.py`, replace `create_splice_table` (lines 134-154) with:

```python
def create_splice_table(splice_list, splice_sum):
    """Four-category splice summary (gt-ag / gc-ag / at-ac / other, count+pct)."""

    def _row(a):
        count_dict = a.get_splice_pair_count_dict()
        gtag = count_dict.get('gt-ag', 0)
        gcag = count_dict.get('gc-ag', 0)
        atac = count_dict.get('at-ac', 0)
        other = sum(
            v for k, v in count_dict.items()
            if k not in ('gt-ag', 'gc-ag', 'at-ac')
        )
        total = gtag + gcag + atac + other
        if total == 0:
            gtagp = gcagp = atacp = otherp = 0.0
        else:
            gtagp = gtag / total * 100
            gcagp = gcag / total * 100
            atacp = atac / total * 100
            otherp = other / total * 100
        return [a.label, gtag, gtagp, gcag, gcagp,
                atac, atacp, other, otherp]

    rows = [_row(a) for a in splice_list] + [_row(splice_sum)]
    return pd.DataFrame(
        rows,
        columns = ['label', 'gt-ag', 'gt-ag_pct',
                   'gc-ag', 'gc-ag_pct',
                   'at-ac', 'at-ac_pct',
                   'other', 'other_pct']
    )


def create_splice_all_table(splice_list, splice_sum):
    """Full splice-pair matrix (one column per observed pair)."""

    sptypes = list(splice_sum.get_splice_pair_count_dict().keys())
    rows = [
        [splice_list[i].label] +
        [splice_list[i].get_splice_pair_count_dict().get(a, 0)
         for a in sptypes]
        for i in range(len(splice_list))
    ] + [
        [splice_sum.label] +
        [splice_sum.get_splice_pair_count_dict().get(a, 0)
         for a in sptypes]
    ]
    return pd.DataFrame(rows, columns = ['label', *sptypes])
```

- [ ] **Step 4: Simplify `html_add_splice_table`**

In `src/lqc/report_html.py`, replace the body of `html_add_splice_table` (lines 282-379) with:

```python
def html_add_splice_table(html_string, splice_table,
                          mean_intron_per_read):

    rowstring_list = []
    total_gtag = 0
    total_gtagp = 0.0
    total_gcag = 0
    total_gcagp = 0.0
    total_atac = 0
    total_atacp = 0.0
    total_other = 0
    total_otherp = 0.0
    for _, row in splice_table.iterrows():
        gtag = row['gt-ag']
        gtagp = row['gt-ag_pct']
        gcag = row['gc-ag']
        gcagp = row['gc-ag_pct']
        atac = row['at-ac']
        atacp = row['at-ac_pct']
        other = row['other']
        otherp = row['other_pct']
        tmprow_list = [
            f'<th scope="row">{row["label"]}</th>',
            f'<td>{gtag}</td>',
            f'<td>{gtagp:.4}</td>',
            f'<td>{gcag}</td>',
            f'<td>{gcagp:.4}</td>',
            f'<td>{atac}</td>',
            f'<td>{atacp:.4}</td>',
            f'<td>{other}</td>',
            f'<td>{otherp:.4}</td>'
        ]
        if row['label'] == "Total":
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + "</tr>"
            )
            total_gtag = gtag
            total_gtagp = gtagp
            total_gcag = gcag
            total_gcagp = gcagp
            total_atac = atac
            total_atacp = atacp
            total_other = other
            total_otherp = otherp
        else:
            rowstring_list.append(
                "<tr>" +
                '\n'.join(tmprow_list) + "</tr>"
            )

    splice_list = [
        f"<li>gt-ag: {total_gtag} ({total_gtagp:.4}%)</li>",
        f"<li>gc-ag: {total_gcag} ({total_gcagp:.4}%)</li>",
        f"<li>at-ac: {total_atac} ({total_atacp:.4}%)</li>",
        f"<li>other: {total_other} ({total_otherp:.4}%)</li>"
    ]
    new_html_string = re.sub(
        r"\{%splice_type_list%\}",
        "<ul>" + '\n'.join(splice_list) + "</ul>",
        html_string
    )
    new_html_string = re.sub(
        r"\{%intron_total_intron_number%\}",
        '{}'.format(
            total_gtag + total_gcag +
            total_atac + total_other
        ),
        new_html_string
    )
    new_html_string = re.sub(
        r"\{%intron_mean_intron_per_read%\}",
        f'{mean_intron_per_read:.4}',
        new_html_string
    )
    new_html_string = re.sub(
        r"\{%splice_table%\}",
        '\n'.join(rowstring_list),
        new_html_string
    )
    return new_html_string
```

- [ ] **Step 5: Update the empty-splice regression test**

In `tests/test_regressions.py`, replace `test_splice_html_table_with_no_introns_does_not_divide_by_zero` (lines 50-56) with:

```python
def test_splice_html_table_with_no_introns_does_not_divide_by_zero():
    s = Splice('chr1')
    total = Splice('Total')
    table = create_splice_table([s], total)
    assert list(table.columns) == [
        'label', 'gt-ag', 'gt-ag_pct', 'gc-ag', 'gc-ag_pct',
        'at-ac', 'at-ac_pct', 'other', 'other_pct',
    ]
    total_row = table.iloc[1]
    assert total_row['other_pct'] == 0.0
    html = html_add_splice_table('<html>{%splice_table%}</html>', table, 0.0)
    assert '{%splice_table%}' not in html
```

- [ ] **Step 6: Export `create_splice_all_table`**

In `src/lqc/__init__.py`, add `create_splice_all_table` to the `from lqc.report_table import (...)` block and to `__all__`.

- [ ] **Step 7: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_table.py tests/test_regressions.py -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add src/lqc/report_table.py src/lqc/report_html.py src/lqc/__init__.py tests/test_report_table.py tests/test_regressions.py
git commit -m "feat: collapse splice table to 4 categories; add splice_all table"
```

---

## Task 7: Fix mismatch type ordering and add `bin_total`

**Files:**
- Modify: `src/lqc/report_table.py`
- Test: `tests/test_report_table.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_table.py`:

```python
def test_mismatch_table_fixed_type_order_and_bin_total():
    a = Mismatch('chr1')
    a.add_mismatch('ct', 0.05)
    a.add_mismatch('ag', 0.05)
    total = Mismatch('Total')
    total.add_mismatch('ct', 0.05)
    total.add_mismatch('ag', 0.05)

    df = create_mismatch_normalized_read_location_table([a], total)
    # fixed canonical order, present types only: ag before ct
    assert list(df.columns) == ['label', 'bin', 'ag', 'ct', 'bin_total']
    row = df[df['label'] == 'chr1'].iloc[0]
    assert row['ag'] == 1
    assert row['ct'] == 1
    assert row['bin_total'] == 2
```

Update the import block at the top of `tests/test_report_table.py` to read:

```python
from lqc import Indel, Mismatch, Splice
from lqc.readstat import ReadStat
from lqc.report_table import (
    create_indel_summary_table,
    create_mapping_table,
    create_mismatch_normalized_read_location_table,
    create_readstat_table,
    create_splice_all_table,
    create_splice_table,
)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_report_table.py::test_mismatch_table_fixed_type_order_and_bin_total -v`
Expected: FAIL (`create_mismatch_normalized_read_location_table` returns columns `['label', 'bin', ...]` in insertion order and has no `bin_total`)

- [ ] **Step 3: Implement the rework**

In `src/lqc/report_table.py`, replace `create_mismatch_normalized_read_location_table` (lines 74-107) with:

```python
MISMATCH_TYPES_IN_ORDER = ['ac', 'ag', 'at', 'ca', 'cg', 'ct',
                           'ga', 'gc', 'gt', 'ta', 'tc', 'tg']


def create_mismatch_normalized_read_location_table(mismatch_list,
                                                   mismatch_sum):
    cuts = [0, 0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1]
    mis_type_bin_counts = [
        mismatch_list[i].get_location_bin_count_by_type(cuts = cuts)
        for i in range(len(mismatch_list))
    ]
    sum_type_bin_count = mismatch_sum.get_location_bin_count_by_type(cuts = cuts)
    mistypes = [
        t for t in MISMATCH_TYPES_IN_ORDER
        if t in sum_type_bin_count
    ]
    if not mistypes:
        return pd.DataFrame(columns = ['label', 'bin'])
    bins = list(
        sum_type_bin_count[mistypes[0]].keys()
    )

    data_list = [
        [mismatch_list[i].label, ibin] +
        [mis_type_bin_counts[i][c][ibin]
         for c in mistypes]
        for i in range(len(mismatch_list))
        for ibin in bins
    ]
    data_list.extend(
        [mismatch_sum.label, ibin] +
        [sum_type_bin_count[c][ibin]
         for c in mistypes]
        for ibin in bins
    )
    mistable = pd.DataFrame(
        data_list,
        columns = ['label', 'bin', *mistypes]
    )
    mistable = mistable.copy()
    mistable['bin_total'] = mistable[mistypes].sum(axis = 1)
    return mistable
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_table.py tests/test_regressions.py -v`
Expected: PASS (the empty-mismatch regression test still returns `['label', 'bin']`)

- [ ] **Step 5: Commit**

```bash
git add src/lqc/report_table.py tests/test_report_table.py
git commit -m "feat: fix mismatch type column order; add per-bin total"
```

---

## Task 8: Share palette and make indel histogram readable

**Files:**
- Modify: `src/lqc/report_figure.py`
- Test: `tests/test_report_figure.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_figure.py`:

```python
from lqc.report_figure import plot_element_total_count, plot_indel_hist_length
from lqc.indel import Indel


def test_plot_readstat_bar_multi_contig_cycles_palette():
    rs = []
    for i in range(5):
        r = ReadStat(f'chr{i}')
        r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
        rs.append(r)
    fig = plot_readstat_bar(rs, 'Read count')
    bars = fig.axes[0].patches
    colors = {p.get_facecolor() for p in bars}
    assert len(colors) > 2
    plt.close('all')


def test_plot_indel_hist_length_readable():
    indels = [Indel('chr1')]
    # one short and one long indel: the long tail must not hide the short bin
    indels[0].add_indel('a', 0.1)
    indels[0].add_indel('t' * 1000, 0.5)
    fig = plot_indel_hist_length(indels, width=5, height=4)
    ax = fig.axes[0]
    assert ax.get_yscale() == 'log'
    plt.close('all')
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_report_figure.py::test_plot_readstat_bar_multi_contig_cycles_palette tests/test_report_figure.py::test_plot_indel_hist_length_readable -v`
Expected: FAIL (first: `len(colors)` is ≤ 2; second: `ax.get_yscale()` is `'linear'`)

- [ ] **Step 3: Use the shared palette in the two bar plotters**

In `src/lqc/report_figure.py`, in `plot_readstat_bar` replace the `axes.bar(...)` call (lines 117-125) with:

```python
    colors = [
        magentas[i % len(magentas)]
        for i in range(len(read_feature))
    ]
    axes.bar(
        [read_feature[i][0]
         for i in range(len(read_feature))],
        [read_feature[i][2]
         for i in range(len(read_feature))],
        width = 0.66,
        color = colors,
        fill = True
    )
```

In `plot_element_total_count`, replace the `axes.bar(...)` call (lines 508-516) with:

```python
    colors = [
        magentas[i % len(magentas)]
        for i in range(len(count_list))
    ]
    axes.bar(
        [count_list[i][0]
         for i in range(len(count_list))],
        [count_list[i][2]
         for i in range(len(count_list))],
        width = 0.66,
        color = colors,
        fill = True
    )
```

- [ ] **Step 4: Make the indel histogram readable**

In `src/lqc/report_figure.py`, add a module-level helper right before `plot_indel_hist_length`, and rewrite that function's facet loop (lines 551-573):

```python
def _indel_length_bins(lengths):
    """1-bp bins up to 10, plus one catch-all tail bin for longer indels."""
    maxlen = max(lengths)
    edges = [i + 0.5 for i in range(0, 11)]
    if maxlen > 10:
        edges.append(maxlen + 0.5)
    return edges


def plot_indel_hist_length(indel_list,
                           width = None,
                           height = None):

    row, col = get_facet_row_col(
        len(indel_list)
    )
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, _ = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(indel_list))):
            lengths = indel_list[ai].get_lengths()
            if lengths:
                fig.axes[ai].hist(
                    lengths,
                    bins = _indel_length_bins(lengths),
                    color = "#7a0177",
                    label = indel_list[ai].label
                )
                fig.axes[ai].set_yscale('log')
                fig.axes[ai].set_title(
                    indel_list[ai].label
                )
                xlimits = fig.axes[ai].get_xlim()
                ylimits = fig.axes[ai].get_ylim()
                fig.axes[ai].text(
                    xlimits[1] * 0.98,
                    ylimits[1] * 0.98,
                    "Median:\n{}".format(int(
                        indel_list[ai].get_median_length()
                    )),
                    ha = "right", va = "top"
                )
            else:
                fig.axes[ai].set_frame_on(False)
                fig.axes[ai].set_axis_off()
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    fig.supxlabel("Element length")
    fig.supylabel("Element count")
    plt.tight_layout()
    return fig
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_figure.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/lqc/report_figure.py tests/test_report_figure.py
git commit -m "feat: share palette across bar plots; log-scale indel histogram"
```

---

## Task 9: Add the three mapping figures

**Files:**
- Modify: `src/lqc/report_figure.py`
- Modify: `src/lqc/__init__.py`
- Test: `tests/test_report_figure.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_report_figure.py`:

```python
from lqc.report_figure import (
    plot_mapping_aligned_fraction_hist,
    plot_mapping_aligned_vs_query,
    plot_mapping_mapq_hist,
)


def _one_mapped_readstat():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
               mapping_quality=60, aligned_length=180)
    return r


def test_plot_mapping_mapq_hist_returns_figure():
    fig = plot_mapping_mapq_hist([_one_mapped_readstat()])
    assert isinstance(fig, Figure)
    plt.close('all')


def test_plot_mapping_aligned_fraction_hist_returns_figure():
    fig = plot_mapping_aligned_fraction_hist([_one_mapped_readstat()])
    assert isinstance(fig, Figure)
    plt.close('all')


def test_plot_mapping_aligned_vs_query_returns_figure():
    fig = plot_mapping_aligned_vs_query([_one_mapped_readstat()])
    assert isinstance(fig, Figure)
    plt.close('all')
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_report_figure.py::test_plot_mapping_mapq_hist_returns_figure -v`
Expected: FAIL with `ImportError: cannot import name 'plot_mapping_mapq_hist'`

- [ ] **Step 3: Implement the three functions**

Append to `src/lqc/report_figure.py` (after `plot_indel_hist_location`, before `plot_mismatch_type_count`):

```python
def plot_mapping_mapq_hist(readstat_list,
                           width = None,
                           height = None):
    row, col = get_facet_row_col(len(readstat_list))
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, _ = plt.subplots(row, col, figsize = (width, height))
    for ai in range(row * col):
        if ai in list(range(len(readstat_list))):
            qualities = readstat_list[ai].get_mapping_qualities()
            maxq = max(qualities) if qualities else 60
            bins = list(range(0, max(maxq, 60) + 2))
            fig.axes[ai].hist(
                qualities,
                bins = bins,
                color = magentas[0],
                label = readstat_list[ai].label
            )
            fig.axes[ai].set_title(readstat_list[ai].label)
            xlimits = fig.axes[ai].get_xlim()
            ylimits = fig.axes[ai].get_ylim()
            fig.axes[ai].text(
                xlimits[1] * 0.98,
                ylimits[1] * 0.98,
                "Median:\n{}".format(readstat_list[ai].get_median_mapping_quality()),
                ha = "right", va = "top"
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()
    fig.supxlabel("Mapping quality")
    fig.supylabel("Read count")
    plt.tight_layout()
    return fig


def plot_mapping_aligned_fraction_hist(readstat_list,
                                       width = None,
                                       height = None):
    row, col = get_facet_row_col(len(readstat_list))
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, _ = plt.subplots(row, col, figsize = (width, height))
    for ai in range(row * col):
        if ai in list(range(len(readstat_list))):
            fractions = readstat_list[ai].get_aligned_fractions()
            fig.axes[ai].hist(
                fractions,
                bins = 50,
                range = (0, 1),
                color = magentas[0],
                label = readstat_list[ai].label
            )
            fig.axes[ai].set_title(readstat_list[ai].label)
            xlimits = fig.axes[ai].get_xlim()
            ylimits = fig.axes[ai].get_ylim()
            fig.axes[ai].text(
                xlimits[1] * 0.98,
                ylimits[1] * 0.98,
                "Median:\n{:.4}".format(
                    readstat_list[ai].get_median_aligned_fraction()
                ),
                ha = "right", va = "top"
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()
    fig.supxlabel("Aligned fraction")
    fig.supylabel("Read count")
    plt.tight_layout()
    return fig


def plot_mapping_aligned_vs_query(readstat_list,
                                  width = None,
                                  height = None):
    row, col = get_facet_row_col(len(readstat_list))
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 3
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, _ = plt.subplots(row, col, figsize = (width, height))
    for ai in range(row * col):
        if ai in list(range(len(readstat_list))):
            qlen = readstat_list[ai].get_lengths()
            alen = readstat_list[ai].get_aligned_lengths()
            if qlen:
                fig.axes[ai].hexbin(
                    qlen, alen,
                    gridsize = 40,
                    bins = 'log',
                    cmap = 'magma_r'
                )
                high = max(qlen)
                fig.axes[ai].plot(
                    [0, high], [0, high],
                    color = '#7a0177', linestyle = ':'
                )
            fig.axes[ai].set_title(readstat_list[ai].label)
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()
    fig.supxlabel("Query length")
    fig.supylabel("Aligned length")
    plt.tight_layout()
    return fig
```

- [ ] **Step 4: Export the three functions**

In `src/lqc/__init__.py`, add `plot_mapping_mapq_hist`, `plot_mapping_aligned_fraction_hist`, and `plot_mapping_aligned_vs_query` to the `from lqc.report_figure import (...)` block and to `__all__`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_figure.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/lqc/report_figure.py src/lqc/__init__.py tests/test_report_figure.py
git commit -m "feat: add mapping quality/aligned-fraction/aligned-vs-query figures"
```

---

## Task 10: Add the "Mapping" HTML section

**Files:**
- Modify: `src/lqc/report_html.py` (add `html_add_mapping`)
- Modify: `src/lqc/template/template.html`
- Modify: `src/lqc/__init__.py`
- Test: `tests/test_report_html.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_html.py`:

```python
from lqc.report_html import html_add_mapping
from lqc.report_table import create_mapping_table


def test_html_add_mapping_substitutes_placeholders():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
                mapping_quality=60, aligned_length=180)
    r2 = ReadStat('chr2')
    r2.add_read(100, insertion=0, deletion=0, mismatch=0, intron=0,
                mapping_quality=30, aligned_length=100)
    total = sum([r1, r2])
    total.label = 'Total'
    table = create_mapping_table([r1, r2], total)

    html = html_add_mapping(get_html_template(), table)
    assert '{%mapping_aligned_fraction_mean%}' not in html
    assert '{%mapping_aligned_fraction_median%}' not in html
    assert '{%mapping_mapq_mean%}' not in html
    assert '{%mapping_mapq_median%}' not in html
    assert '{%mapping_table%}' not in html
    assert '<tr class="table-secondary">' in html
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_report_html.py::test_html_add_mapping_substitutes_placeholders -v`
Expected: FAIL with `ImportError: cannot import name 'html_add_mapping'`

- [ ] **Step 3: Implement `html_add_mapping`**

Append to `src/lqc/report_html.py` (after `html_add_splice_table`, before the `########` separator at line 382):

```python
def html_add_mapping(html_string, mapping_table):
    rowstring_list = []
    total_aligned_fraction_mean = 0.0
    total_aligned_fraction_median = 0.0
    total_mapq_mean = 0.0
    total_mapq_median = 0.0
    for _, row in mapping_table.iterrows():
        tmprow_list = [
            '<th scope="row">{}</th>'.format(row['label']),
            '<td>{}</td>'.format(row['read_count']),
            '<td>{}</td>'.format(row['query_base']),
            '<td>{}</td>'.format(row['aligned_base']),
            '<td>{:.4}</td>'.format(float(row['aligned_fraction_mean'])),
            '<td>{:.4}</td>'.format(float(row['aligned_fraction_median'])),
            '<td>{:.4}</td>'.format(float(row['mapq_mean'])),
            '<td>{:.4}</td>'.format(float(row['mapq_median'])),
            '<td>{}</td>'.format(row['reads_aligned_fraction_lt_0.9']),
            '<td>{}</td>'.format(row['reads_fully_aligned'])
        ]
        if row['label'] == 'Total':
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + '</tr>'
            )
            total_aligned_fraction_mean = float(row['aligned_fraction_mean'])
            total_aligned_fraction_median = float(row['aligned_fraction_median'])
            total_mapq_mean = float(row['mapq_mean'])
            total_mapq_median = float(row['mapq_median'])
        else:
            rowstring_list.append('<tr>' + '\n'.join(tmprow_list) + '</tr>')

    new_html_string = re.sub(
        r"\{%mapping_aligned_fraction_mean%\}",
        f'{total_aligned_fraction_mean:.4}',
        html_string
    )
    new_html_string = re.sub(
        r"\{%mapping_aligned_fraction_median%\}",
        f'{total_aligned_fraction_median:.4}',
        new_html_string
    )
    new_html_string = re.sub(
        r"\{%mapping_mapq_mean%\}",
        f'{total_mapq_mean:.4}',
        new_html_string
    )
    new_html_string = re.sub(
        r"\{%mapping_mapq_median%\}",
        f'{total_mapq_median:.4}',
        new_html_string
    )
    new_html_string = re.sub(
        r"\{%mapping_table%\}",
        '\n'.join(rowstring_list),
        new_html_string
    )
    return new_html_string
```

- [ ] **Step 4: Add the nav entry**

In `src/lqc/template/template.html`, insert this `<li>` immediately after the "Read" nav item (after the `</li>` that closes the `#collapseReadstat` toggle, around lines 31-35):

```html
            <li class="nav-item text-center">
              <div class="text-light" type="button" data-bs-toggle="collapse" data-bs-target="#collapseMapping" aria-expanded="false" aria-controls="collapseMapping" data-bs-parent="#toggleGroup">
                <p class="fs-5 fw-bold">Mapping</p>
              </div>
            </li>
```

- [ ] **Step 5: Add the accordion section**

In `src/lqc/template/template.html`, insert the whole "Mapping" section between the closing `</div>` of `collapseReadstat` (the `</div>` right before the blank line that precedes `collapseInsertion`, around line 268) and the `collapseInsertion` opening `<div>` (line 271):

```html
          <div class="accordion-collapse collapse p-3" id="collapseMapping" data-bs-parent="#toggleGroup">
            <h1 id="h_mapping">Mapping</h1>
            <div class="d-flex flex-column gap-3">
              <div class="card mb-3">
                <div class="row g-0">
                  <div class="col-md-5">
                    <div class="card-body">
                      <h3 class="card-title">Mapping quality</h3>
                      <p class="card-text">The mean mapping quality is
                      {%mapping_mapq_mean%} and the median is
                      {%mapping_mapq_median%}.</p>
                    </div>
                  </div>
                  <div class="col-md-7">
                    <img class="card-img" src="fig/mapping_hist_mapq.Total.png"/>
                  </div>
                </div>
              </div>
              <div class="card">
                <div class="card-body">
                  <h3 class="card-title">Mapping quality by contig</h3>
                </div>
                <img class="card-img-bottom" src="fig/mapping_hist_mapq.png" />
              </div>
              <div class="card mb-3">
                <div class="row g-0">
                  <div class="col-md-5">
                    <div class="card-body">
                      <h3 class="card-title">Aligned fraction</h3>
                      <p class="card-text">The fraction of each read that is
                      aligned (aligned length / read length). Mean
                      {%mapping_aligned_fraction_mean%}, median
                      {%mapping_aligned_fraction_median%}.</p>
                    </div>
                  </div>
                  <div class="col-md-7">
                    <img class="card-img" src="fig/mapping_hist_aligned_fraction.Total.png"/>
                  </div>
                </div>
              </div>
              <div class="card">
                <div class="card-body">
                  <h3 class="card-title">Aligned fraction by contig</h3>
                </div>
                <img class="card-img-bottom" src="fig/mapping_hist_aligned_fraction.png" />
              </div>
              <div class="card mb-3">
                <div class="row g-0">
                  <div class="col-md-5">
                    <div class="card-body">
                      <h3 class="card-title">Aligned vs query length</h3>
                      <p class="card-text">Each read's aligned query length
                      against its full query length, with the y=x fully-aligned
                      reference line.</p>
                    </div>
                  </div>
                  <div class="col-md-7">
                    <img class="card-img" src="fig/mapping_scatter_aligned_vs_query.Total.png"/>
                  </div>
                </div>
              </div>
              <div class="card">
                <div class="card-body">
                  <h3 class="card-title">Aligned vs query length by contig</h3>
                </div>
                <img class="card-img-bottom" src="fig/mapping_scatter_aligned_vs_query.png" />
              </div>
              <div class="card">
                <div class="card-header">Statistic Summary</div>
                <div class="card-body">
                  <table class="table table-striped text-center">
                    <thead>
                      <tr>
                        <th scope="col">Contig</th>
                        <th scope="col">Read count</th>
                        <th scope="col">Query base</th>
                        <th scope="col">Aligned base</th>
                        <th scope="col">Aligned fraction mean</th>
                        <th scope="col">Aligned fraction median</th>
                        <th scope="col">MAPQ mean</th>
                        <th scope="col">MAPQ median</th>
                        <th scope="col">Reads &lt;0.9 aligned</th>
                        <th scope="col">Fully aligned</th>
                      </tr>
                    </thead>
                    <tbody>
                      {%mapping_table%}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>
          </div>

```

- [ ] **Step 6: Export `html_add_mapping`**

In `src/lqc/__init__.py`, add `html_add_mapping` to the `from lqc.report_html import (...)` block and to `__all__`.

- [ ] **Step 7: Run tests to verify they pass**

Run: `uv run pytest tests/test_report_html.py -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add src/lqc/report_html.py src/lqc/template/template.html src/lqc/__init__.py tests/test_report_html.py
git commit -m "feat: add Mapping section to HTML report"
```

---

## Task 11: Wire everything into the CLI

**Files:**
- Modify: `src/lqc/cli.py`

- [ ] **Step 1: Extend the imports and `o_files`**

In `src/lqc/cli.py`, add to the `from lqc import (...)` block: `create_mapping_table`, `create_splice_all_table`, `html_add_mapping`, `plot_mapping_aligned_fraction_hist`, `plot_mapping_aligned_vs_query`, `plot_mapping_mapq_hist`.

In the `o_files` dict (around line 236-350), add these entries (next to the other table/fig keys):

```python
        't_mapping': os.path.join(
            o_dirs['table'], 'mapping.txt'
        ),
        't_splice_all': os.path.join(
            o_dirs['table'], 'splice_all.txt'
        ),
        'f_mapping_mapq_hist': os.path.join(
            o_dirs['fig'], 'mapping_hist_mapq'
        ),
        'f_mapping_aligned_fraction_hist': os.path.join(
            o_dirs['fig'], 'mapping_hist_aligned_fraction'
        ),
        'f_mapping_aligned_vs_query': os.path.join(
            o_dirs['fig'], 'mapping_scatter_aligned_vs_query'
        ),
```

- [ ] **Step 2: Build the new tables**

In `src/lqc/cli.py`, right after `t_splice = create_splice_table(...)` (line 452), add:

```python
    t_mapping = create_mapping_table(l_readstat, sreadstat)
    t_splice_all = create_splice_all_table(l_splice, ssplice)
```

- [ ] **Step 3: Write the new tables**

In `src/lqc/cli.py`, after the `t_splice.to_csv(...)` call (lines 521-524), add:

```python
    t_mapping.to_csv(
        o_files['t_mapping'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
    t_splice_all.to_csv(
        o_files['t_splice_all'],
        sep = '\t', index = False,
        float_format = FLOAT_FORMAT
    )
```

- [ ] **Step 4: Generate the mapping figures**

In `src/lqc/cli.py`, after the `splice_type` figure block (lines 704-714) and before the `message = 'Output figures finished.'` line, add:

```python
    # mapping metrics
    message = 'Output figures: mapping metrics.'
    logger.debug(message)
    filelabel = 'f_mapping_mapq_hist'
    generate_multiple_figs(
        plot_mapping_mapq_hist,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )
    filelabel = 'f_mapping_aligned_fraction_hist'
    generate_multiple_figs(
        plot_mapping_aligned_fraction_hist,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )
    filelabel = 'f_mapping_aligned_vs_query'
    generate_multiple_figs(
        plot_mapping_aligned_vs_query,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )
```

- [ ] **Step 5: Add the HTML step and JSON entry**

In `src/lqc/cli.py`, after `new_html_string = html_add_readstat_table(html_string, t_readstat)` (line 728-730), add:

```python
    new_html_string = html_add_mapping(
        new_html_string, t_mapping
    )
```

In the `html_add_data` call (lines 770-779), add `'mapping': t_mapping` to the tables dict:

```python
        {
            'readstat': t_readstat,
            'mapping': t_mapping,
            'insertion': t_insertion,
            'deletion': t_deletion,
            'mismatch': t_mismatch,
            'splice': t_splice,
        }
```

- [ ] **Step 6: Run the CLI integration tests**

Run: `uv run pytest tests/test_cli.py -v`
Expected: PASS (`test_main_report_is_self_contained`, `test_tables_identical_across_thread_counts`, etc.)

- [ ] **Step 7: Commit**

```bash
git add src/lqc/cli.py
git commit -m "feat: wire mapping metrics and splice_all into CLI output"
```

---

## Task 12: Update docs and verify end-to-end

**Files:**
- Modify: `docs/reporting-and-output.md`

- [ ] **Step 1: Update `docs/reporting-and-output.md`**

In the "Summary tables (columns)" section, update the read_stat bullet to the new column list (add `aligned_base`, rename the three `*_per_kb` to `*_per_query_kb`, add `*_per_aligned_kb`), add bullets for `mapping.txt` and `splice_all.txt`, and in the output-layout code block add `mapping.txt`, `splice_all.txt`, and the three `mapping_*` figure names. In the "Figures" section, list the three new `plot_mapping_*` figures.

- [ ] **Step 2: Run lint and the full unit suite**

Run: `uv run ruff check`
Expected: no errors.

Run: `uv run pytest`
Expected: all tests PASS.

- [ ] **Step 3: Smoke-run against the real BAM**

Run: `bash tmp/analyze_chr22.sh`
Expected: completes; `tmp/out-src/table/mapping.txt` exists, its `Total` row has `aligned_fraction_median` ≈ 0.8404, `mapq_median` = 13, `reads_aligned_fraction_lt_0.9` = 194381; `splice.txt` has 9 columns; `splice_all.txt` exists; no `{%` remains in `LQC_report.html`.

*(Note: the earlier `mapq_median = 1` expectation was corrected to `13` during
execution — independently re-derived from the BAM's MAPQ distribution median.)*

Run to confirm the exact figures are emitted:

Run: `ls tmp/out-src/fig/mapping_hist_mapq*.png tmp/out-src/fig/mapping_hist_aligned_fraction*.png tmp/out-src/fig/mapping_scatter_aligned_vs_query*.png`
Expected: 6 `.png` (3 metrics × Total/contig) plus PDFs.

- [ ] **Step 4: Commit**

```bash
git add docs/reporting-and-output.md
git commit -m "docs: document mapping metrics, splice_all, and column changes"
```

---

## Self-Review (run by the plan author)

- **Spec coverage:** Task 1-2 = data collection; Task 3 = rounding (#4); Task 4 = mapping table; Task 5 = #2 normalization; Task 6 = #1 splice; Task 7 = #3 mismatch; Task 8 = #5 palette + #6 indel hist; Task 9 = new figures; Task 10 = HTML; Task 11 = wiring/exports; Task 12 = docs + verify. All ten spec points are covered.
- **Type/name consistency:** `mapping_quality`/`aligned_length` field names match across `ReadRecord`, `record_from_read`, `add_read`, `process_record`, and all tests. `FLOAT_FORMAT` is imported in `cli.py` only where used. The `mapping_*` token names in `template.html` match exactly the five tokens `html_add_mapping` substitutes.
- **Executing agents note:** Tasks are ordered so the suite stays green at every commit boundary (each new getter/function is added before its first use, and each table/function rework updates its test in the same task).