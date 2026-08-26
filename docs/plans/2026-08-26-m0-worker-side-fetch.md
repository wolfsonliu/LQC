# M0: Worker-Side Fetch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace LQC's global read materialization + IPC pickling with balanced, worker-side coordinate-window fetching so a whole-genome BAM (24 contigs, ~1.57 M reads) completes within a few GB of RAM.

**Architecture:** `plan_tasks` reads the index and emits `(index, contig, start, end)` coordinate windows; each `mp.Pool` worker runs `stat_region`, which opens the BAM itself, fetches its window, assigns each read to the window containing its `reference_start`, and delegates accumulation to the existing `stat_records`. `reduce_blocks_to_contigs` is unchanged, and window indices preserve current per-contig traversal order so `read.cs`/tables stay deterministic and byte-identical.

**Tech Stack:** Python 3.11+, pysam, `multiprocessing` (unchanged stack).

---

## File Structure

- **Modify `src/lqc/stat.py`** — the parallel-processing core.
  - `record_from_read` (lines 48–76): use `read.query_length` instead of `len(read.query_sequence)`.
  - Add `plan_tasks(bam_file, contigs, target_reads=10000)` (after `prefetch_records`).
  - Add `stat_region(task, bam_file, genome_file, method, cs_dir=None)` (after `stat_records`).
  - `prefetch_records`, `stat_records`, `reduce_blocks_to_contigs`, `stat_element_from_bam_by_contig` stay untouched (still used by tests / library path / the profile harness).
- **Modify `src/lqc/cli.py`** — switch the orchestration hot path.
  - Import `plan_tasks, stat_region` instead of `prefetch_records, stat_records`.
  - In `main()`, replace the `prefetch_records` + `_chunk` + `mp.Pool.map(stat_records)` block with `plan_tasks` + `mp.Pool.map(stat_region)`.
  - Delete the now-unused `_chunk` function.
- **Modify `tests/conftest.py`** — add a `cs_bam_boundary` fixture (two reads, one spanning the coordinate-window midpoint).
- **Modify `tests/test_stat.py`** — new tests for `plan_tasks` and `stat_region`; update the `record_from_read` fake.

---

## Task 0: Capture baseline outputs (byte-parity reference)

This must run on the **current HEAD before any code change** so the end-to-end
diff in Task 5 has a reference. (`tmp/` is gitignored; `tmp/full_bam/` is a
read-only data mount, so write outputs under `tmp/baseline/`.)

- [ ] **Step 1: Run the current CLI on chr22**

```bash
mkdir -p tmp/baseline
.venv/bin/lqc -b tmp/full_bam/ENCFF417VHJ.sorted.bam -o tmp/baseline/chr22 -c chr22 -t 1
```

Expected: exit 0; `tmp/baseline/chr22/table/`, `fig/`, `LQC_report.html` written.

- [ ] **Step 2: Record the fixture reference `read.cs` fields**

```bash
.venv/bin/lqc -b tests/data/cs_test.bam -o tmp/baseline/fixture -c chr1 -t 1 --output-cs 2>/dev/null || true
```

(If `tests/data/cs_test.bam` does not exist, skip — the tiny fixture used by the
unit tests is built in `tests/conftest.py` at test time; Task 3's equivalence
tests already pin the stat content).

---

## Task 1: `record_from_read` reads `query_length`, not `len(query_sequence)`

**Files:**
- Modify: `src/lqc/stat.py:48-76`
- Test: `tests/test_stat.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_stat.py`:

```python
def test_record_from_read_cs_uses_query_length_not_sequence():
    from lqc.stat import record_from_read

    class FakeRead:
        is_reverse = False
        reference_start = 0
        mapping_quality = 60
        query_alignment_length = 10
        query_name = 'read1'
        query_length = 10

        @property
        def query_sequence(self):
            raise AssertionError(
                'cs path must use query_length, not materialize query_sequence'
            )

        def get_tag(self, tag):
            assert tag == 'cs'
            return ':10'

    record = record_from_read(FakeRead(), 'chr1', 'cs')
    assert record.query_length == 10
```

- [ ] **Step 2: Run it to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_stat.py::test_record_from_read_cs_uses_query_length_not_sequence -v`
Expected: FAIL — `AssertionError: cs path must use query_length, not materialize query_sequence` (the old code calls `len(read.query_sequence)`).

- [ ] **Step 3: Update the existing MD fake so it keeps working after the fix**

In `tests/test_stat.py`, inside `test_record_from_read_md_captures_mapping_fields`, add `query_length = 10` to `FakeRead` and assert it:

```python
    class FakeRead:
        is_reverse = False
        reference_start = 5
        query_sequence = 'ACGTACGTAC'
        query_length = 10          # NEW: the md path also reads this
        mapping_quality = 42
        query_alignment_length = 10
        query_name = 'read1'
        cigarstring = '10M'
        reference_end = 15

        def get_tag(self, tag):
            assert tag == 'MD'
            return '10'

    record = record_from_read(FakeRead(), 'chr1', 'md')
    assert record.mapping_quality == 42
    assert record.aligned_length == 10
    assert record.query_length == 10   # NEW
```

- [ ] **Step 4: Implement the fix**

In `src/lqc/stat.py`, change **both** occurrences of
`query_length=len(read.query_sequence),` (lines ~56 and ~70) to:

```python
            query_length=read.query_length,
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_stat.py -v`
Expected: PASS (both the new test and the updated MD test).

- [ ] **Step 6: Commit**

```bash
git add src/lqc/stat.py tests/test_stat.py
git commit -m "perf: use read.query_length instead of materializing query_sequence"
```

---

## Task 2: `plan_tasks` — balanced coordinate windows

**Files:**
- Modify: `src/lqc/stat.py`
- Test: `tests/test_stat.py`

- [ ] **Step 1: Add `plan_tasks` to the test import**

In `tests/test_stat.py`, change the top import to:

```python
from lqc.stat import (
    plan_tasks,
    prefetch_records,
    reduce_blocks_to_contigs,
    stat_element_from_bam_by_contig,
    stat_region,
    stat_records,
)
```

(Importing `plan_tasks`/`stat_region` before they exist will make collection fail — expected for TDD.)

- [ ] **Step 2: Write the failing test**

Append to `tests/test_stat.py`:

```python
def test_plan_tasks_splits_by_target_reads(cs_bam):
    # cs_bam: chr1 (LN=1,000,000) with 2 mapped reads.
    # target_reads=1  -> ceil(2/1)=2 windows of 500,000 bp each.
    assert plan_tasks(cs_bam, ['chr1'], target_reads=1) == [
        (0, 'chr1', 0, 500000),
        (1, 'chr1', 500000, 1000000),
    ]
    # target_reads >= mapped -> a single full-contig window.
    assert plan_tasks(cs_bam, ['chr1'], target_reads=100) == [
        (0, 'chr1', 0, 1000000),
    ]
    # A contig with zero mapped reads is skipped.
    assert plan_tasks(cs_bam, ['chr1', 'chr999'], target_reads=1) == [
        (0, 'chr1', 0, 500000),
        (1, 'chr1', 500000, 1000000),
    ]
```

- [ ] **Step 3: Run it to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_stat.py::test_plan_tasks_splits_by_target_reads -v`
Expected: FAIL — `ImportError: cannot import name 'plan_tasks'`.

- [ ] **Step 4: Implement `plan_tasks`**

In `src/lqc/stat.py`, add after `prefetch_records` (after line 155):

```python
def plan_tasks(bam_file, contigs, target_reads=10000):
    """Return ``[(index, contig, start, end), ...]`` coordinate windows.

    Splits each contig's coordinate span into windows of roughly
    ``target_reads`` mapped reads so a large contig can be processed by
    several workers. Windows are emitted in requested-contig order, then
    ascending coordinate order, so concatenating per-window results preserves
    the current single-contig traversal order.
    """
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    try:
        stats = {
            stat.contig: stat.mapped
            for stat in bam.get_index_statistics()
        }
        tasks = []
        index = 0
        for contig in contigs:
            mapped = stats.get(contig, 0)
            if mapped == 0:
                continue
            length = bam.get_reference_length(contig)
            n_windows = max(1, (mapped + target_reads - 1) // target_reads)
            for w in range(n_windows):
                start = (w * length) // n_windows
                end = ((w + 1) * length) // n_windows
                tasks.append((index, contig, start, end))
                index += 1
    finally:
        bam.close()
    return tasks
```

**Note:** `plan_tasks`/`stat_region` rely on a sorted, indexed BAM
(`get_index_statistics()` + windowed `fetch`), consistent with AGENTS.md
constraint #4. The single-contig library path `stat_element_from_bam_by_contig`
stays unchanged as the index-free/linear option.

- [ ] **Step 5: Run the test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_stat.py::test_plan_tasks_splits_by_target_reads -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/lqc/stat.py tests/test_stat.py
git commit -m "feat: plan balanced per-contig coordinate-window tasks"
```

---

## Task 3: `stat_region` — worker-side fetch with boundary de-duplication

**Files:**
- Modify: `src/lqc/stat.py`
- Modify: `tests/conftest.py` (add `cs_bam_boundary` fixture)
- Test: `tests/test_stat.py`

- [ ] **Step 1: Add the boundary fixture**

In `tests/conftest.py`, append:

```python
@pytest.fixture()
def cs_bam_boundary(tmp_path):
    """chr1 with two reads: one spanning the coordinate-window midpoint.

    With target_reads=1 the contig splits into windows [0,500000) and
    [500000,1000000). 'span' starts at 499990 with a 20M cigar, so it
    overlaps the boundary; it must be assigned to the FIRST window only.
    """
    bam_path = tmp_path / 'boundary.cs.bam'
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': 'chr1', 'LN': 1000000}],
    }
    with pysam.AlignmentFile(bam_path, 'wb', header=header) as out:
        span = pysam.AlignedSegment()
        span.query_name = 'span'
        span.flag = 0
        span.reference_id = 0
        span.reference_start = 499990
        span.mapping_quality = 60
        span.cigar = [(0, 20)]
        span.query_sequence = 'ACGT' * 5
        span.query_qualities = pysam.qualitystring_to_array('I' * 20)
        span.set_tag('cs', ':20')
        out.write(span)

        far = pysam.AlignedSegment()
        far.query_name = 'far'
        far.flag = 0
        far.reference_id = 0
        far.reference_start = 900000
        far.mapping_quality = 60
        far.cigar = [(0, 10)]
        far.query_sequence = 'ACGTACGTAC'
        far.query_qualities = pysam.qualitystring_to_array('I' * 10)
        far.set_tag('cs', ':10')
        out.write(far)
    pysam.index(str(bam_path))
    return str(bam_path)
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_stat.py`:

```python
def test_stat_region_matches_serial(cs_bam):
    # Two windows over chr1 (both reads fall in the first window; the second
    # is empty). The reduced result must equal the serial single-contig path.
    tasks = plan_tasks(cs_bam, ['chr1'], target_reads=1)
    blocks = [
        stat_region(task, cs_bam, None, 'cs')
        for task in tasks
    ]
    reduced = reduce_blocks_to_contigs(blocks, ['chr1'])
    serial = stat_element_from_bam_by_contig(cs_bam, None, 'chr1', 'cs')
    _assert_stat_tuples_equal(reduced[0], serial)
    assert reduced[0][0].get_read_count() == 2


def test_stat_region_assigns_boundary_read_once(cs_bam_boundary):
    tasks = plan_tasks(cs_bam_boundary, ['chr1'], target_reads=1)
    blocks = [
        stat_region(task, cs_bam_boundary, None, 'cs')
        for task in tasks
    ]
    # 'span' (reference_start=499990, 20M) crosses the 500000 boundary; it must
    # be counted in the first window and skipped in the second.
    assert [b.readstat.get_read_count() for b in blocks] == [1, 1]
    reduced = reduce_blocks_to_contigs(blocks, ['chr1'])
    serial = stat_element_from_bam_by_contig(cs_bam_boundary, None, 'chr1', 'cs')
    _assert_stat_tuples_equal(reduced[0], serial)
```

- [ ] **Step 3: Run them to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_stat.py::test_stat_region_matches_serial tests/test_stat.py::test_stat_region_assigns_boundary_read_once -v`
Expected: FAIL — `ImportError: cannot import name 'stat_region'`.

- [ ] **Step 4: Implement `stat_region`**

In `src/lqc/stat.py`, add after `stat_records` (after line 208):

```python
def stat_region(task, bam_file, genome_file, method, cs_dir=None):
    """Process one ``(index, contig, start, end)`` coordinate-window task.

    Opens the BAM itself, fetches the window, and assigns each read to the
    window containing its ``reference_start`` (a read that only overlaps the
    window's left edge belongs to the previous window and is skipped). The
    selected records are then accumulated by ``stat_records``.
    """
    index, contig, start, end = task
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    records = []
    for read in bam.fetch(contig, start, end):
        if read.reference_start < start:
            continue
        records.append(record_from_read(read, contig, method))
    bam.close()
    return stat_records(
        (index, contig, records), genome_file, method, cs_dir
    )
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_stat.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/lqc/stat.py tests/conftest.py tests/test_stat.py
git commit -m "feat: add stat_region worker with coordinate-window fetch"
```

---

## Task 4: Wire `cli.py` to the new path and delete `_chunk`

**Files:**
- Modify: `src/lqc/cli.py`

No new behavior to TDD-pin here — this is a pure internal refactor whose
correctness is pinned by the existing integration tests (`test_main_*`,
`test_tables_identical_across_thread_counts`,
`test_output_cs_identical_across_thread_counts`) which now run through the new
path. Verify with those tests plus the full suite.

- [ ] **Step 1: Change the import**

In `src/lqc/cli.py`, change line 56 from:

```python
from lqc.stat import prefetch_records, reduce_blocks_to_contigs, stat_records
```

to:

```python
from lqc.stat import plan_tasks, reduce_blocks_to_contigs, stat_region
```

- [ ] **Step 2: Delete the now-unused `_chunk` function**

Remove lines 82–87 (the whole `def _chunk(...)` block).

- [ ] **Step 3: Replace the run-jobs block in `main()`**

Replace from `per_contig = prefetch_records(...)` through the closing of the
`with mp.Pool(...)` block (lines 410–429) with:

```python
    tasks = plan_tasks(args['bam_file'], contigs)

    with mp.Pool(args['thread']) as p:
        block_results = p.map(
            partial(
                stat_region,
                bam_file = args['bam_file'],
                genome_file = args['genome_fasta'],
                method = bam_type,
                cs_dir = o_dirs['base'] if args['output_cs'] else None
            ),
            tasks
        )
```

(The surrounding `logger.info('Element statistic process starts.')` and the
`result = reduce_blocks_to_contigs(block_results, contigs)` line stay as-is.)

- [ ] **Step 4: Run lint and the full suite**

```bash
.venv/bin/python -m ruff check
.venv/bin/python -m pytest -q
```

Expected: `All checks passed!` and `113 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/cli.py
git commit -m "feat: worker-side fetch in CLI (plan_tasks + stat_region, drop _chunk)"
```

---

## Task 5: Byte-parity and performance verification on real data

**Files:** none (verification only). Outputs under `tmp/` (gitignored).

- [ ] **Step 1: chr22 byte-parity vs. Task 0 baseline**

```bash
mkdir -p tmp/verify
.venv/bin/lqc -b tmp/full_bam/ENCFF417VHJ.sorted.bam -o tmp/verify/chr22 -c chr22 -t 1
diff -r tmp/baseline/chr22 tmp/verify/chr22 && echo "BYTE-PARITY OK"
```

Expected: `BYTE-PARITY OK` (tables, figures, and `LQC_report.html` identical).

- [ ] **Step 2: thread determinism on chr22**

```bash
.venv/bin/lqc -b tmp/full_bam/ENCFF417VHJ.sorted.bam -o tmp/verify/chr22_t4 -c chr22 -t 4
diff -r tmp/verify/chr22 tmp/verify/chr22_t4 && echo "THREAD-DETERMINISM OK"
```

Expected: `THREAD-DETERMINISM OK`.

- [ ] **Step 3: full-genome run with wall-time + peak RSS**

```bash
bash tmp/run_lqc_full.sh tmp/full_bam/ENCFF417VHJ.sorted.bam tmp/verify/full 4
```

Expected: exits 0 and completes (previously it OOM'd >7 GB without finishing).
Record the `tmp/verify/full.time.txt` "Maximum resident set size" — target
**≤ ~2 GB** (see spec §6).

- [ ] **Step 4: Record the result in this plan's comments (no code change)**

Note the new wall time and peak RSS next to the baseline (15.62 s / 221 MB for
chr22; full run previously OOM'd). If byte-parity fails, stop and investigate
before proceeding to M1.

---

## Self-Review (run before handoff)

- **Spec coverage:** §4.1 of the design spec — `plan_tasks` (Task 2),
  `stat_region` worker-side fetch + boundary dedup (Task 3), CLI rewiring +
  `_chunk` removal (Task 4), `query_length` micro-fix (Task 1),
  `read.cs` ordering preserved via window indices (Task 3/4), byte-parity +
  perf measurement (Task 5). Covered.
- **Placeholder scan:** no TBD/TODO; every code step shows complete code and
  exact commands with expected output.
- **Type consistency:** `plan_tasks` returns `list[tuple[int, str, int, int]]`
  and `stat_region(task, bam_file, genome_file, method, cs_dir=None)` —
  consistent across Tasks 2, 3, and 4; the test calls `stat_region(task, cs_bam, None, 'cs')` match.