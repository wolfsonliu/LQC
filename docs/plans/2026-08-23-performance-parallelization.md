# Performance & Parallelization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `-t/--thread` actually speed up LQC on any BAM (single- or multi-contig) and remove the redundant second full parse behind `--output-cs`, with byte-identical outputs.

**Architecture:** Refactor `stat.py` so a read is reduced to a tiny `ReadRecord` and accumulated by a pure `process_record`; the CLI prefetches all records, chunks them into `threads` blocks, `Pool.map`s a `stat_records` worker, and reduces blocks back to per-contig tuples with the existing `sum()`/`__add__` machinery. `--output-cs` is emitted by the same workers into ordered per-chunk temp files and merged once.

**Tech Stack:** Python 3.11+, pysam, `multiprocessing`, matplotlib (Agg), existing `NamedTuple`/`typing`, uv + pytest + ruff.

**Design spec:** `docs/specs/2026-08-23-performance-parallelization-design.md`

---

## File Structure

- **Modify `src/lqc/cs.py`** — intron donor/acceptor via slicing (no regex); a shared `_get_coordinate_list`; new `get_indel_mismatches()` grouped getter; delete dead `_get_read_length`.
- **Modify `src/lqc/utils.py`** — `convert_complement` via `str.translate`.
- **Modify `src/lqc/stat.py`** — introduce `ReadRecord`, `record_from_read`, `_build_cs`, `process_record`, `prefetch_records`, `stat_records`, `reduce_blocks_to_contigs`; re-implement `stat_element_from_bam_by_contig` on top of them (same public signature).
- **Modify `src/lqc/cli.py`** — prefetch → chunk → `Pool.map` → reduce; `-t` default `min(cpu_count(), 4)`; single-pass `--output-cs` merge; remove `stat_bam`.
- **Modify `tests/test_cs.py`** — coverage for `get_indel_mismatches`.
- **Modify `tests/test_stat.py`** — parallel-vs-serial equivalence test.
- **Modify `tests/test_cli.py`** — determinism tests across thread counts.

---

## Task 1: CS micro-optimizations (`cs.py`, `utils.py`)

**Files:**
- Modify: `src/lqc/cs.py:375-382`, `:481-494`, `:591-619`
- Modify: `src/lqc/utils.py:6-14`
- Test: `tests/test_cs.py`

- [ ] **Step 1: Write the failing test for the grouped getter**

Append to `tests/test_cs.py` (after `test_deletion_counts`):

```python
def test_get_indel_mismatches_groups(cs_obj):
    insertions, deletions, mismatches = cs_obj.get_indel_mismatches(
        coordinate='normalized_read'
    )
    assert [a[2] for a in insertions] == ['+']
    assert [a[2] for a in deletions] == ['-']
    assert [a[2] for a in mismatches] == ['*']
    assert insertions[0][3] == 'tt'
    assert deletions[0][3] == 'ccc'
    assert mismatches[0][3] == 'ag'
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_cs.py::test_get_indel_mismatches_groups -v`
Expected: FAIL with `AttributeError: 'CS' object has no attribute 'get_indel_mismatches'`

- [ ] **Step 3: Implement the grouped getter and shared coordinate list**

In `src/lqc/cs.py`, replace the `_get_marks` method body (lines 591-619) with:

```python
    def _get_coordinate_list(self, coordinate):
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        if coordinate == 'contig':
            return self.get_contig_position()
        elif coordinate == 'relative_contig':
            return self.get_relative_position()
        elif coordinate == 'read':
            return self.get_read_location()
        else:
            return self.get_normalized_read_location()

    def _get_marks(self, mark, coordinate):
        assert mark in ['~', '*', ':', '+', '-'], \
            "Mark should be in [~, *, :, +, -]."
        return [
            a for a in self._get_coordinate_list(coordinate)
            if a[2] == mark
        ]

    def get_indel_mismatches(self, coordinate='normalized_read'):
        """Return (insertions, deletions, mismatches) in one pass over the
        coordinate list, avoiding three separate full-list scans."""
        insertions = []
        deletions = []
        mismatches = []
        for a in self._get_coordinate_list(coordinate):
            mark = a[2]
            if mark == '+':
                insertions.append(a)
            elif mark == '-':
                deletions.append(a)
            elif mark == '*':
                mismatches.append(a)
        return insertions, deletions, mismatches
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_cs.py::test_get_indel_mismatches_groups -v`
Expected: PASS

- [ ] **Step 5: Replace intron regex with slicing**

In `src/lqc/cs.py`, replace lines 375-382 inside `_compute_counts`:

```python
                splice_pair_count[
                    re.sub('[0-9]+', '-', value)
                ] += 1
                splice_sites = re.sub(
                    '[0-9]+', ' ', value
                ).strip().split()
                splice_l_count[splice_sites[0]] += 1
                splice_r_count[splice_sites[-1]] += 1
```

with:

```python
                donor = value[:2]
                acceptor = value[-2:]
                splice_pair_count[f'{donor}-{acceptor}'] += 1
                splice_l_count[donor] += 1
                splice_r_count[acceptor] += 1
```

- [ ] **Step 6: Delete the dead `_get_read_length`**

Delete the whole `_get_read_length` method (lines 481-494). Leave `get_read_length` (line 496-497) untouched.

- [ ] **Step 7: Rewrite `convert_complement` with `str.translate`**

In `src/lqc/utils.py`, replace lines 6-10:

```python
def convert_complement(string):
    ntpair = {'a': 't', 'c': 'g',
              'g': 'c', 't': 'a',
              'n': 'n', '-': '-'}
    return ''.join([ntpair[c] for c in string])
```

with:

```python
_COMPLEMENT_TABLE = str.maketrans(
    {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', '-': '-'}
)


def convert_complement(string):
    return string.translate(_COMPLEMENT_TABLE)
```

(Note: `str.translate` leaves characters outside `{a,c,g,t,n,-}` unchanged instead of raising `KeyError`; that is intentional and safe — valid cs/MD values only use those six lowercase characters. All existing tests and real data are unaffected.)

- [ ] **Step 8: Run the characterization tests and lint**

Run:
```bash
uv run pytest tests/test_cs.py tests/test_utils.py -v
uv run ruff check
```
Expected: all tests PASS; ruff reports no errors.

- [ ] **Step 9: Commit**

```bash
git add src/lqc/cs.py src/lqc/utils.py tests/test_cs.py
git commit -m "perf: drop intron regex, add grouped indel/mismatch getter, translate complement"
```

---

## Task 2: Refactor `stat.py` into record/prefetch/reduce API

**Files:**
- Modify: `src/lqc/stat.py` (imports + full body)
- Test: `tests/test_stat.py`

This introduces the parallel primitives but changes **no** public behavior: `stat_element_from_bam_by_contig` keeps its signature and results.

- [ ] **Step 1: Write the failing equivalence test**

Append to `tests/test_stat.py`:

```python
from lqc.stat import (
    prefetch_records,
    reduce_blocks_to_contigs,
    stat_element_from_bam_by_contig,
    stat_records,
)


def _assert_stat_tuples_equal(a, b):
    rs1, ins1, dele1, mis1, spl1 = a
    rs2, ins2, dele2, mis2, spl2 = b
    assert rs1.get_read_count() == rs2.get_read_count()
    assert rs1.get_total_base() == rs2.get_total_base()
    assert rs1.get_lengths() == rs2.get_lengths()
    assert rs1.get_insertions() == rs2.get_insertions()
    assert rs1.get_deletions() == rs2.get_deletions()
    assert rs1.get_mismatches() == rs2.get_mismatches()
    assert rs1.get_introns() == rs2.get_introns()
    assert ins1.get_indel_count() == ins2.get_indel_count()
    assert dele1.get_indel_count() == dele2.get_indel_count()
    assert mis1.get_type_count() == mis2.get_type_count()
    assert spl1.get_splice_pair_count_dict() == spl2.get_splice_pair_count_dict()


def test_prefetch_and_reduce_matches_serial(cs_bam, cs_bam_with_indel):
    for bam in (cs_bam, cs_bam_with_indel):
        serial = stat_element_from_bam_by_contig(
            bam_file=bam, genome_file=None, contig='chr1', method='cs'
        )
        contig, records = prefetch_records(bam, ['chr1'], 'cs')[0]
        half = (len(records) + 1) // 2
        chunks = [records[:half], records[half:]]
        block_results = [
            stat_records((i, contig, chunk), None, 'cs')
            for i, chunk in enumerate(chunks)
            if chunk
        ]
        reduced = reduce_blocks_to_contigs(
            [
                (c, rs, ins, dele, mis, spl)
                for c, rs, ins, dele, mis, spl, _ in block_results
            ],
            ['chr1'],
        )
        _assert_stat_tuples_equal(reduced[0], serial)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_stat.py::test_prefetch_and_reduce_matches_serial -v`
Expected: FAIL with `ImportError: cannot import name 'prefetch_records' from 'lqc.stat'`

- [ ] **Step 3: Add imports to `stat.py`**

Replace the top of `src/lqc/stat.py` (lines 1-8):

```python
import pysam

from lqc.cs import CS
from lqc.indel import Indel
from lqc.mismatch import Mismatch
from lqc.readstat import ReadStat
from lqc.splice import Splice
from lqc.utils import bam_or_sam, convert_complement, convert_reverse_complement
```

with:

```python
import os
from collections import defaultdict
from typing import NamedTuple, Optional

import pysam

from lqc.cs import CS
from lqc.indel import Indel
from lqc.mismatch import Mismatch
from lqc.readstat import ReadStat
from lqc.splice import Splice
from lqc.utils import bam_or_sam, convert_complement, convert_reverse_complement


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
    query_name: str
    cs_string: Optional[str] = None
    cigar: Optional[str] = None
    md_string: Optional[str] = None
    query_sequence: Optional[str] = None
    reference_end: Optional[int] = None


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
            query_name=read.query_name,
            cs_string=read.get_tag('cs'),
        )
    return ReadRecord(
        method='md',
        contig=contig,
        start_pos=read.reference_start,
        strand=strand,
        query_length=len(read.query_sequence),
        query_name=read.query_name,
        cigar=read.cigarstring,
        md_string=read.get_tag('MD'),
        query_sequence=read.query_sequence,
        reference_end=read.reference_end,
    )


def _build_cs(record, genome):
    if record.method == 'cs':
        return CS.from_cs_tag_string(
            record.cs_string,
            contig=record.contig,
            start_pos=record.start_pos,
            strand=record.strand,
        )
    ref_seq = genome.fetch(
        record.contig, record.start_pos, record.reference_end
    )
    return CS.from_cigar_string(
        record.cigar,
        record.md_string,
        record.query_sequence,
        ref_seq,
        contig=record.contig,
        start_pos=record.start_pos,
        strand=record.strand,
    )


def process_record(record, genome,
                   readstat, insertion, deletion, mismatch, splice):
    """Accumulate one read into the five stat objects and return its CS."""
    cs = _build_cs(record, genome)
    strand = record.strand
    readstat.add_read(
        length=record.query_length,
        insertion=cs.get_insertion_count(),
        deletion=cs.get_deletion_count(),
        mismatch=cs.get_mismatch_count(),
        intron=cs.get_intron_count(),
    )
    insertions, deletions, mismatches = cs.get_indel_mismatches(
        coordinate='normalized_read'
    )
    for a, _, _, d in insertions:
        insertion.add_indel(
            d if strand == '+' else convert_reverse_complement(d), a
        )
    for a, _, _, d in deletions:
        deletion.add_indel(
            d if strand == '+' else convert_reverse_complement(d), a
        )
    for a, _, _, d in mismatches:
        mismatch.add_mismatch(
            d if strand == '+' else convert_complement(d), a
        )
    for s_str, s_ct in cs.get_splice_pair_count().items():
        splice.add_splice_pair_count_dict({
            s_str if strand == '+' else convert_reverse_complement(s_str): s_ct
        })
    return cs


def prefetch_records(bam_file, contigs, method):
    """Return ``[(contig, [ReadRecord, ...]), ...]`` in BAM traversal order."""
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    per_contig = []
    for contig in contigs:
        records = [
            record_from_read(read, contig, method)
            for read in bam.fetch(contig)
        ]
        per_contig.append((contig, records))
    bam.close()
    return per_contig


def stat_records(task, genome_file, method, cs_dir=None):
    """Process one ``(index, contig, records)`` task; return per-chunk stats.

    When ``cs_dir`` is set, also write this chunk's read.cs lines to a
    per-chunk temp file and return its path for ordered merging.
    """
    index, contig, records = task
    readstat = ReadStat(contig)
    insertion = Indel(contig)
    deletion = Indel(contig)
    mismatch = Mismatch(contig)
    splice = Splice(contig)

    genome = None
    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)

    cs_path = None
    cs_fh = None
    if cs_dir is not None:
        cs_path = os.path.join(cs_dir, f'.readcs-{index:08d}.tmp')
        cs_fh = open(cs_path, 'w')

    for record in records:
        cs = process_record(
            record, genome,
            readstat, insertion, deletion, mismatch, splice
        )
        if cs_fh is not None:
            for low, high, mark, value in cs.get_contig_position():
                cs_fh.write(
                    f'{record.query_name}\t{record.contig}\t'
                    f'{low}\t{high}\t{mark}\t{value}\n'
                )

    if cs_fh is not None:
        cs_fh.close()
    if genome is not None:
        genome.close()

    return contig, readstat, insertion, deletion, mismatch, splice, cs_path


def reduce_blocks_to_contigs(block_results, contigs):
    """Fold per-chunk ``(contig, rs, ins, dele, mis, spl)`` into per-contig tuples."""
    by_contig = defaultdict(list)
    for contig, rs, ins, dele, mis, spl in block_results:
        by_contig[contig].append((rs, ins, dele, mis, spl))
    result = []
    for contig in contigs:
        blocks = by_contig.get(contig)
        if not blocks:
            continue
        rss, inss, deles, miss, spls = zip(*blocks)
        result.append((
            sum(rss), sum(inss), sum(deles), sum(miss), sum(spls)
        ))
    return result
```

- [ ] **Step 4: Re-implement `stat_element_from_bam_by_contig`**

Replace the whole existing `stat_element_from_bam_by_contig` body (lines 11-119) with:

```python
def stat_element_from_bam_by_contig(bam_file,
                                    genome_file,
                                    contig,
                                    method='cs'):
    assert method in ['cs', 'MD', 'both'],\
        "method should be either: cs, MD, both."
    readstat = ReadStat(contig)
    insertion = Indel(contig)
    deletion = Indel(contig)
    mismatch = Mismatch(contig)
    splice = Splice(contig)

    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)

    genome = None
    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)

    for read in bam.fetch(contig):
        record = record_from_read(read, contig, method)
        process_record(
            record, genome,
            readstat, insertion, deletion, mismatch, splice
        )

    bam.close()
    if genome is not None:
        genome.close()

    return readstat, insertion, deletion, mismatch, splice
```

Keep the trailing `########################################` (line 122) if removing the function left two separators; ensure exactly one blank line separates module items.

- [ ] **Step 5: Run the full stat test module**

Run: `uv run pytest tests/test_stat.py -v`
Expected: all tests PASS, including the new `test_prefetch_and_reduce_matches_serial`.

- [ ] **Step 6: Run regression + lint**

Run:
```bash
uv run pytest -q
uv run ruff check
```
Expected: full suite green; ruff reports no errors.

- [ ] **Step 7: Commit**

```bash
git add src/lqc/stat.py tests/test_stat.py
git commit -m "refactor: extract record/prefetch/reduce primitives from stat"
```

---

## Task 3: Wire data-parallel execution into the CLI

**Files:**
- Modify: `src/lqc/cli.py`
- Test: `tests/test_cli.py`

- [ ] **Step 1: Write the failing determinism test**

Append to `tests/test_cli.py`:

```python
def test_tables_identical_across_thread_counts(cs_bam, tmp_path):
    out1 = tmp_path / 'o1'
    out2 = tmp_path / 'o2'
    assert main(['-b', cs_bam, '-o', str(out1), '-c', 'chr1', '-t', '1']) == 0
    assert main(['-b', cs_bam, '-o', str(out2), '-c', 'chr1', '-t', '2']) == 0
    for name in (
        'read_stat.txt', 'insertion.txt', 'deletion.txt',
        'mismatch.txt', 'splice.txt',
    ):
        assert (out1 / 'table' / name).read_bytes() == \
               (out2 / 'table' / name).read_bytes()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_cli.py::test_tables_identical_across_thread_counts -v`
Expected: PASS today (`-t 2` and `-t 1` are both serial over one contig, so tables match). This is a characterization test that must remain green through the parallel rewrite.

- [ ] **Step 3: Add imports and helpers to `cli.py`**

In `src/lqc/cli.py`, add to the existing `from lqc import (...)` block (lines 19-50) — remove `stat_element_from_bam_by_contig` and `write_readcs` is kept for now — and add a new import line after it:

```python
from lqc.stat import prefetch_records, reduce_blocks_to_contigs, stat_records
```

Add `_chunk` near `get_stat_list` (after line 82):

```python
def _chunk(records, n):
    n = max(1, int(n))
    if not records:
        return []
    size = (len(records) + n - 1) // n
    return [records[i:i + size] for i in range(0, len(records), size)]
```

- [ ] **Step 4: Delete `stat_bam`**

Delete the `stat_bam` function (lines 57-63).

- [ ] **Step 5: Replace the per-contig pool with prefetch/chunk/map/reduce**

Replace lines 375-392 (the `# run jobs by contigs` block) with:

```python
    # run jobs by chunked reads (data-parallel, independent of contig count)
    message = 'Element statistic process starts.'
    logger.info(message)
    per_contig = prefetch_records(
        args['bam_file'], contigs, bam_type
    )
    tasks = []
    for contig, records in per_contig:
        if not records:
            continue
        for chunk in _chunk(records, args['thread']):
            tasks.append((len(tasks), contig, chunk))

    with mp.Pool(args['thread']) as p:
        block_results = p.map(
            partial(
                stat_records,
                genome_file = args['genome_fasta'],
                method = bam_type,
                cs_dir = None
            ),
            tasks
        )

    result = reduce_blocks_to_contigs(
        [
            (contig, rs, ins, dele, mis, spl)
            for contig, rs, ins, dele, mis, spl, _ in block_results
        ],
        contigs
    )
    message = 'Element statistic process finished.'
    logger.info(message)
```

- [ ] **Step 6: Change `-t` default**

Replace lines 169-174 with:

```python
    parser.add_argument(
        '-t', '--thread',
        help = 'threads to be used in calculation',
        type = int,
        default = min(os.cpu_count() or 1, 4)
    )
```

- [ ] **Step 7: Run tests and lint**

Run:
```bash
uv run pytest -q
uv run ruff check
```
Expected: full suite green (including `test_tables_identical_across_thread_counts` and the existing end-to-end `--output-cs` test, which still goes through `write_readcs`); ruff clean.

- [ ] **Step 8: Commit**

```bash
git add src/lqc/cli.py tests/test_cli.py
git commit -m "feat: parallelize stat over read chunks; default threads to min(cpu,4)"
```

---

## Task 4: Single-pass `--output-cs`

**Files:**
- Modify: `src/lqc/cli.py`
- Test: `tests/test_cli.py`

- [ ] **Step 1: Write the failing read.cs determinism test**

Append to `tests/test_cli.py`:

```python
def test_output_cs_identical_across_thread_counts(cs_bam, tmp_path):
    out1 = tmp_path / 'o1'
    out2 = tmp_path / 'o2'
    main(['-b', cs_bam, '-o', str(out1), '-c', 'chr1', '-t', '1', '--output-cs'])
    main(['-b', cs_bam, '-o', str(out2), '-c', 'chr1', '-t', '2', '--output-cs'])
    assert (out1 / 'read.cs').read_bytes() == (out2 / 'read.cs').read_bytes()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_cli.py::test_output_cs_identical_across_thread_counts -v`
Expected: FAIL or is meaningless while `--output-cs` still runs the old `write_readcs` second pass and `stat_records` ignores `cs_dir`. (It may actually pass content-wise; the point is to introduce the assertion that gates the new path.)

- [ ] **Step 3: Route `cs_dir` into the workers**

In the pool block added in Task 3, change the `partial` call:

```python
            partial(
                stat_records,
                genome_file = args['genome_fasta'],
                method = bam_type,
                cs_dir = None
            ),
```

to:

```python
            partial(
                stat_records,
                genome_file = args['genome_fasta'],
                method = bam_type,
                cs_dir = o_dirs['base'] if args['output_cs'] else None
            ),
```

- [ ] **Step 4: Replace the `write_readcs` block with ordered merge**

Replace lines 462-470 (the current `if args['output_cs']:` block) with:

```python
    if args['output_cs']:
        message = 'Output processed cs tags.'
        logger.info(message)
        with open(o_files['cs'], 'w') as out:
            out.write(
                'read_name\tcontig\tlow\thigh\tcs_mark\tcs_value\n'
            )
            for *_, cs_path in block_results:
                with open(cs_path, 'r') as fh:
                    out.write(fh.read())
                os.remove(cs_path)
        message = 'Output processed cs tags finished.'
        logger.debug(message)
    else:
        pass
```

- [ ] **Step 5: Drop the now-unused `write_readcs` import**

In `src/lqc/cli.py`, remove `write_readcs,` from the `from lqc import (...)` block (lines 19-50). Leave `write_readcs` in `src/lqc/utils.py` for library callers.

- [ ] **Step 6: Run tests and lint**

Run:
```bash
uv run pytest -q
uv run ruff check
```
Expected: full suite green; `test_output_cs_identical_across_thread_counts` green; ruff clean.

- [ ] **Step 7: Manual determinism check on real data (optional, ~8 min)**

Run:
```bash
bash tmp/analyze_chr22.sh --output-cs
bash tmp/analyze_chr22.sh --output-cs -t 2
diff <(cat tmp/out-src/table/*.txt) <(cat tmp/out-src/table/*.txt)   # tables unchanged
```
Expected: run completes; `read.cs` and tables byte-identical across `-t 4` and `-t 2`. (Do this if the ~200 MB BAM is present.)

- [ ] **Step 8: Commit**

```bash
git add src/lqc/cli.py src/lqc/utils.py tests/test_cli.py
git commit -m "perf: single-pass --output-cs emitted from stat workers"
```

---

## Self-Review

**Spec coverage** — Stage 1 micro-opts (get_tag is folded into Task 2's `record_from_read`; regex removal, grouped getter, translate, dead-code removal in Task 1); Stage 2 parallelization (Task 2 + Task 3); Stage 3 single-pass `--output-cs` (Task 4); `-t` default change (Task 3). Success criteria (bit-identical tables/read.cs, determinism, ruff+pytest) are encoded as the Task 3/4 tests and the Task 7 verification steps.

**Placeholder scan** — no TBD/TODO/"similar to"/"add tests for the above" left; every code step shows full code and exact commands.

**Type consistency** — `ReadRecord` field names (`cs_string`, `query_sequence`, `reference_end`, etc.) are used identically in `record_from_read`, `_build_cs`, and `process_record`. `stat_records` returns a 7-tuple `(contig, readstat, insertion, deletion, mismatch, splice, cs_path)`; both `reduce_blocks_to_contigs` call sites and the read.cs merge unpack it in that exact order. `get_indel_mismatches` returns `(insertions, deletions, mismatches)` and is consumed in that order in `process_record`.