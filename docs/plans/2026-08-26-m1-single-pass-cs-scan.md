# M1 — Single-pass `cs` scan + fused count/location Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fuse the per-read `cs`-tag processing from ~6 passes (regex tokenize → `_compute_counts` → read-location → strand mirror → normalize → filter) down to 2 cheap passes, while keeping every output byte-identical.

**Architecture:** Optimize only the `cs`-tag path in `src/lqc/cs.py` (design spec §4.2, root cause P1). (1) Replace the regex-based `cs_to_list` tokenizer with a manual single-pass scan. (2) In `CS.__init__`, fold `_compute_counts` and the `+`/`-`/`*` normalized-read-location computation into one pass over the element list, caching the result in `self._indel_mismatch_normalized`. (3) `get_indel_mismatches(coordinate='normalized_read')` serves straight from that cache. All other `CS` public methods keep their existing lazy behavior; the `MD`/CIGAR path is untouched end-to-end (it only shares the now-fused `__init__`, with byte-identical results).

**Tech Stack:** Python 3.11+ (stdlib only; `numpy` not used in this milestone), pytest, ruff. Run test/lint via `.venv/bin/python -m pytest` / `.venv/bin/python -m ruff check` (the `uv` cache is read-only in this sandbox — never `uv run`).

---

## File Structure

- **Modify `src/lqc/cs.py`** — the only production file touched.
  - `cs_to_list` (lines 8-35): manual tokenizer (Task 1).
  - `CS.__init__` (lines 337-355) + `_compute_counts` (357-405): fused pass, renamed `_compute_counts_and_indel_mismatch_locations` (Task 2).
  - `get_indel_mismatches` (594-608): normalize fast path (Task 2).
- **Modify `tests/test_cs.py`** — add tokenizer-parity and fused-location tests.
- **No other file changes.** `src/lqc/stat.py`, `indel.py`, `mismatch.py`, `readstat.py` are consumers and are unchanged (their `add_*` methods only read the `[low, high, mark, value]` items — verified in `process_record`: `for a, _, _, d in insertions:`).

### Byte-parity invariants (must hold exactly)

1. `cs_to_list(s)` returns the same `list[[ref_low, ref_high, mark, value], ...]` as the old regex implementation for every `s`.
2. `get_indel_mismatches('normalized_read')` returns the same three lists as the old lazy path:
   - strand `+`: `[read_low/read_len, read_high/read_len, mark, value]`.
   - strand `-`: `[(read_len-read_high)/read_len, (read_len-read_low)/read_len, mark, value]`.
   These use Python 3 `int/int` (float) division, so the float64 values are bit-identical.
3. The 11 precomputed count attributes are unchanged in value.

---

## Task 1: Manual single-pass `cs_to_list`

**Files:**
- Modify: `src/lqc/cs.py:8-35`
- Test: `tests/test_cs.py`

The value character set `[0-9a-z]` is disjoint from the mark set `{:,*,+,-,~}`, so the tokenizer can scan one mark then one maximal digit/lowercase run with no regex.

- [ ] **Step 1: Write the failing parity tests**

Add to `tests/test_cs.py` (below the existing `test_cs_to_list`), including a reference copy of the OLD regex implementation used only for comparison:

```python
import re

# NOTE: CS_STRING is already defined at module scope (line 54 of test_cs.py);
# do NOT redefine it here.


def _ref_cs_to_list(cs_string):
    """Reference copy of the pre-M1 regex tokenizer, for byte-parity only."""
    pos = 0
    cslenfuncs = {
        ':': int,
        '*': lambda x: 1,
        '+': lambda x: 0,
        '-': len,
        '~': lambda x: int(re.sub('[a-z]', '', x)),
    }
    cs_mark = re.sub('[0-9a-z]', ' ', cs_string).strip().split()
    cs_value = re.sub(r'[:*\-+~]', ' ', cs_string).strip().split()
    cslist = []
    for a, b in zip(cs_mark, cs_value, strict=True):
        low = pos
        pos += cslenfuncs[a](b)
        high = pos
        cslist.append([low, high, a, b])
    return cslist


def test_cs_to_list_all_marks():
    assert cs_to_list(CS_STRING) == [
        [0, 10, ':', '10'],
        [10, 11, '*', 'ag'],
        [11, 16, ':', '5'],
        [16, 16, '+', 'tt'],
        [16, 19, ':', '3'],
        [19, 119, '~', 'gt100ag'],
        [119, 123, ':', '4'],
        [123, 126, '-', 'ccc'],
        [126, 128, ':', '2'],
    ]


PARITY_STRINGS = [
    ':29*ga:19',
    CS_STRING,
    '*ag',
    ':10',
    '+atcg',
    '-c',
    '~gt100ag',
    ':4~ct636ac:5',
    ':29*ga:19*at:61~ct140ac:45*gc',
]


@pytest.mark.parametrize('cs_string', PARITY_STRINGS)
def test_cs_to_list_matches_reference(cs_string):
    assert cs_to_list(cs_string) == _ref_cs_to_list(cs_string)


def test_cs_to_list_matches_reference_on_fixtures(cs_test_data_records):
    for record in cs_test_data_records:
        cs_string = record[2]
        assert cs_to_list(cs_string) == _ref_cs_to_list(cs_string), record[0]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_cs.py -q`
Expected: `test_cs_to_list_all_marks` and `test_cs_to_list_matches_reference*` PASS (the old regex `cs_to_list` already produces identical output); this establishes the baseline before refactoring. NOTE: these are parity tests, so they are expected to PASS against the old code too — the point is they keep passing after the rewrite.

- [ ] **Step 3: Replace `cs_to_list` with the manual tokenizer**

Replace the body of `cs_to_list` in `src/lqc/cs.py` (keep the docstring):

```python
def cs_to_list(cs_string):
    '''convert the cs_strings into a list with values and
       genomic positions (0-based) [low, high, cs tag, cs value].'''
    # Manual single-pass scan. Marks {:,*,+,-,~} are disjoint from the value
    # character set ([0-9a-z]), so read one mark then one maximal value run.
    cslist = []
    pos = 0
    i = 0
    n = len(cs_string)
    while i < n:
        mark = cs_string[i]
        i += 1
        j = i
        while j < n and (
            ('0' <= cs_string[j] <= '9') or ('a' <= cs_string[j] <= 'z')
        ):
            j += 1
        value = cs_string[i:j]
        i = j
        low = pos
        if mark == ':':
            pos += int(value)
        elif mark == '*':
            pos += 1
        elif mark == '-':
            pos += len(value)
        elif mark == '~':
            # intron reference length = the digit run (drop donor/acceptor letters)
            pos += int(''.join(c for c in value if '0' <= c <= '9'))
        elif mark == '+':
            pass
        else:
            pass
        cslist.append([low, pos, mark, value])
    return cslist
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_cs.py -q`
Expected: all pass. Also run ruff: `.venv/bin/python -m ruff check src/lqc/cs.py tests/test_cs.py` → no errors.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/cs.py tests/test_cs.py
git commit -m "perf: manual single-pass cs_to_list tokenizer (drop regex)"
```

---

## Task 2: Fused counts + normalized indel/mismatch locations in `CS.__init__`

**Files:**
- Modify: `src/lqc/cs.py:337-405` (`__init__` + `_compute_counts` → fused)
- Modify: `src/lqc/cs.py:594-608` (`get_indel_mismatches` fast path)
- Test: `tests/test_cs.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_cs.py`:

```python
def test_indel_mismatch_normalized_eager_matches_lazy(cs_obj):
    # cached attribute is in cs-tag order (for CS_STRING: *, +, -)
    cached = cs_obj._indel_mismatch_normalized
    assert [a[2] for a in cached] == ['*', '+', '-']
    ins, dels, mism = cs_obj.get_indel_mismatches(coordinate='normalized_read')
    # the grouped getter must still match the single-mark lazy getters
    assert ins == cs_obj.get_insertions(coordinate='normalized_read')
    assert dels == cs_obj.get_deletions(coordinate='normalized_read')
    assert mism == cs_obj.get_mismatches(coordinate='normalized_read')
    # regroup the cached items by mark; must reproduce the grouped getter exactly
    by_mark = {'+': [], '-': [], '*': []}
    for a in cached:
        by_mark[a[2]].append(a)
    assert by_mark['+'] == ins
    assert by_mark['-'] == dels
    assert by_mark['*'] == mism


def test_get_indel_mismatches_reverse_strand_mirrors():
    rev = CS.from_cs_tag_string(
        CS_STRING, contig='chr1', start_pos=1000, strand='-'
    )
    ins, dels, mism = rev.get_indel_mismatches(coordinate='normalized_read')
    # read_len == 27; *ag at read_low=10, read_high=11 -> mirror to (16,17)
    assert len(mism) == 1
    assert mism[0][0] == (27 - 11) / 27
    assert mism[0][1] == (27 - 10) / 27
    assert mism[0][2] == '*'
    assert mism[0][3] == 'ag'
    # +tt at read_low=16, read_high=18 -> mirror to (9,11)
    assert len(ins) == 1
    assert ins[0][0] == (27 - 18) / 27
    assert ins[0][1] == (27 - 16) / 27
    # -ccc at read_low=25, read_high=25 -> mirror to (2,2)
    assert len(dels) == 1
    assert dels[0][0] == (27 - 25) / 27
    assert dels[0][1] == (27 - 25) / 27


def test_indel_mismatch_other_coordinates_unchanged(cs_obj):
    # non-normalized coordinates still go through the lazy path
    ins_c, dels_c, mism_c = cs_obj.get_indel_mismatches(coordinate='contig')
    assert [a[2] for a in ins_c] == ['+']
    assert dels_c[0][2] == '-'
    assert mism_c[0][2] == '*'
    assert mism_c[0][0] == 1010  # start_pos 1000 + ref_low 10
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_cs.py -q`
Expected: `test_indel_mismatch_normalized_eager_matches_lazy` FAILS with `AttributeError: 'CS' object has no attribute '_indel_mismatch_normalized'`. The other two should pass against current code (they only exercise existing behavior).

- [ ] **Step 3: Implement the fused `__init__` pass**

Replace `__init__` (lines 337-355) and `_compute_counts` (357-405) in `src/lqc/cs.py` with:

```python
    def __init__(self, cs_tag_list, contig, start_pos, strand):
        super().__init__()
        self._cs = cs_tag_list
        self._start_pos = start_pos
        self._contig = contig
        self._strand = strand
        self._read_location = None
        self._normalized_read_location = None
        (self._read_length,
         self._intron_count,
         self._splice_pair_count,
         self._splice_l_count,
         self._splice_r_count,
         self._mismatch_count,
         self._mismatch_type_count,
         self._insertion_count,
         self._insertion_length,
         self._deletion_count,
         self._deletion_length,
         self._indel_mismatch_normalized) = \
            self._compute_counts_and_indel_mismatch_locations()

    def _compute_counts_and_indel_mismatch_locations(self):
        """One pass over the cs list: precomputed metrics AND the normalized-read
        coordinates of ``+``/``-``/``*`` elements (design spec §4.2)."""
        read_pos = 0
        intron_count = 0
        splice_pair_count = Counter()
        splice_l_count = Counter()
        splice_r_count = Counter()
        mismatch_count = 0
        mismatch_type_count = Counter()
        insertion_count = 0
        insertion_length = 0
        deletion_count = 0
        deletion_length = 0
        raw_coords = []   # (read_low, read_high, mark, value) for +,-,*
        for _, _, mark, value in self._cs:
            read_low = read_pos
            if mark == ':':
                read_pos += int(value)
            elif mark == '~':
                intron_count += 1
                donor = value[:2]
                acceptor = value[-2:]
                splice_pair_count[f'{donor}-{acceptor}'] += 1
                splice_l_count[donor] += 1
                splice_r_count[acceptor] += 1
            elif mark == '+':
                read_pos += len(value)
                insertion_count += 1
                insertion_length += len(value)
                raw_coords.append((read_low, read_pos, mark, value))
            elif mark == '*':
                read_pos += 1
                mismatch_count += 1
                mismatch_type_count[f'{value[0]}{value[1]}'] += 1
                raw_coords.append((read_low, read_pos, mark, value))
            elif mark == '-':
                deletion_count += 1
                deletion_length += len(value)
                raw_coords.append((read_low, read_pos, mark, value))
            else:
                pass
        read_len = read_pos
        strand = self._strand
        normalized = []
        for read_low, read_high, mark, value in raw_coords:
            if strand == '+':
                low_n = read_low / read_len
                high_n = read_high / read_len
            else:
                low_n = (read_len - read_high) / read_len
                high_n = (read_len - read_low) / read_len
            normalized.append([low_n, high_n, mark, value])
        return (read_len,
                intron_count,
                splice_pair_count,
                splice_l_count,
                splice_r_count,
                mismatch_count,
                mismatch_type_count,
                insertion_count,
                insertion_length,
                deletion_count,
                deletion_length,
                normalized)
```

- [ ] **Step 4: Add the `normalized_read` fast path to `get_indel_mismatches`**

Replace `get_indel_mismatches` (lines 594-608) with:

```python
    def get_indel_mismatches(self, coordinate='normalized_read'):
        """Return (insertions, deletions, mismatches) in one pass, avoiding
        three separate full-list scans. ``normalized_read`` serves from the
        ``__init__``-cached attribute (no recompute)."""
        insertions = []
        deletions = []
        mismatches = []
        if coordinate == 'normalized_read':
            source = self._indel_mismatch_normalized
        else:
            source = self._get_coordinate_list(coordinate)
        for a in source:
            mark = a[2]
            if mark == '+':
                insertions.append(a)
            elif mark == '-':
                deletions.append(a)
            elif mark == '*':
                mismatches.append(a)
        return insertions, deletions, mismatches
```

- [ ] **Step 5: Run the full `cs` test module**

Run: `.venv/bin/python -m pytest tests/test_cs.py -q`
Expected: all pass (including the new Task 2 tests and the pre-existing `test_get_indel_mismatches_groups`, `test_element_count`, `test_cs_and_cigar_md_produce_identical_lists`).

- [ ] **Step 6: Run ruff**

Run: `.venv/bin/python -m ruff check src/lqc/cs.py tests/test_cs.py`
Expected: no errors (the reference helper `_ref_cs_to_list` in the test uses `re`, which is imported at test module top — ensure `import re` is present in `tests/test_cs.py`).

- [ ] **Step 7: Commit**

```bash
git add src/lqc/cs.py tests/test_cs.py
git commit -m "perf: fuse cs counts + normalized indel/mismatch locations, cache get_indel_mismatches"
```

---

## Task 3: Full-suite + byte-parity + performance verification

**Files:** none (verification only).

- [ ] **Step 1: Full test suite**

Run: `.venv/bin/python -m pytest -q`
Expected: `118 passed` (no new tests were added to the running total except the parametrized parity cases; if the count differs, confirm every test passes and note the new total).

- [ ] **Step 2: chr22 byte-parity vs baseline**

Run:
```bash
rm -rf tmp/m1/chr22 && mkdir -p tmp/m1/chr22
.venv/bin/lqc -b tmp/full_bam/ENCFF417VHJ.sorted.bam --subsample-contigs chr22 -o tmp/m1/chr22 -t 1 --output-cs
```
(If the CLI flag for a single contig is `--contig chr22` or similar, use the exact flag from `lqc --help`; the goal is a `-t 1` chr22 run producing `table/`, `fig/`, `LQC_report.html`, `read.cs`.)

Then compare the deterministic artifacts against `tmp/verify/chr22` (the accepted M0 baseline):
```bash
for f in table/*.txt read.cs LQC_report.html; do
  cmp tmp/m1/chr22/$f tmp/verify/chr22/$f && echo "OK $f" || echo "DIFF $f"
done
for f in fig/*.png; do
  cmp tmp/m1/chr22/$f tmp/verify/chr22/$f && echo "OK $f" || echo "DIFF $f"
done
```
Expected: every line `OK`. `fig/*.pdf` embeds a creation timestamp and is excluded.

- [ ] **Step 3: chr22 stat-phase timing (CPU win)**

Run: `.venv/bin/python tmp/profile_lqc.py` (or the equivalent single-contig `stat_element_from_bam_by_contig` profile) on chr22.
Expected: per-read wall time drops (the ~35% location recompute + ~27% regex tokenization are removed). Record the before/after numbers in a comment here.

- [ ] **Step 4: Record results in this plan (no code change)**

Note the byte-parity outcome and the new timing next to the M0 baseline (chr22 stat ~16.7 s single-threaded). If byte-parity fails, STOP and investigate before proceeding to M2.

---

## Self-Review (run before handoff)

- **Spec coverage:** §4.2 of the design spec — manual tokenizer (Task 1 covers §4.2 step 1 + the `~`-intron `re.sub` removal), same-loop count accumulation (Task 2 covers step 2), `+`/`-`/`*` raw read-coord recording + post-loop normalize/strand-mirror (Task 2 covers steps 3-4), `get_indel_mismatches('normalized_read')` served from cache (Task 2 Step 4 covers "Public-API preservation"), float semantics unchanged (`int/int`, Task 2 Step 3), `MD`/CIGAR scope left alone (no `convert_cigar_md_to_cs_list`/`cigar_to_list`/`md_to_list` touched). Covered.
- **Placeholder scan:** no TBD/TODO; every code step shows complete code; every test shows complete code; commands include expected output.
- **Type consistency:** `_compute_counts_and_indel_mismatch_locations` returns a 12-tuple `(read_len, intron_count, splice_pair_count, splice_l_count, splice_r_count, mismatch_count, mismatch_type_count, insertion_count, insertion_length, deletion_count, deletion_length, normalized)` — matches the 12-name unpack in `__init__` exactly; `_indel_mismatch_normalized` is a `list[[low_n, high_n, mark, value]]` assigned to `self` and read by `get_indel_mismatches`. Consistent.
- **Conflicts with Task 2 Step 1 imports:** `tests/test_cs.py` must import `pytest` (already does) and now `re` and `CS_STRING`. The existing file already defines `CS_STRING` at module scope (line 54); the plan re-declares it in the test snippet for self-containment — during implementation, do NOT duplicate the `CS_STRING` definition; reuse the existing one.