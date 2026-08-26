# M2 — numpy columnar storage, cheap merges, cached aggregates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the boxed per-read/per-event Python lists and their `deepcopy`-based `__add__`/`sum()` merges with numpy columnar storage and `np.concatenate`, so the full-genome run completes in a few GB of RAM (fixes the P2 `'Total'`-sum OOM at 8.4 GB) while keeping every output byte-identical.

**Architecture:** Optimize the four accumulator classes (`readstat.py`, `indel.py`, `mismatch.py`, `splice.py`) only. Accumulation still appends to plain Python lists (cheap, chunk-bounded after M0); a private `_finalize()` converts them to numpy arrays lazily and drops the boxed list. `__add__` concatenates finalized arrays with `np.concatenate` (never `deepcopy`). Every public getter keeps its signature and returns Python scalar/list (`.tolist()` at the boundary) so `report_table`/`report_figure` stay **zero-change**; histogram getters keep returning numpy `(edges, hist)` exactly as today. Labels stay strings; `'Total'` sentinel and `relabel` in `cli.main` are untouched.

**Tech Stack:** Python 3.11+, numpy (already a dependency), pytest, ruff. Run tests/lint via `.venv/bin/python -m pytest` / `.venv/bin/python -m ruff check` — never `uv run`.

### Precision & byte-parity invariants (must hold exactly)

1. Count/length columns are `int32`; normalized locations are `float64` (bit-identical to the Python `float` produced by `int / int` — both are IEEE-754 correctly-rounded doubles).
2. `float` divisions reproduce the old `a/b` Python semantics: `arr_col / arr_col` on int32 promotes to float64 and is bit-identical to Python `int/int` for these magnitudes.
3. `np.histogram` on an in-order float64 array produces the same `(edges, hist)` as `np.histogram` on the old in-order list (same values, same order).
4. Empty (zero-read/zero-event) objects must work: `_finalize()` on an empty boxed list yields shape `(0, N)`; `__add__` concatenating empties works; list getters return `[]`.
5. `__radd__` keeps returning `self` for `other == 0` (that is how `sum(gen)` starts and how single-block contigs in `reduce_blocks_to_contigs` stay correct).

### Drop-list hardening note
`_finalize()` sets the private boxed/list attribute to `None` after building its array. This is safe because, in the pipeline, `add_*` runs only during chunk accumulation (before reduce) and getters/`__add__` run only on reduced/final objects. The cumulative `sum()` builds every intermediate result via `__add__`, so intermediates hold numpy arrays — this is where the multi-GB boxed lists used to blow up, and they are now gone.

---

## File Structure

- **Modify `src/lqc/readstat.py`** (Task 1): `_reads` → lazy `int32 (n,7)` array; memoized reversed-sorted lengths for `get_length_NL`; vectorized aligned-fraction counters; `__add__` concatenates.
- **Modify `src/lqc/indel.py`** (Task 2): `_indels` → `(iidx int32, ilen int32, loc float64)` parallel arrays; direct `np.histogram`; `__add__` concatenates + re-maps string index.
- **Modify `src/lqc/mismatch.py`** (Task 3): `_type_locations` → per-type `float64` arrays; `__add__` concatenates per type.
- **Modify `src/lqc/splice.py`** (Task 4): drop `deepcopy` from `__add__`.
- **Modify tests**: `tests/test_readstat.py`, `tests/test_indel.py`, `tests/test_mismatch.py`, `tests/test_splice.py` — add empty-object merge/ordering tests.

---

## Task 1: `ReadStat` numpy columnar storage

**Files:**
- Modify: `src/lqc/readstat.py`
- Test: `tests/test_readstat.py`

- [ ] **Step 1: Add the empty/merge tests**

Append to `tests/test_readstat.py`:

```python
def test_add_empty_and_nonempty_via_sum():
    empty = ReadStat('e')
    full = ReadStat('f')
    full.add_read(100, 1, 2, 3, 4, mapping_quality=60, aligned_length=90)
    full.add_read(200, 5, 6, 7, 8, mapping_quality=30, aligned_length=200)
    total = sum([empty, full])
    assert total.get_read_count() == 2
    assert total.get_lengths() == [100, 200]
    assert total.get_total_base() == 300
```

- [ ] **Step 2: Run tests to verify they pass against old code**

Run: `.venv/bin/python -m pytest tests/test_readstat.py -q`
Expected: all pass (this is a characterization test — the current boxed-list implementation already satisfies it).

- [ ] **Step 3: Implement the numpy version**

In `src/lqc/readstat.py`: replace `from copy import deepcopy` with `import numpy as np`, and replace `__init__`, `add_read`, `get_lengths`, `get_insertions`, `get_deletions`, `get_mismatches`, `get_introns`, `get_length_normalized_insertions`, `get_length_normalized_deletions`, `get_length_normalized_mismatches`, `get_length_normalized_introns`, `get_mapping_qualities`, `get_aligned_lengths`, `get_aligned_fractions`, `get_length_NL`, `get_read_count_with_aligned_fraction_below`, `get_read_count_fully_aligned`, and `__add__` with:

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
        self._reads_array = None
        self._sorted_lengths = None

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
        self._reads_array = None
        self._sorted_lengths = None

    def _finalize(self):
        """Convert the boxed row list to an ``(n, 7)`` int32 array (lazily) and
        drop the boxed list. Columns: length, insertion, deletion, mismatch,
        intron, mapping_quality, aligned_length."""
        if self._reads_array is None:
            self._reads_array = np.asarray(
                self._reads, dtype = np.int32
            ).reshape(-1, 7)
            self._reads = None
        return self._reads_array
```

(Continue the replacements; place methods where the originals were, keep all other methods byte-identical.)

```python
    def get_read_count(self):
        return self._read_count

    def get_lengths(self):
        return self._finalize()[:, 0].tolist()

    def get_insertions(self):
        return self._finalize()[:, 1].tolist()

    def get_deletions(self):
        return self._finalize()[:, 2].tolist()

    def get_mismatches(self):
        return self._finalize()[:, 3].tolist()

    def get_introns(self):
        return self._finalize()[:, 4].tolist()

    def get_length_normalized_insertions(self):
        arr = self._finalize()
        return (arr[:, 1] / arr[:, 0]).tolist()

    def get_length_normalized_deletions(self):
        arr = self._finalize()
        return (arr[:, 2] / arr[:, 0]).tolist()

    def get_length_normalized_mismatches(self):
        arr = self._finalize()
        return (arr[:, 3] / arr[:, 0]).tolist()

    def get_length_normalized_introns(self):
        arr = self._finalize()
        return (arr[:, 4] / arr[:, 0]).tolist()

    def get_mapping_qualities(self):
        return self._finalize()[:, 5].tolist()

    def get_aligned_lengths(self):
        return self._finalize()[:, 6].tolist()

    def get_aligned_fractions(self):
        arr = self._finalize()
        length = arr[:, 0]
        aligned = arr[:, 6]
        fraction = np.where(length > 0, aligned / length, 0.0)
        return fraction.tolist()

    def get_read_count_with_aligned_fraction_below(self, threshold = 0.9):
        arr = self._finalize()
        length = arr[:, 0]
        aligned = arr[:, 6]
        fraction = np.where(length > 0, aligned / length, 0.0)
        return int(np.sum(fraction < threshold))

    def get_read_count_fully_aligned(self):
        arr = self._finalize()
        return int(np.sum(arr[:, 6] == arr[:, 0]))

    def _sorted_lengths_desc(self):
        if self._sorted_lengths is None:
            self._sorted_lengths = sorted(
                self.get_lengths(), reverse = True
            )
        return self._sorted_lengths

    def get_length_NL(self, percent):
        assert percent >= 0 and percent <= 100,\
            "percent value should be between 0 and 100."
        lengths = self._sorted_lengths_desc()
        basesum = 0
        previous_basesum = 0
        percentbase = self._total_base * percent / 100
        for i in range(self._read_count):
            length = lengths[i]
            previous_basesum = basesum
            basesum += length
            if (previous_basesum <= percentbase) \
               and (basesum >= percentbase):
                N = length
                L = i + 1
                break
            else:
                continue
        return N, L

    def __add__(self, other):
        assert isinstance(other, type(self)),\
            'wrong object to add'
        sumReadStat = type(self)(
            f'{self.label} {other.label}'
        )
        sumReadStat._reads_array = np.concatenate([
            self._finalize(), other._finalize()
        ])
        sumReadStat._reads = None
        sumReadStat._read_count = self._read_count + other._read_count
        sumReadStat._total_base = self._total_base + other._total_base
        sumReadStat._total_aligned_base = self._total_aligned_base + other._total_aligned_base
        return sumReadStat
```

IMPORTANT implementation notes:
- Do NOT change any other method (`get_total_base`, `get_mean_length`, `get_min/max_*`, `get_mean_*`, `*_per_base`, `*_per_aligned_base`, `_get_median`, `get_median_*`, `get_read_count_with_n_*`, `__repr__`, `__str__`, `__radd__`, `__len__`). They already call the list-returning getters above (which now `tolist()`).
- `__add__` produces an array-only result (`_reads = None`), so it must NOT be called on an object that will be `add_read`-into afterwards (the pipeline never does this; merged objects are final).

- [ ] **Step 4: Run tests + ruff**

Run: `.venv/bin/python -m pytest tests/test_readstat.py -q` and `.venv/bin/python -m ruff check src/lqc/readstat.py tests/test_readstat.py`
Expected: all pass; ruff clean (no unused `deepcopy` import — it is removed).

- [ ] **Step 5: Commit**

```bash
git add src/lqc/readstat.py tests/test_readstat.py
git commit -m "perf: numpy columnar ReadStat with concatenate merge + memoized N/L"
```

---

## Task 2: `Indel` numpy parallel arrays

**Files:**
- Modify: `src/lqc/indel.py`
- Test: `tests/test_indel.py`

- [ ] **Step 1: Add empty/merge tests**

Append to `tests/test_indel.py`:

```python
def test_indel_empty_plus_nonempty_merge():
    empty = Indel('e')
    full = Indel('f')
    full.add_indel('a', 0.1)
    full.add_indel('b', 0.5)
    total = sum([empty, full])
    assert total.get_total_count() == 2
    assert total.get_lengths() == [1, 1]
    assert total.get_locations() == [0.1, 0.5]
    assert total.get_indel_count() == {'a': 1, 'b': 1}


def test_indel_locations_order_preserved():
    a = Indel('a')
    b = Indel('b')
    a.add_indel('x', 0.25)
    a.add_indel('y', 0.9)
    b.add_indel('z', 0.33)
    total = a + b
    assert total.get_locations() == [0.25, 0.9, 0.33]
```

- [ ] **Step 2: run to verify they pass against old code** — `.venv/bin/python -m pytest tests/test_indel.py -q` (characterization).

- [ ] **Step 3: Implement**

In `src/lqc/indel.py`, replace `__init__`, `add_indel`, `get_indel_count`, `get_total_count`, `get_total_length`, `get_lengths`, `get_locations`, `get_location_histogram`, `convert_reverse_complement`, and `__add__` with:

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
        self._indel_strings = []
        self._indel_index = {}
        self._indels = []
        self._indel_arrays = None  # (iidx int32, ilen int32, loc float64)

    def add_indel(self, indel,
                  normalized_read_location):
        iidx = self._indel_index.get(indel)
        if iidx is None:
            iidx = len(self._indel_strings)
            self._indel_strings.append(indel)
            self._indel_index[indel] = iidx
        self._indels.append([
            iidx, len(indel),
            normalized_read_location
        ])
        self._indel_arrays = None

    def _finalize(self):
        if self._indel_arrays is None:
            rows = np.asarray(self._indels, dtype = np.float64).reshape(-1, 3)
            self._indel_arrays = (
                rows[:, 0].astype(np.int32, copy = False),
                rows[:, 1].astype(np.int32, copy = False),
                rows[:, 2],
            )
            self._indels = None
        return self._indel_arrays

    def get_indel_count(self):
        indel_count = Counter()
        iidx, _, _ = self._finalize()
        for i in iidx.tolist():
            indel_count[self._indel_strings[i]] += 1
        return indel_count

    def get_total_count(self):
        return int(self._finalize()[0].size)

    def get_total_length(self):
        return int(self._finalize()[1].sum())

    def get_lengths(self):
        return self._finalize()[1].tolist()

    def get_locations(self):
        return self._finalize()[2].tolist()

    def get_location_histogram(self, cuts = None):
        if cuts is None:
            cuts = [i / 10 for i in range(11)]
        hist, edges = np.histogram(
            self._finalize()[2], bins = cuts, density = False
        )
        return edges, hist

    def convert_reverse_complement(self):
        newIndel = type(self)(self.label)
        iidx, ilen, loc = self._finalize()
        newIndel._indel_arrays = (iidx.copy(), ilen.copy(), loc.copy())
        newIndel._indels = None
        new_strings = [
            convert_reverse_complement(indel)
            for indel in self._indel_strings
        ]
        newIndel._indel_strings = new_strings
        newIndel._indel_index = {
            indel: i for i, indel in enumerate(new_strings)
        }
        return newIndel

    def __add__(self, other):
        assert isinstance(other, type(self)), 'wrong object to add'
        newIndel = type(self)(
            f'{self.label} {other.label}'
        )
        siidx, silen, sloc = self._finalize()
        oiidx, oilen, oloc = other._finalize()
        # string table: keep self's ordering, append other's new strings once
        new_strings = list(self._indel_strings)
        string_index = {
            indel: i for i, indel in enumerate(new_strings)
        }
        other_new_idx = []
        for indel in other._indel_strings:
            if indel in string_index:
                other_new_idx.append(string_index[indel])
            else:
                string_index[indel] = len(new_strings)
                new_strings.append(indel)
                other_new_idx.append(len(new_strings) - 1)
        remapped_oiidx = np.asarray(
            [other_new_idx[i] for i in oiidx.tolist()],
            dtype = np.int32,
        )
        newIndel._indel_arrays = (
            np.concatenate([siidx, remapped_oiidx]),
            np.concatenate([silen, oilen]),
            np.concatenate([sloc, oloc]),
        )
        newIndel._indels = None
        newIndel._indel_strings = new_strings
        newIndel._indel_index = string_index
        return newIndel
```

IMPORTANT:
- `get_median_length`, `get_longest_indel`, `get_most_abundant_indel`, `get_location_bin_count`, `__radd__`, `__str__`, `__repr__` are UNCHANGED (they call the list/count getters above, which now come from arrays).
- `get_total_length` must return a Python `int` (hence `int(...)`), matching the old `sum(...)`.
- `get_mean_length` is UNCHANGED.

- [ ] **Step 4: Run tests + ruff** — `.venv/bin/python -m pytest tests/test_indel.py -q` and `.venv/bin/python -m ruff check src/lqc/indel.py tests/test_indel.py` → all pass.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/indel.py tests/test_indel.py
git commit -m "perf: numpy columnar Indel with concatenate merge"
```

---

## Task 3: `Mismatch` per-type numpy arrays

**Files:**
- Modify: `src/lqc/mismatch.py`
- Test: `tests/test_mismatch.py`

- [ ] **Step 1: Add empty/merge/order tests**

Append to `tests/test_mismatch.py`:

```python
def test_mismatch_empty_plus_nonempty():
    empty = Mismatch('e')
    full = Mismatch('f')
    full.add_mismatch('ac', 0.1)
    full.add_mismatch('gt', 0.9)
    total = sum([empty, full])
    assert total.get_total_count() == 2
    assert total.get_locations() == [0.1, 0.9]


def test_mismatch_merge_order_by_type():
    a = Mismatch('a')
    b = Mismatch('b')
    a.add_mismatch('ac', 0.1)
    a.add_mismatch('gt', 0.2)
    a.add_mismatch('ac', 0.3)
    b.add_mismatch('gt', 0.4)
    total = a + b
    # get_locations iterates types in first-seen order: ac then gt
    assert total.get_locations() == [0.1, 0.3, 0.2, 0.4]
    assert total.get_type_count() == {'ac': 2, 'gt': 2}
```

- [ ] **Step 2: run to verify they pass against old code** — `.venv/bin/python -m pytest tests/test_mismatch.py -q` (characterization).

- [ ] **Step 3: Implement**

In `src/lqc/mismatch.py`, replace `__init__`, `add_mismatch`, `_get_bin_count`, `get_location_bin_count`, `get_location_histogram`, `get_location_bin_count_by_type`, `get_locations`, `convert_complement`, and `__add__` with:

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
        self._mismatch_types = []
        self._mismatch_index = {}
        self._type_count = Counter()
        self._type_locations = defaultdict(list)
        self._type_location_arrays = None  # {type: float64 array}

    def add_mismatch(self, mismatch,
                     normalized_read_location):
        if mismatch not in self._mismatch_index:
            self._mismatch_index[mismatch] = len(self._mismatch_types)
            self._mismatch_types.append(mismatch)
        self._type_count[mismatch] += 1
        self._type_locations[mismatch].append(
            normalized_read_location
        )
        self._type_location_arrays = None

    def _finalize(self):
        if self._type_location_arrays is None:
            self._type_location_arrays = {
                ty: np.asarray(locs, dtype = np.float64)
                for ty, locs in self._type_locations.items()
            }
            self._type_locations = None
        return self._type_location_arrays

    def _get_bin_count(self, value_list, cuts):
        hist, edges = np.histogram(
            value_list, bins = cuts, density = False
        )
        bin_count = Counter()
        for i in range(len(hist)):
            if i < (len(hist) - 1):
                label = f'[{edges[i]},{edges[i+1]})'
            else:
                label = f'[{edges[i]},{edges[i+1]}]'
            bin_count[label] = hist[i]
        return bin_count

    def get_location_bin_count(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        bin_count = Counter()
        arrays = self._finalize()
        for ty in self._mismatch_types:
            bin_count += self._get_bin_count(
                arrays[ty], cuts = cuts
            )
        return bin_count

    def get_location_histogram(self, cuts = None):
        if cuts is None:
            cuts = [i / 10 for i in range(11)]
        arrays = self._finalize()
        total_hist = None
        edges = None
        for ty in self._mismatch_types:
            hist, edges = np.histogram(
                arrays[ty], bins = cuts, density = False
            )
            if total_hist is None:
                total_hist = hist.astype(np.float64)
            else:
                total_hist += hist
        if total_hist is None:
            total_hist, edges = np.histogram(
                [], bins = cuts, density = False
            )
        return edges, total_hist

    def get_location_bin_count_by_type(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        type_bin_count_dict = defaultdict(Counter)
        arrays = self._finalize()
        for ty in self._mismatch_types:
            type_bin_count_dict[ty] += self._get_bin_count(
                arrays[ty], cuts = cuts
            )
        return type_bin_count_dict

    def get_locations(self):
        arrays = self._finalize()
        locations = []
        for ty in self._mismatch_types:
            locations.extend(arrays[ty].tolist())
        return locations

    def convert_complement(self):
        newMis = type(self)(self.label)
        new_types = [
            convert_complement(a)
            for a in self._mismatch_types
        ]
        newMis._mismatch_types = new_types
        newMis._mismatch_index = {
            ty: i for i, ty in enumerate(new_types)
        }
        newMis._type_count = Counter({
            new_ty: self._type_count[old_ty]
            for old_ty, new_ty in zip(
                self._mismatch_types, new_types, strict = True
            )
        })
        arrays = self._finalize()
        newMis._type_location_arrays = {
            new_ty: arrays[old_ty].copy()
            for old_ty, new_ty in zip(
                self._mismatch_types, new_types, strict = True
            )
        }
        newMis._type_locations = None
        return newMis

    def __add__(self, other):
        assert isinstance(other, type(self)), 'wrong object to add'
        newMis = type(self)(f'{self.label} {other.label}')
        new_types = list(self._mismatch_types)
        type_index = {ty: i for i, ty in enumerate(new_types)}
        for ty in other._mismatch_types:
            if ty not in type_index:
                type_index[ty] = len(new_types)
                new_types.append(ty)
        s_arrs = self._finalize()
        o_arrs = other._finalize()
        empty = np.array([], dtype = np.float64)
        new_arrays = {}
        for ty in new_types:
            new_arrays[ty] = np.concatenate([
                s_arrs.get(ty, empty), o_arrs.get(ty, empty)
            ])
        newMis._mismatch_types = new_types
        newMis._mismatch_index = type_index
        newMis._type_count = self._type_count + other._type_count
        newMis._type_location_arrays = new_arrays
        newMis._type_locations = None
        return newMis
```

IMPORTANT:
- `get_total_count`, `get_type_count`, `__repr__`, `__str__`, `__radd__` are UNCHANGED (`get_type_count` returns `Counter(self._type_count)` as before).
- `get_type_count` returning a fresh `Counter` is kept; it must remain a plain Python `Counter`, not numpy.

- [ ] **Step 4: Run tests + ruff** — `.venv/bin/python -m pytest tests/test_mismatch.py -q` and `.venv/bin/python -m ruff check src/lqc/mismatch.py tests/test_mismatch.py` → all pass.

- [ ] **Step 5: Commit**

```bash
git add src/lqc/mismatch.py tests/test_mismatch.py
git commit -m "perf: numpy per-type Mismatch with concatenate merge"
```

---

## Task 4: `Splice` — drop `deepcopy` from `__add__`

**Files:**
- Modify: `src/lqc/splice.py`
- Test: `tests/test_splice.py`

- [ ] **Step 1: Add a test pinning merge behavior**

Append to `tests/test_splice.py`:

```python
def test_add_counts_accumulate_without_aliasing():
    a = Splice('a')
    b = Splice('b')
    a.add_splice_pair('gt-ag')
    b.add_splice_pair('ct-ac')
    b.add_splice_pair('gt-ag')
    total = a + b
    assert total.get_splice_pair_count_dict() == {'gt-ag': 2, 'ct-ac': 1}
    # the operands are untouched (no aliasing)
    assert a.get_splice_pair_count_dict() == {'gt-ag': 1}
    assert b.get_splice_pair_count_dict() == {'ct-ac': 1, 'gt-ag': 1}
```

- [ ] **Step 2: run to verify it passes against old code** — `.venv/bin/python -m pytest tests/test_splice.py -q` (characterization).

- [ ] **Step 3: Implement**

In `src/lqc/splice.py`, remove the now-unused `from copy import deepcopy` and replace `__add__` with:

```python
    def __add__(self, other):
        assert isinstance(other, type(self)),\
            'wrong object to add'
        new_dict = (
            self.get_splice_pair_count_dict() +
            other.get_splice_pair_count_dict()
        )
        newSp = type(self)(
            f'{self.label} {other.label}'
        )
        newSp.add_splice_pair_count_dict(new_dict)
        return newSp
```

(`Counter.__add__` already returns a new Counter; values are `int`, keys are `str`, so no deep copy is needed.)

- [ ] **Step 4: Run tests + ruff** — `.venv/bin/python -m pytest tests/test_splice.py -q` and `.venv/bin/python -m ruff check src/lqc/splice.py tests/test_splice.py` → all pass (no unused import).

- [ ] **Step 5: Commit**

```bash
git add src/lqc/splice.py tests/test_splice.py
git commit -m "perf: drop deepcopy from Splice.__add__"
```

---

## Task 5: Full-suite + byte-parity + memory verification

**Files:** none (verification only).

- [ ] **Step 1: Full test suite** — `.venv/bin/python -m pytest -q` → all pass (133 + the ~8 new tests ≈ **141+**; confirm every test passes and note the exact total). `.venv/bin/python -m ruff check src tests` → clean.

- [ ] **Step 2: chr22 byte-parity** — `.venv/bin/lqc -b tmp/full_bam/ENCFF417VHJ.sorted.bam -c chr22 -o tmp/m2/chr22 -t 1 --output-cs`, then `cmp` all `table/*.txt`, `read.cs`, `LQC_report.html`, and `fig/*.png` against `tmp/verify/chr22`. Expected: every deterministic artifact `OK` (PDFs excluded — creation timestamp).

- [ ] **Step 3: full-genome run + peak RSS** — `bash tmp/run_lqc_full.sh tmp/full_bam/ENCFF417VHJ.sorted.bam tmp/m2/full 4`. Expected: exits 0 and `All done!` (this is THE P2 fix — previously OOM-killed at 8.4 GB in the `sum()`/`'Total'` step). Read `tmp/m2/full.time.txt` "Maximum resident set size" → target **≤ ~2 GB**.

- [ ] **Step 4: Record results in this plan** — note byte-parity outcome, new peak RSS, and wall time. If byte-parity fails or the run still OOMs, STOP and investigate before proceeding to M3.

---

## Self-Review (run before handoff)

- **Spec coverage:** §4.3 — ReadStat int32 (n,7) + memoized N/L (Task 1), Indel three parallel arrays + direct histogram + concat merge (Task 2), Mismatch per-type float64 arrays + concat (Task 3), Splice deepcopy removal (Task 4), float64 precision + int32 counts (all tasks), getters return `.tolist()` so report_table/figure stay zero-change (all tasks), merge ordering preserved (concat keeps order; Task 5 byte-parity). Covered.
- **Placeholder scan:** complete code in every step; no TBD/TODO. Merged objects are final (array-only) — the plan's tests only exercise add-then-merge and merge ordering, never add-after-merge.
- **Type consistency:** `_finalize()` return types are consistent per class (`(n,7)` int32; `(iidx, ilen, loc)`; `dict[str, float64 array]`); `get_*` boundaries return Python `int`/`float`/`list`; `get_type_count` returns `Counter`. Consistent.
- **Key risk to re-verify in review:** empty-object `reshape(-1, N)` (Task 1/2), `__radd__` returning `self` for `0 + x` (unchanged), and `Mismatch.get_locations`/`get_type_count` ordering (first-seen type order must match old `_mismatch_types` order — preserved because `_mismatch_types` is still a list).
---

## Results (recorded post-implementation)

- **Branch:** `perf/m2-numpy-columnar-storage`; commits `1981c06` (ReadStat), `2613cef` (zero-length divide guard), `9f171e0` (Indel), `7c1176b` (Mismatch), `9e871e6` (Splice). Merged to `main` via `--ff-only`.
- **Tests/lint:** full suite **140 passed** (133 baseline + 7 new), `.venv/bin/python -m ruff check src tests` clean.
- **chr22 byte-parity:** 70/70 deterministic artifacts match `tmp/verify/chr22` (7 `table/*.txt`, `read.cs`, `LQC_report.html`, 61 `fig/*.png`).
- **Full-genome run** (`-t 4`, `tmp/m2/full`): **exit 0 — completes** (previously SIGKILL at >8.4 GB). Wall **8:26.86**, stat 2.6 min / reduce+Total 19 s / figures 5.3 min. Peak RSS **6,658,812 KB ≈ 6.35 GB**.
- **Honest finding:** the merge/`'Total'` deep-copy OOM is fixed (the run now completes), and peak RSS dropped from the fatal >8.4 GB to 6.35 GB — but §6's "≤ ~2 GB working set" target is NOT met. The residual peak lives in the stat phase: worker-side boxed accumulation plus the `mp.Pool.map` result transfer of all chunks into the parent before `reduce_blocks_to_contigs` finalizes them (lazy `_finalize()` only converts at merge). Finalizing inside `stat_region` before returning (numpy pickle) is a candidate follow-up, outside this milestone's merge-scope.
