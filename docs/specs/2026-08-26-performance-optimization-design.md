# LQC Whole-Genome Performance Optimization — Design

- **Date:** 2026-08-26
- **Status:** Approved (design level; implementation plans follow in `docs/plans/`)
- **Scope:** performance/memory optimization for whole-genome long-read RNA-seq
  inputs (24 references, ~1.6 M reads). No user-facing feature changes.
- **Related:** profiling harness `tmp/profile_lqc.py`, run wrapper `tmp/run_lqc_full.sh`.

## 1. Problem Statement

LQC's current read-processing path does not scale to a whole-genome, `cs`-tagged
PacBio CCS BAM. The motivating input is `tmp/full_bam/ENCFF417VHJ.sorted.bam`:

| property | value |
| --- | --- |
| size | 561 MB (sorted, `.bai` indexed) |
| references | 24 (`chr1`..`chr22`, `chrX`, `chrY`) |
| mapped reads | 1,574,822 |
| `chr1` (largest) | 406,210 reads |
| `chrY` (smallest) | 125 |
| tags | `cs` only (PacBio CCS, predominantly reverse-strand) |

Two measurements establish the problem empirically.

### 1.1 CPU: read processing is dominated by redundant work

A single-contig, single-process profile of `stat_element_from_bam_by_contig` on
`chr22` (43,339 reads) gives **15.62 s wall / 221 MB RSS**, with cumulative time:

| hotspot | cum time | share |
| --- | --- | --- |
| `get_indel_mismatches` → location chain (`_get_coordinate_list` → `get_normalized_read_location` → `get_read_location` → `_get_element_read_location`) | ~5.5 s | **35%** |
| `from_cs_tag_string` → `cs_to_list` (2× `re.sub` + `split` + boxed lists) | ~4.3 s | **27%** |
| `CS._compute_counts` | ~2.0 s | 13% |
| `re.sub` primitive (265,998 calls incl. the per-intron `~` length lambda) | ~1.9 s | 12% |

The single largest cost is computing a full normalized `[0,1]` read location for
**every** `cs` element — including every `:` match run — through four stacked
list passes, only to keep the `+`/`-`/`*` entries (see §4.2).

### 1.2 Memory: the full-genome run exhausts RAM

`prefetch_records` first materializes every read into a `ReadRecord` in one list,
single-threaded, before any parallelism; `cli.main` then re-serializes the entire
task list through `mp.Pool.map`. Measured prefetch cost is ~0.86 KB/read
(`chr1` = 406,210 reads → 422 MB RSS), i.e. ≥ 1.35 GB for the full set before the
pool pass, fork workers, and downstream `sum()` deep-copies even run.

A full `-t 4` run on an 11.8 GB machine climbed past **7 GB RSS** and made no
progress past the prefetch/pool phase in ~8 minutes; it was killed and memory
recovered to 3.9 GB. The current architecture cannot process this input on a
modest-RAM machine.

### 1.3 Root-cause summary (P0–P3)

- **P0** — global materialization + IPC serialization of all reads (`prefetch_records` + `Pool.map` pickle).
- **P1** — redundant normalized-location computation + multi-regex `cs` tokenization (CPU).
- **P2** — per-read (`ReadStat._reads`) and per-event (`Indel._indels`, `Mismatch._type_locations`) boxed Python lists, and `__add__`/`sum()` deep-copies them repeatedly (memory + merge time).
- **P3** — report-phase polish: per-contig figures for 24 contigs, plus micro-fixes.

## 2. Goals & Non-Goals

### 2.1 Goals

1. Whole-genome input completes on a machine with a few GB of available RAM
   (target: comfortably under 2× input reading working set; see §6).
2. Read-processing CPU reduced ~2× via P1 (single-pass parse + fused location).
3. Outputs stay **byte-identical** to today's implementation across `table/*.txt`,
   `fig/*`, `LQC_report.html`, and `read.cs`.
4. Deterministic: results independent of `-t/--thread` (existing tests guard this).

### 2.2 Non-goals

- No change to CLI argument semantics, summary-table column names/formats, HTML
  template structure, or output file layout.
- No change to the `cs`/`MD` precedence rule; `MD` inputs keep working (unchanged,
  see §4.2 "Scope").
- No new runtime dependencies (`numpy` is already a dependency) and no change to
  the multiprocessing framework (still `multiprocessing`).
- No changes to `CS`'s public methods consumed by `--output-cs` or library callers.

## 3. Proposed Solution (Roadmap)

Four milestones, applied in dependency order `M0 → M1 → M2 → M3`. Each milestone
keeps outputs byte-identical and is independently verifiable (§5).

| milestone | maps to | what changes | why first |
| --- | --- | --- | --- |
| M0 | P0 | worker-side `fetch`; delete global prefetch + pool pickle | unblocks whole-genome runs (hard blocker) |
| M1 | P1 | single-pass `cs` scan + fused count/location | ~2× stat CPU, low risk, same data model |
| M2 | P2 | numpy columnar storage for `ReadStat`/`Indel`/`Mismatch`; cheap merges | memory + merge time; larger surface, done later |
| M3 | P3 | per-contig figure policy + micro-fixes | polish, not blocking |

M1 depends on M0 solely for comfortable end-to-end validation on the full
dataset; the two are logically independent. M2 depends on M0/M1 being stable
because it touches the accumulator classes and their `__add__` that M0's reduce
step relies on.

## 4. Detailed Design

### 4.1 M0 — Worker-side fetch (remove global materialization + IPC serialization)

**Current flow** (`src/lqc/cli.py` → `src/lqc/stat.py`):

```
prefetch_records(bam, contigs, method)   # builds [(contig, [ReadRecord, ...]), ...], all reads
  └─ for contig in contigs: [record_from_read(r) for r in bam.fetch(contig)]
_chunk(records, threads)                 # per-contig slices -> tasks
mp.Pool(threads).map(partial(stat_records, ...), tasks)   # pickles ALL tasks
reduce_blocks_to_contigs(...)
```

**Target flow:**

```
plan_tasks(bam, contigs, threads)        # [(index, contig, start, end), ...] only
mp.Pool(threads).map(partial(stat_region, ...), tasks)
  └─ stat_region: opens its own AlignmentFile, fetch(contig, start, end),
     builds transient ReadRecords, accumulates a StatBlock
reduce_blocks_to_contigs(...)            # unchanged
```

- `ReadRecord` becomes a **transient** object inside a worker; it is no longer
  collected into a global list, and is never pickled across the pool boundary.
- Each worker opens its own read-only `AlignmentFile` handle (open-per-worker is
  cheap and avoids sharing a handle across processes).

**Chunking = balanced coordinate windows (decision: option (b)).**

- `task[0]` `index` preserves global order; `(contig, start, end)` is a
  half-open `[start, end)` window on the reference.
- Per-contig mapped counts come from `AlignmentFile.get_index_statistics()`
  (same data as `samtools idxstats`). A target chunk size (~10k reads) is chosen;
  each contig's coordinate span is split into roughly
  `ceil(mapped / target)` windows of ~equal read count. Contigs smaller than the
  target become a single task.
- Tasks are emitted in **requested-contig order, then ascending coordinate**
  within a contig, so `bam.fetch(contig, start, end)` returns reads in the same
  order the current `prefetch_records` produces. This preserves the current
  read.cs and table row ordering and keeps the thread-determinism tests green.

**`read.cs` ordering (already-audited invariant).**

`stat_region` writes its `.readcs-NNNN.tmp` as today; `cli.main` merges in
`block_results` order. Because `index` ordering equals current traversal order,
the merged `read.cs` is byte-identical. No change to the merge code.

**Micro-fix folded in:** `record_from_read` computes `query_length =
read.query_length` instead of `len(read.query_sequence)` (the latter materializes
the full sequence string purely to count it; measured 1.6× on that op and it
allocates the sequence per read).

**Load balancing note.** `chr1` (406k) vs `chrY` (125) differ by 3000×; any
per-contig tasking leaves `chr1` as a straggler. Coordinate windows let a large
contig occupy several workers; the pool's `map` still applies at the end, so the
last window's shortfall is bounded by one chunk per worker — acceptable.

### 4.2 M1 — Single-pass `cs` scan + fused count/location

**Current call chain for one read (`process_record`):**

```
_build_cs
  CS.from_cs_tag_string → cs_to_list (2× re.sub + 2× split; boxed [low,high,mark,value] lists)
  CS.__init__ → _compute_counts (second full pass over the boxed list)
readstat.add_read(<counts from getters>)
cs.get_indel_mismatches('normalized_read')          # ← the 35% cost
  _get_coordinate_list → get_normalized_read_location → get_read_location
    → _get_element_read_location (builds read coords for EVERY element)
```

The last step computes `[low, high]` read coordinates and divides by `read_len`
for **all** elements (including `:` match runs and `~` introns), wraps them in
three more lists, and then `get_indel_mismatches` keeps only `+`/`-`/`*`.
`process_record` itself uses only the normalized `low` and the element `value`
(sequence).

**Target:** one pass over the `cs` string produces (a) the parsed element list,
(b) the precomputed counts, and (c) the raw read coordinates of `+`/`-`/`*`
elements; normalization (divide + strand mirror) happens once afterwards, on
indels/mismatches only.

Concretely, in `CS.__init__`:

1. **Tokenizer (replaces `cs_to_list`).** Manual single-pass scan: the mark set
   `{:,*,+,-,~}` is disjoint from the value character set (digits + `a-z`), so a
   state-free loop advances a cursor, reads one mark then one maximal digit/lower
   run. No `re`, no intermediate full-string `re.sub`/`split`. Each element is
   recorded as a `(mark, value)` tuple. The intron `~` length is derived by
   reading the digits directly (replacing the per-intron
   `re.sub('[a-z]', '', value)` lambda).
2. **Same loop accumulates the counts** currently computed by `_compute_counts`
   (read length, intron/mismatch/insertion/deletion counts and lengths,
   splice-pair/site counters, mismatch-type counter) into the same instance
   attributes, so all existing getters keep working unchanged.
3. **Same loop records, for `+`/`-`/`*` only**, the raw read-space coordinates
   `(read_low, read_high, mark, value)` into a scratch list. `read_low/read_high`
   are the running read position before/after the element (per the existing
   advance rules: `:` and `+` and `*` advance the read position by their read
   length; `-` and `~` do not).
4. **After the loop**, `read_len` is known. For each recorded `+`/`-`/`*`:
   - strand `+`: `low_n = read_low / read_len`, `high_n = read_high / read_len`;
   - strand `-`: mirror exactly as today — `low_n = (read_len - read_high) / read_len`,
     `high_n = (read_len - read_low) / read_len`.
   Store the result once in a **new private attribute**
   `self._indel_mismatch_normalized` as a `[low_n, high_n, mark, value]` list for
   the `+`/`-`/`*` elements only.

**Public-API preservation.** `get_indel_mismatches(coordinate='normalized_read')`
serves from `self._indel_mismatch_normalized` (no recompute). The existing
`get_normalized_read_location()` keeps its current lazy **full-list** behavior
(all elements, including `:` and `~`) for library callers, so nothing that
depends on it changes. Every other `CS` public method (`get_relative_position`,
`get_contig_position`, `get_read_location`, `get_cs_tag_string`,
`count_*_in_read_location_bin`, …) keeps its existing lazy behavior, so
`--output-cs` (`get_contig_position`) and library callers are untouched. The
non-`normalized_read` coordinate variants of `get_indel_mismatches` continue down
the existing lazy path.

**Float semantics.** The float division and strand mirror are computed with the
same expressions as the current code (`int / int` in Python 3), so the produced
`float64` values are bit-identical; this is what makes the assertion-byte-parity
test in §5 meaningful.

**Scope.** This milestone optimizes the `cs`-tag path only. The `MD`/CIGAR path
(`convert_cigar_md_to_cs_list` → `CS.from_cigar_string`) is left unchanged and is
not a dependency of this dataset; optimizing it is tracked separately.

### 4.3 M2 — numpy columnar storage, cheap merges, cached aggregates

Goal: eliminate the boxed per-read / per-event lists and the `__add__`/`sum()`
deep-copies that create the "multi-GB event lists" flagged in `cli.py`.

**Common contract (applies to all four classes):**

- Accumulation-time storage may still append to plain Python lists (cheap, and
  chunk-bounded after M0); conversion to numpy happens lazily via a private
  `_finalize()`/`_as_array()` invoked at merge/finalize time.
- `__add__`/`sum()` merges use `np.concatenate` (Indel/Mismatch rebuild a string
  index once), never `deepcopy`.
- Every public getter keeps its signature and return type. Where current code
  returns `list[int]`/`list[float]`, callers that rely on Python-list semantics
  (e.g. `report_figure` doing `sorted(...)` or list concatenation) receive
  `.tolist()`; `report_table`/`report_figure` remain **zero-change**.
- Labels stay strings; `sum()` continues to concatenate labels and the `'Total'`
  sentinel is applied by the existing shallow-copy `relabel` in `cli.main`.

**`ReadStat`** (`src/lqc/readstat.py`)

- `_reads` → `np.int32` array of shape `(n, 7)` with columns
  `[length, insertion, deletion, mismatch, intron, mapping_quality, aligned_length]`.
- Derived metrics (`get_length_NL(50)`, medians, means, per-base rates) are
  computed on the array and **memoized once**; `get_length_NL` currently calls
  `get_lengths()` + `sorted(...)` on every call, and is invoked separately for N50
  and L50 by the table and figure layers.

**`Indel`** (`src/lqc/indel.py`)

- `_indels` → three parallel arrays: `_indel_iidx` (`int32`), `_indel_len`
  (`int32`), `_indel_loc` (`float64`). `_indel_strings` (a bounded list of
  distinct sequences) and `_indel_index` are retained.
- `get_location_histogram` drops `np.fromiter(...)` and calls
  `np.histogram(self._indel_loc, bins=cuts)` directly.
- `__add__` concatenates the arrays and re-maps the string index once (no deepcopy).

**`Mismatch`** (`src/lqc/mismatch.py`)

- `_type_locations` (`defaultdict(list)` of per-type floats) → per-type numpy
  `float64` arrays (or one concatenated array plus per-type offsets). `_type_count`
  (a `Counter`) is retained.
- `get_location_histogram` vectors over each type array with `np.histogram`.
- `__add__` merges per-type arrays by concatenation (no per-type `list` rebuild).

**`Splice`** (`src/lqc/splice.py`)

- Already compact (a `Counter` of splice-pair strings). Replace `__add__`'s
  `deepcopy` of the Counter with additive `Counter.__add__`/`update` (a deepcopy
  of the underlying dict of a bounded Counter is unnecessary).

**Precision decision (approved): `float64`.** Normalized locations are `float64`
to preserve bit-identical histogram output. Using `float32` could shift values
across the 0.25/0.5/0.75 bin edges and break byte-parity. The count/length
columns stay `int32` (max read length ≪ 2^31).

**Merge ordering (adjacent to M0's reduce).** `reduce_blocks_to_contigs` and the
`sum(...)` for `'Total'` must concatenate in the same contig/chunk order as today;
this is already the case and must be re-verified byte-for-byte after M2.

### 4.4 M3 — Report-phase polish (lower priority)

- **Per-contig figures.** `cli.generate_multiple_figs` emits per-contig + total +
  aggregate figures for ~15 plot types; for 24 contigs that is ≈ 780 PNG/PDF
  files. The cost has not yet been measured end-to-end (the full run dies in M0).
  Options, to be decided in the M3 implementation plan after measuring:
  1. keep current behavior (no change);
  2. add an opt-in/opt-out for per-contig figures (default off for ≥ N contigs);
  3. render per-contig figures in a thread pool (matplotlib under `Agg`, no
     shared state).
- **Micro-fixes:** `convert_reverse_complement` uses a dedicated
  reverse-complement translate table (one `translate` instead of
  translate + reverse slice); `convert_complement`/`convert_reverse_complement`
  remain pure helpers.

## 5. Verification Strategy

Each milestone must pass, before being considered done:

1. `uv run ruff check` (via a sandbox-compatible invocation, e.g.
   `.venv/bin/python -m ruff check`) and the full pytest suite.
2. **Byte-parity** against the pre-change implementation:
   - On the small committed fixture (`tests/data`): run `table/*.txt`,
     `LQC_report.html`, and `read.cs` and `diff` them byte-for-byte.
   - On a real contig (e.g. `chr22` of the full BAM): diff `table/`, `fig/` (PNG
     byte-identical), `LQC_report.html`, and `read.cs`.
   - On the full BAM (post-M0, it should actually complete): spot-check or full
     diff of `table/` and `read.cs`.
3. **Thread determinism:** existing `test_tables_identical_across_thread_counts`
   and `test_output_cs_identical_across_thread_counts` stay green (they already
   compare `mapping.txt` and `splice_all.txt` byte-for-byte across `-t 1/-t 2`).
4. **Performance regression:** re-run `tmp/profile_lqc.py` (single-contig) and
   `tmp/run_lqc_full.sh` (full, `-t 4`) and record wall + peak RSS; compare to the
   15.62 s / 221 MB (chr22) and >7 GB (full, killed) baselines.

The profiling harness (`tmp/profile_lqc.py`) and run wrapper (`tmp/run_lqc_full.sh`)
are gitignored; keep them as the reproducible measurement tooling.

## 6. Success Criteria

- **Correctness:** byte-identical `table/`, `fig/`, `LQC_report.html`, `read.cs`
  vs. current implementation on fixture + `chr22` + full BAM (the non-negotiable bar).
- **Memory:** full-genome run completes with peak RSS comfortably within a few GB
  (target **≤ ~2 GB** working set, vs. the current >7 GB that OOMs); `ReadStat`
  no longer holds 1.57 M boxed lists and `Indel`/`Mismatch` no longer hold
  ~20 M boxed event lists.
- **CPU:** read-processing throughput roughly **2×** on the `cs` path
  (P1 removes the ~35% location recompute and most of the ~27% tokenization cost);
  merge/`Total` aggregation no longer performs deep-copy-of-growing-list work.
- **Determinism:** `-t 1` ≡ `-t 4` for all outputs (existing tests).

## 7. Open Questions (each resolved before its milestone's implementation plan)

1. **Chunk target size / load balancing (M0).** Default ~10k reads; calibrate
   against `get_index_statistics()` so the largest contig's last window does not
   straggle. (Leaning: 8–16k reads per window.)
2. **`get_index_statistics()` fidelity (M0).** Confirm its per-reference counts
   match `samtools idxstats` and that `fetch(contig, start, end)` yields reads in
   deterministic, index order; fall back to `samtools idxstats` parsing if not.
3. **M2 array-reduction boundary.** Whether `ReadStat` getters return `.tolist()`
   or numpy views is finalized during M2 by auditing every consumer; the design
   default is `.tolist()` at the boundary to keep `report_table`/`report_figure`
   unchanged.
4. **M3 figure policy.** Decide among keep / gated / parallel after measuring the
   per-contig figure cost post-M0; not yet measured.
5. **`MD`-path tokenizer.** Out of scope for this design; if the same single-pass
   approach is wanted for `convert_cigar_md_to_cs_list`, it becomes a follow-up.

## 8. Constraints Trace (AGENTS.md)

- Python 3.11+ — numpy usage is stdlib-compatible; no newer syntax.
- `cs` wins over `MD` — detection/precedence untouched (`utils.check_bam_with_cs_or_md`).
- Sorted + indexed input required — M0 relies on `fetch(contig, start, end)` (requires index).
- `stat_element_from_bam_by_contig(..., method in ['cs','MD','both'])` — unchanged public wrapper.
- `label` is `str`; `'Total'` reserved — M2 merge keeps the existing `relabel`.
- Summary-table column names are a contract — tables/figure layer untouched.
- 0-based half-open coordinates — M0 windows and M1 positions follow the existing convention.
- matplotlib `Agg` — M3 figure work must not introduce a GUI backend.
- GPLv3+ — no new third-party code; numpy already a dependency.
- `ruff check` — all changes must pass the frozen ruleset.

## 9. Out of Scope

- Changing `-t/--thread` semantics, summary-table schema, HTML template, or output layout.
- Migrating off `multiprocessing` (e.g. to dask/ray), or adding Cython/Rust extensions.
- Optimizing the `MD`/CIGAR→`cs` conversion path.
- Any change to `CS` public API or `--output-cs` file format.