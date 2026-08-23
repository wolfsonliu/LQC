# Performance & Parallelization Design

> **Status:** approved for implementation
> **Branch/context:** post `src/lqc` migration (src layout)
> **Baseline:** `tmp/data/ENCFF417VHJ.chr22.sorted.bam` — 1 contig (chr22), 355,250 reads, ~4m23s with `-t 4` (tables+figures+HTML only)

## Problem Statement

The tool is effectively single-threaded on its dominant input shape, despite a `-t/--thread` flag.

1. **Parallelism is per-contig, not per-read.** `mp.Pool(thread).map(stat_bam, contigs)` (`src/lqc/cli.py:378-389`) fans work only over contigs. A single-contig BAM (e.g. chr22) therefore schedules exactly one task → one worker. Evidence: an earlier `/usr/bin/time -v` run with `-t 4` reported ~98% CPU (≈1 core). The per-read loop `for read in bam.fetch(contig)` (`src/lqc/stat.py:34`) is pure-Python and is the critical path (~3 minutes of the run).
2. **Per-read overhead** in the loop is dominated by: cs-string parsing (`CS.from_cs_tag_string` → `cs_to_list` + `_compute_counts`, two passes over the element list), lazy derivation of normalized read coordinates (`get_normalized_read_location`, cached but built once per read), three sequential `_get_marks` filters over that coordinate list (`get_insertions/get_deletions/get_mismatches` at `src/lqc/stat.py:76-101`), and per-event string construction on reverse-strand reads (`convert_reverse_complement` / `convert_complement` per event).
3. **`_compute_counts` uses two regex substitutions per intron** (`re.sub('[0-9]+', ...)` twice, `src/lqc/cs.py:375-382`) to recover donor/acceptor, though both are fixed-width slices of the intron token (`value[:2]`, `value[-2:]`). With ~800k intron events in the baseline this is measurable.
4. **`--output-cs` re-parses the whole BAM a second time.** `write_readcs` (`src/lqc/utils.py:76-137`) re-opens the BAM and rebuilds a `CS` object per read, then writes a ~3.5 GB, 6-column `read.cs` through a per-line `writelines` generator that allocates a list and a format string per element.

## Proposed Solution

Three stages, all output-preserving (integer counts are summed identically):

- **Stage 1 — serial micro-optimizations** (low risk, measured first).
- **Stage 2 — data-parallel map/reduce** ("prefetch → chunk → `Pool.map` → reduce") so `-t N` scales independent of contig count.
- **Stage 3 — single-pass `--output-cs`** (fold into the stat pass, remove the redundant re-parse).

## Detailed Design

### Stage 1 — Serial micro-optimizations

| # | Change | Location | Rationale |
|---|---|---|---|
| 1.1 | Replace `next(a[1] for a in read.tags if a[0]=='cs')` with `read.get_tag('cs')` (and `get_tag('MD')`) | `src/lqc/stat.py:38-41`; `src/lqc/utils.py:99-102`, `:117-120` | `read.tags` materializes and scans the whole tag tuple per read; `get_tag` is a direct lookup. |
| 1.2 | In `_compute_counts`, derive intron donor/acceptor/pair via slicing instead of two regex subs | `src/lqc/cs.py:375-382` | `value[:2]`, `value[-2:]`, and `f"{value[:2]}-{value[-2:]}"` are exact; drops 2 regex calls per intron (~1.6M calls on baseline). |
| 1.3 | Compute the three mark-filtered lists in one pass (group by mark) instead of three `_get_marks` scans; expose cached `get_insertions/get_deletions/get_mismatches` results | `src/lqc/cs.py:591-619` | The normalized coordinate list is computed once; filtering it three times is redundant. |
| 1.4 | Replace `convert_complement`'s per-character dict-lookup `join` with `str.translate` (precomputed 256-byte table over the same `acgtn-` alphabet), reverse via slicing | `src/lqc/utils.py:6-14` | Per-event string build over millions of events; `translate` is a single C-level pass and is byte-identical to today. |
| 1.5 | Remove dead `_get_read_length` (duplicate of the cached `get_read_length`) | `src/lqc/cs.py:481-494` | Clarity; it re-implements the cached counter. |

**Stage 1 exit criterion:** identical table outputs on the chr22 baseline and the `tests/` fixtures; a recorded wall-clock delta.

### Stage 2 — Data-parallel map/reduce

Replace the per-contig pool with per-read-chunk parallelism:

1. **Prefetch (main process, single pass).** Iterate `bam.fetch()` once, extracting only lightweight fields per read: `(contig, start_pos, strand, query_length, query_name, cs_or_md_tag_string)`. This is record-decoding plus `get_tag` — far cheaper than `CS` parsing. (For the MD path, the worker is responsible for `genome.fetch`, which must stay inside the worker to avoid serializing reference windows.)
2. **Chunk.** Split the ordered prefetch list into `threads` contiguous blocks (one per worker). Blocks preserve BAM order, which keeps every order-sensitive output (e.g. length lists feeding cumulative-length plots) identical to serial.
3. **Map.** `Pool(threads).map(worker, blocks)`, where `worker` is `stat_element_from_bam_by_contig` refactored to accept a read tuple/block instead of opening the BAM itself for the `cs` path. Each worker returns the existing `(ReadStat, Indel, Indel, Mismatch, Splice)` 5-tuple.
4. **Reduce.** Sum the per-chunk tuples with the existing `sum()` / `__add__` machinery (`src/lqc/cli.py:405-425` is reused verbatim). Chunk count equals thread count, so the O(n²) list concat in `ReadStat.__add__` is bounded and small in practice.

**Correctness/determinism argument.** Every emitted metric is either an integer count (summed identically by chunking) or a statistic derived from the *merged* per-read list (mean/median/N50/L50), which is order-preserving because blocks are ordered and concatenated in order. Fractional columns in the summary tables are ratios of integer sums, so they reproduce bit-for-bit. Figure data comes from the same merged objects.

**Behavior change.** `-t/--thread` default becomes `min(os.cpu_count(), 4)` (from `1`). The flag keeps the same meaning; documented in `README.md` usage and `src/lqc/cli.py` help.

**Memory note.** The prefetch list is ≈355k × a small tuple (~tens of MB) for the cs-path baseline, and it lives until `main()` returns. For the MD path each record additionally carries the full `query_sequence`, so whole-genome MD input needs proportionally more memory.

### Stage 3 — Single-pass `--output-cs`

1. Extend each worker to also emit its chunk's `read.cs` lines (the `get_contig_position()` records for the reads it parsed) into a per-chunk temp file under the output dir, in BAM order.
2. Main process concatenates the chunk temp files **in chunk order** into `read.cs`, then removes the temp files. This preserves the current `read.cs` schema, per-contig row order, and row contents for the analyzed contigs, while eliminating the second full parse (the old `write_readcs` iterated the entire BAM).
3. Buffer writes per chunk (write whole lines in memory per chunk, or use a buffered writer) instead of the current per-element `writelines` generator + per-line list/format allocation.
4. Optional `--readcs-gzip` flag to gzip the (large) `read.cs`; default stays uncompressed to preserve byte-compatibility for existing consumers.

## User Experience / Behavior Considerations

- No CLI flag is removed; `-t`, `-c`, `--output-cs`, `--output-pickle` keep their meaning.
- `read.cs` schema is unchanged; it now contains only the analyzed contigs in requested order (the CLI warns when references are omitted), instead of the old whole-BAM dump.
- Empty-contig and header-present-but-unmapped cases are handled identically to today (the existing `read_count > 0` filter at `src/lqc/cli.py:396-403` still applies after reduce).
- Logging keeps the existing start/finish lines; add one debug line per stage for diagnosis.

## Implementation Notes / Trade-offs

- **Why "prefetch → chunk → map" over coordinate-window `fetch`:** coordinate windows risk double-counting reads that span window boundaries and suffer load imbalance; the prefetch approach makes determinism and single-parsing trivially correct.
- **Why chunk count = thread count:** keeps the reduce step O(threads) and preserves order; more chunks would only add pickle/reduce overhead without more parallelism (CPU-bound).
- The MD path (`genome.fetch` inside worker) is exercised by existing fixtures; no FASTA serialization is introduced.
- `ReadStat.__add__`'s `deepcopy(self._reads + other._reads)` is left as-is for now; moving `ReadStat` to columnar/numpy storage is explicitly **out of scope** for this design.

## Out of Scope (follow-up, not this change)

- Default contig list derived from `bam.references` instead of the hardcoded `chr1..chrY` (`src/lqc/cli.py:120-127`) — the 23 "contig not in BAM" warnings.
- Optional PDF output, skipping redundant single-contig `Total`/contig figures, lazy `matplotlib` import.
- Progress reporting and per-contig timing.
- `ReadStat` columnar storage / O(n²) `__add__` elimination for many-contig genomes.

## Success Criteria

1. **Bit-identical outputs:** summary tables (`table/*.txt`) and `read.cs` for the chr22 baseline and all `tests/` fixtures are byte-identical to the pre-change run.
2. **Determinism:** two consecutive runs of the parallel build produce byte-identical `table/` and `read.cs`.
3. **Speedup:** with `-t 4`, the chr22 baseline wall-clock drops measurably (target ≥2.5× for the element-statistics phase; the exact measured ratio is reported alongside the change).
4. **Quality gate:** `uv run ruff check` passes and the full `uv run pytest` suite is green with **no** behavior change to existing assertions beyond the documented `-t` default.
5. **Back-compat:** `uv run lqc --version`, the `-c`/`--output-pickle` flag paths, and the MD-tag fallback all still work.

## Dependencies

- None external. Relies on the existing `__add__`/`__radd__` reduce and `sum()` composition already used by `src/lqc/cli.py:405-425`.
- Order of work: Stage 1 → Stage 2 → Stage 3 (each independently shippable and verifiable).

## Open Questions

None unresolved at design time. Decisions made above: parallel approach = prefetch/chunk/map/reduce; chunk count = thread count; `-t` default = `min(cpu_count(), 4)`; `--readcs-gzip` optional and off by default.