# LQC Mapping Metrics & Output Polish (Design)

Date: 2026-08-23 · Branch: `improve` · Status: approved (Approach A: extend `ReadStat`; scope "all-in": #13 + bundled P0/P1 polish + P0 #2 aligned-base normalization)

## Problem Statement

The report has no view of **whether and how well each read mapped**. A probe of the canonical
`ENCFF417VHJ.chr22` BAM (355250 primary alignments) shows a large, currently invisible QC
signal:

- MAPQ is coarse for long reads (minimap2): `0`→4,435 and `1`→140,841 (~40%), 61 distinct
  values.
- **Aligned fraction** (`query_alignment_length / query_length`) has mean 0.687, median 0.840;
  194,381 reads (~54.7%) have <90% of their query aligned (soft-clipped/unmapped tails).

Two consequences follow. First, those metrics should be collected and reported. Second, the
existing `total_base` is the *full query* length, so the per-kb error rates are diluted by the
~31% unaligned portion — the root of the `mean_X_per_read_per_kb` vs `X_per_kb` ambiguity
(P0 #2). This design adds the mapping metrics end-to-end and, in the same pass, fixes the
output-quality issues that were catalogued in the brainstorming review so the new tables and
figures land on clean conventions instead of inheriting the old ones.

## Proposed Solution

Adopt **Approach A**: extend the existing `ReadStat` object (rather than a new stat object) to
carry per-read `mapping_quality` and `aligned_length`. Emit a new `table/mapping.txt`, three
new figures, and a new HTML "Mapping" section. In the same change, apply the agreed polish:
(1) unify float formatting, (2) unify the color palette, (3) collapse `splice.txt` to the four
canonical categories, (4) fix `mismatch.txt` column ordering + add a per-bin total, (5) make
the indel-length histogram readable, and (6) resolve the per-kb normalization ambiguity by
renaming `X_per_kb`→`X_per_query_kb` and adding `X_per_aligned_kb`.

## Detailed Design

### 1. Data collection (`stat.py`)

- Extend `ReadRecord` with two required `int` fields, placed after `query_length` (all
  following fields already have defaults):
  - `mapping_quality`
  - `aligned_length` = `read.query_alignment_length` (aligned query bases from CIGAR; excludes
    introns and soft/hard clips).
- `record_from_read` sets both for the `cs` and the `MD` branches identically
  (`read.mapping_quality`, `read.query_alignment_length`). No cs/MD-specific logic: both are
  CIGAR-derived and available on every alignment.

### 2. `ReadStat` extensions (`readstat.py`)

- `__init__`: add `self._total_aligned_base = 0`.
- `add_read(self, length, insertion, deletion, mismatch, intron, mapping_quality=0,
  aligned_length=0)`:
  - `_reads` row grows from 5 to 7 elements:
    `[length, insertion, deletion, mismatch, intron, mapping_quality, aligned_length]`.
  - new keyword params default to `0` so existing 5-arg callers/tests keep working.
- New getters (div-by-zero guarded where relevant):
  - `get_mapping_qualities()`, `get_aligned_lengths()`, `get_aligned_fractions()`
    (= `aligned/length` per read).
  - `get_total_aligned_base()`.
  - `get_mean_mapping_quality()`, `get_median_mapping_quality()`.
  - `get_mean_aligned_fraction()`, `get_median_aligned_fraction()`.
  - `get_read_count_with_aligned_fraction_below(threshold)` (default `0.9`).
  - `get_read_count_fully_aligned()` (`aligned_length == length`, integer equality).
  - `insertions_per_aligned_base()`, `deletions_per_aligned_base()`,
    `mismatches_per_aligned_base()` (total errors / aligned base; 0 when no aligned base).
- `__add__`: also sum `_total_aligned_base`.
- `__repr__`: add one line with mean aligned fraction and median MAPQ.

### 3. Summary tables (`report_table.py`)

**New** `create_mapping_table(readstat_list, readstat_sum)` → `table/mapping.txt`, columns:

```
label, read_count, query_base, aligned_base,
aligned_fraction_mean, aligned_fraction_median,
mapq_mean, mapq_median,
reads_aligned_fraction_lt_0.9, reads_fully_aligned
```

Per-contig rows plus a `'Total'` row (same shape as the other `create_*` functions).

**Rework** `create_readstat_table` (P0 #2) — final `read_stat.txt` column order:

```
label, read_count, total_base, aligned_base,
read_length_mean, read_length_median, read_length_N50, read_length_L50,
mean_insertion_per_read, mean_insertion_per_read_per_kb,
    insertion_per_query_kb, insertion_per_aligned_kb,
mean_deletion_per_read, mean_deletion_per_read_per_kb,
    deletion_per_query_kb, deletion_per_aligned_kb,
mean_mismatch_per_read, mean_mismatch_per_read_per_kb,
    mismatch_per_query_kb, mismatch_per_aligned_kb,
mean_intron_per_read, mean_intron_per_read_per_kb
```

where:
- `X_per_query_kb` **renames** the old `X_per_kb` (= total X / total query base × 1000) to
  disambiguate it from the other two.
- `X_per_aligned_kb` **is new** (= total X / aligned base × 1000).
- `mean_X_per_read_per_kb` and `mean_intron_per_read` are **unchanged** (internal contract,
  AGENTS.md #8).

**Rework** `create_splice_table` (P0 #1) → `table/splice.txt` becomes the four-category table
that already matches the HTML and figures:

```
label, gt-ag, gt-ag_pct, gc-ag, gc-ag_pct, at-ac, at-ac_pct, other, other_pct
```

and a **new** `create_splice_all_table` writes the full pair matrix (today's ~257-column
shape) to `table/splice_all.txt`, so no data is lost.

**Rework** `create_mismatch_normalized_read_location_table` (P0 #3): order the type columns in
the fixed canonical order `ac, ag, at, ca, cg, ct, ga, gc, gt, ta, tc, tg` (matching
`plot_mismatch_type_count`'s `miscolors`), and append a `bin_total` column (row sum).

### 4. Figures (`report_figure.py`)

New, each reusing `get_facet_row_col` / `determine_figure_size` and wired through
`generate_multiple_figs` (all / `Total` / per-contig):

- `plot_mapping_mapq_hist`: MAPQ histogram, integer bins 0–60 (extend to the observed max if
  a value > 60 appears), median annotated.
- `plot_mapping_aligned_fraction_hist`: aligned fraction histogram, 50 bins over [0, 1],
  median annotated.
- `plot_mapping_aligned_vs_query`: x = query length, y = aligned length, **hexbin density**
  (gridsize ~40) plus a y=x reference line (a raw scatter of 355k points is unreadable).

Bundled figure fixes:

- Introduce a shared categorical palette constant and use it in the two bar plotters that
  currently hard-code a 2-color list (`plot_readstat_bar`, `plot_element_total_count`), and in
  the new figures.
- `plot_indel_hist_length`: use explicit bins (1-bp steps for the short end, one catch-all tail
  bin) and/or a `log` y-axis so the insertion long tail (median 1, mean 28.5) is readable.

### 5. Wiring (`cli.py`)

- `o_files` additions:
  - `t_mapping` → `table/mapping.txt`
  - `f_mapping_mapq_hist` → `fig/mapping_hist_mapq`
  - `f_mapping_aligned_fraction_hist` → `fig/mapping_hist_aligned_fraction`
  - `f_mapping_aligned_vs_query` → `fig/mapping_scatter_aligned_vs_query`
  - `t_splice_all` → `table/splice_all.txt`
- Build `t_mapping = create_mapping_table(l_readstat, sreadstat)` and write it (same `to_csv`
  sep=`\t`, `index=False`); write `t_splice_all` alongside.
- Generate the three new figures via `generate_multiple_figs`.
- Extend the HTML assembly chain with `html_add_mapping(...)` (placed after the read-stat
  call), and add `'mapping'` to the tables dict passed to `html_add_data`.

### 6. HTML (`report_html.py` + `template.html`)

- New `html_add_mapping(html_string, mapping_table)`: renders the mapping summary table
  (`Total` row styled `table-secondary`) and substitutes the tokens
  `{%mapping_aligned_fraction_mean%}`, `{%mapping_aligned_fraction_median%}`,
  `{%mapping_mapq_mean%}`, `{%mapping_mapq_median%}`, `{%mapping_table%}`.
- Because `create_splice_table` now returns the count+percent columns directly, simplify
  `html_add_splice_table` to read those columns instead of recomputing counts/percentages from
  the full pair matrix (its computed `gt-ag/gc-ag/at-ac/other` values must equal the ones the
  table now carries). `html_add_mismatch_table` is unaffected by the type-column reorder and
  the added `bin_total` column (it looks types up by name and never reads `bin_total`).
- `template.html`: add a "Mapping" accordion section + nav entry (between "Read" and
  "Insertion") with cards for MAPQ histogram, aligned-fraction histogram, and the
  aligned-vs-query scatter (each `Total` + by-contig), plus the Statistic Summary table.
  Figure `src` references follow the existing `fig/<name>.png` / `fig/<name>.Total.png`
  convention so `inline_figures` handles them unchanged.

### 7. Exports (`__init__.py`)

Export `create_mapping_table`, `create_splice_all_table`, the three `plot_mapping_*` functions,
and `html_add_mapping`; keep `__all__` in sync.

### 8. Conventions foundation

- Add a small formatting helper (single source of truth for table float rounding) and the
  shared color palette constant. The new mapping table/figures use them; the bundled table and
  figure reworks adopt them too (this is the "round consistently across txt and HTML" fix,
  P0 #4).

### 9. Column-name contract (AGENTS.md #8)

- `mean_mismatch_per_read_per_kb`, `mean_insertion_per_read_per_kb`,
  `mean_deletion_per_read_per_kb`, `mean_intron_per_read` stay intact.
- Three **documented breaking changes** for external TSV/JSON consumers: `insertion_per_kb`,
  `deletion_per_kb`, `mismatch_per_kb` → `*_per_query_kb`; `splice.txt` columns become the
  four-category set (full matrix moves to `splice_all.txt`). None of these columns are read
  internally by `report_html.py`/`cli.py`. Recorded here and to be noted in the release notes /
  `docs/reporting-and-output.md`.

### 10. Documentation

- Update `docs/reporting-and-output.md`: add `mapping.txt` + `splice_all.txt` to the output
  layout and the summary-table column lists; rename the `X_per_kb` columns in the read_stat
  doc to `*_per_query_kb` and document `*_per_aligned_kb`; document the new `mapping_*`
  figures and `mapping_*` HTML tokens.

## User Experience

- `table/mapping.txt` and the "Mapping" HTML section surface, per contig and total: aligned
  vs query base, aligned-fraction mean/median, MAPQ mean/median, and counts of low-alignment
  and fully-aligned reads.
- Known edge cases:
  - Reads with `query_length == 0` (should not occur for mapped records) are guarded in every
    division.
  - A BAM entirely without secondary/supplementary (the chr22 example) and one with them behave
    the same: both metrics are per-alignment-record, consistent with the current read-counting
    semantics (no flag filtering is introduced in this design).
  - `aligned_base == 0` → per-aligned-kb columns are `0` (never NaN/Inf).

## Implementation Notes / Trade-offs

- **Pro:** Approach A reuses `StatBlock`, `sum`, `relabel`, and `get_stat_list` untouched;
  memory grows ~2 ints/read (~few MB at 355k reads). Tag-agnostic (works for cs and MD).
- **Trade-off:** `ReadStat`'s semantics broaden slightly from "length + error counts" to
  "length + error counts + mapping"; accepted to avoid threading a sixth stat object through
  the whole pipeline.
- **Trade-off:** `X_per_kb → X_per_query_kb` and the `splice.txt` reshaped columns are
  breaking for external consumers, but pre-1.0 and confined to columns not read internally;
  `splice_all.txt` preserves the full splice data.
- **Explicitly NOT being done:** MAPQ re-normalization; GC content (#16); intron-length
  distribution (#14); homopolymer analysis (#15); the top-level dashboard (P2 #12); any
  secondary/supplementary flag filtering.

## Success Criteria

1. `uv run ruff check` and `uv run pytest` pass, including new unit tests for the new
   `ReadStat` getters, `create_mapping_table`, and `create_splice_all_table`; existing tests
   that call `add_read`/`ReadRecord` still pass (defaults/compat checked).
2. Re-running `tmp/analyze_chr22.sh` produces `table/mapping.txt` whose Total row matches the
   probe within rounding: `aligned_fraction_median ≈ 0.840`, `mapq_median = 1`,
   `reads_aligned_fraction_lt_0.9 ≈ 194381`; `read_stat.txt` contains `aligned_base`,
   `*_per_query_kb`, and `*_per_aligned_kb`; `splice.txt` has 9 columns (4 categories + pct)
   and `splice_all.txt` has the full pair matrix.
3. `LQC_report.html` remains self-contained and offline-safe, has no unreplaced `{%...%}`
   tokens, and renders the new "Mapping" section (including base64-inlined figures).
4. `fig/` gains `mapping_hist_mapq{,.Total,.<contig>}`, `mapping_hist_aligned_fraction{...}`,
   `mapping_scatter_aligned_vs_query{...}` PNG+PDF.

## Dependencies

None beyond the existing stack (pysam, pandas, matplotlib, numpy). This builds on the current
`ReadStat`/`StatBlock` pipeline and the HTML token-substitution finalize pass.

## Decisions (resolved, no placeholders)

- Storage: **A** (extend `ReadStat`), not B (new stat object).
- Scope: **all-in** — core #13 + bundled P0/P1 polish + P0 #2 aligned-base normalization.
- Aligned fraction definition: `query_alignment_length / query_length` (per read), `int`-exact
  "fully aligned" when equal.
- Normalization trio kept: `mean_X_per_read_per_kb` (contract), `X_per_query_kb` (renamed),
  `X_per_aligned_kb` (new). Intron columns untouched this round.
- Splice: `splice.txt` → 4 categories; full matrix → `splice_all.txt` (always written).
- Mapping-aligned-vs-query plot: hexbin, not scatter.