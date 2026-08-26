# Changelog

All notable changes to LQC are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and LQC adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.8] - 2026-08-27

### Fixed
- The indel-length histogram (`plot_indel_hist_length`) now bins one integer width per length (1 through the observed max) instead of a short-end 1-bp + catch-all tail, so long tails no longer render as one flat plateau bar.
- Guard against zero-length reads in `ReadStat` (divide-by-zero warning) and against zero-read-length `cs` tags in the normalized-location pass.
- Make worker `stat_region` teardown exception-safe and pin multi-contig task ordering.

### Performance
- Worker-side fetch (M0): each worker opens its own read-only BAM handle and processes a balanced per-contig coordinate window, replacing the global read materialization + `Pool.map` pickle pass; whole-genome inputs now complete within a few GB instead of exhausting RAM.
- Single-pass `cs` scan (M1): a manual tokenizer replaces the regex-based `cs_to_list`, and count + normalized indel/mismatch location computation are fused into one pass and cached.
- Numpy columnar storage (M2): `ReadStat`/`Indel`/`Mismatch` accumulate into numpy arrays and merge via `np.concatenate` (no deepcopy); `Splice.__add__` drops its deepcopy.
- Micro-fix (M3): `convert_reverse_complement` uses a single translate table.
- Outputs (`table/*.txt`, `fig/*`, HTML, `read.cs`) stay byte-identical to 0.0.7.

## [0.0.7] - 2026-08-26

### Added
- Mapping metrics: per-read mapping quality and aligned length are now captured, a mapping summary table is emitted, and a dedicated "Mapping" section (MAPQ, aligned-fraction, and aligned-vs-query histograms) is embedded in the HTML report.
- A `splice_all` table alongside the four-category splice table, including an explicit "other" bucket.
- Aligned-base normalization rates in the read-stat table.

### Changed
- TSV float columns are rounded to 4 significant digits, matching the HTML report's formatting.
- The indel-length histogram uses a log scale and bar plots share a common color palette; the mismatch-type column order now ends with a per-bin total.
- The splice table is collapsed to four categories.

### Fixed
- MAPQ bins are capped at 60 and the mapping histograms guard against empty reads.
- Splice-table per-contig fallback when references are omitted.

## [0.0.6] - 2026-08-23

### Added
- Self-contained offline HTML report: Bootstrap CSS/JS inlined, summary tables embedded as JSON, and figures embedded as base64 data URIs.
- Comprehensive pytest suite: unit tests for `ReadStat`, `Indel`, `Mismatch`, `Splice`, `utils`, and the report modules, plus integration tests over a tiny indexed BAM and the real chr22 BAM.

### Changed
- `-t/--thread` now defaults to `min(cpu_count, 4)` (was `1`).
- `read.cs` (from `--output-cs`) now contains only the analyzed contigs in requested order instead of the whole BAM; the CLI logs a warning when references are omitted.

### Fixed
- Recorded deletions from the cs `del` operand instead of the insertion operand.
- `--output-cs` now writes `read.cs` atomically and removes per-chunk temp files even on failure.

### Performance
- `-t/--thread` now parallelizes over read chunks (previously one task per contig), so single- and multi-contig BAMs both benefit; table and `read.cs` outputs are byte-identical across thread counts.
- `--output-cs` is emitted in a single pass by the stat workers instead of a second full BAM parse.
- cs parser micro-optimizations: removed intron donor/acceptor regexes, grouped insertion/deletion/mismatch extraction into one filtered pass, and translated complements via a lookup table.

### Internal / Packaging
- Moved the package to a `src/lqc` layout.
- Replaced `setup.py` with `pyproject.toml` (hatchling) + `uv`, added a `lqc` console entry point, and exposed `__version__`.
- Adopted `ruff` as the lint gate (frozen ruleset); migrated the test suite from `unittest` to pytest.