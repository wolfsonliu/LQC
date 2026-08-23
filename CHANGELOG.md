# Changelog

All notable changes to LQC are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and LQC adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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