# Input & Tag Conventions

> **When to read:** before touching `src/lqc/cs.py`, `src/lqc/stat.py`, or `src/lqc/utils.py`, or when
> diagnosing a parsing bug. Full per-method formulas stay in code; this doc records the
> conventions that are easy to get wrong.

## Input requirements

- The alignment file must end in `.bam` or `.sam` (case-insensitive check in
  `bam_or_sam`).
- BAM must be **sorted and indexed** — `pysam.AlignmentFile.fetch(contig)` requires the
  index. *why: random access per contig · expire: if an index-free scan replaces `fetch`.*
- The file must carry `cs` tags or `MD` tags. `check_bam_with_cs_or_md` scans the first
  reads and returns `None` if neither is present.
  - `cs` only → parse via `cs`, no reference FASTA needed.
  - `MD` only → parse via `MD`+`CIGAR`, **requires `--genome-fasta`** (splice reconstruction
    needs the reference sequence).
  - both → `cs` wins (case downgraded to `'cs'`).
  *why: `cs` encodes splice donor/acceptor, so genome can be skipped · expire: if tag
  support changes.*

## `cs` tag format (minimap2 `--cs=long`)

The tag is a run-length string; `cs_to_list` expands it into `[low, high, mark, value]`
items where positions are **0-based half-open on the reference**.

| Mark | Meaning | Length consumed (reference) | Example |
| --- | --- | --- | --- |
| `:` | match | `int(value)` | `:29` |
| `*` | substitution (1 bp) | 1 | `*ga` |
| `+` | insertion (relative to reference) | 0 | `+cccc` |
| `-` | deletion (relative to reference) | `len(value)` | `-at` |
| `~` | intron/splice | numeric prefix | `~ct140ac` → 140 bp, `ct` donor, `ac` acceptor |

Read length is reconstructed as the sum of `:` + `*` + `+` contributions (deletions and
introns consume reference but not read). See `_compute_counts` in `src/lqc/cs.py`.

> **Why** this table is a hard constraint: every downstream count (insertions, deletions,
> mismatches, splice pairs) derives from these five tokens. **When:** any change to `cs.py`
> parsing. **Expire:** only if a different cs dialect is adopted.

## `MD` + `CIGAR` fallback

When only `MD` is present, `convert_cigar_md_to_cs_list` rebuilds an equivalent cs list
from the CIGAR string, MD string, query sequence, and reference sequence. It is expected
to produce output **identical** to the cs path — asserted in `tests/test_cs.py`.

## Strand & complement handling

- Reverse-strand reads are corrected by `convert_reverse_complement` (and
  `convert_complement`), defined in `src/lqc/utils.py`.
- The `CS` object stores a `strand` (`'+'`/`'-'`) and uses `_start_pos`/`_contig` to
  translate relative → contig coordinates (`get_contig_position`).

## Coordinate conventions (easy to get wrong)

- Parsed/tag coordinates are **0-based half-open** (`[low, high)`).
- `get_relative_position()` returns tag-relative intervals; `get_contig_position()` adds
  `_start_pos`. Any 1-based/closed presentation must happen at the report boundary only.
  *why: mixing bases off-by-one is the top source of subtle bugs · expire: only if the
  whole pipeline migrates in one coordinated change.*

## `read.cs` output (`--output-cs`)

`--output-cs` writes one line per cs element into `read.cs` with six columns:
`read name, contig, low, high, cs mark, cs value`. Keep this schema if you extend the
option. *why: users may consume `read.cs` directly · expire: if the file format is
deliberately versioned/renamed.*

`read.cs` now contains only the reads of the analyzed contigs (`-c`, or the default
chromosome list), in requested-contig order — not every reference in the BAM. The CLI
logs a warning when it omits any header reference. `src/lqc/utils.py` `write_readcs`
still dumps the entire BAM in file order for library callers.

## How to verify a parsing change

`tests/test_cs.py` compares the cs path against the CIGAR+MD path over the committed
`tests/data/cs_test.test_data` fixture, and checks element counts on one long read. Extend
that fixture (`scripts/generate_cs_test_data.py`) rather than hand-rolling asserts.