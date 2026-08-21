# Architecture & Data Flow

> **When to read:** you need the big picture, are adding a feature, or are tracing a bug
> that crosses modules. **Expire:** replace when the pipeline or module split changes.

## Pipeline (who calls whom)

`lqc/cli.py` orchestrates the whole run in one linear sequence; the package is a library of
stat classes and report renderers that `lqc/cli.py` composes.

1. Parse CLI args, set up logging (`lqc/cli.py`).
2. `check_bam_with_cs_or_md()` decides the parse method (`cs`, `MD`, or `both`→`cs`).
3. `mp.Pool` maps `stat_element_from_bam_by_contig()` over the requested contigs.
4. Each worker returns a 5-tuple `(ReadStat, Indel, Indel, Mismatch, Splice)`.
5. `sum()` combines per-contig objects into one `'Total'` object per type.
6. `report_table` turns objects into pandas DataFrames; `report_figure` renders PNG/PDF;
   `report_html` injects tables into the HTML template.

Each of the four stat classes implements `__add__`/`__radd__`, so Python's `sum()` (which
seeds with `0`) produces a merged object. That is the only reason `lqc/cli.py` can do
`sum(l_readstat)`.

> **Why** this layering exists: the element classes hold no I/O, `stat.py` is the only
> module that opens BAM/FASTA, and `report_*` modules only format/plot already-computed
> objects. **Expire:** if the boundary is broken (e.g. a report module opening a BAM).

## Module map

| Module | Role | Key symbols |
| --- | --- | --- |
| `lqc/cs.py` | Parse a `cs` tag (or CIGAR+MD+seqs) into an element list; `CS` object. Single source of truth for tag parsing. | `cs_to_list`, `cigar_to_list`, `convert_cigar_md_to_cs_list`, `CS` |
| `lqc/readstat.py` | Per-contig read metrics (count, lengths, N50/L50, per-read element means). | `ReadStat` |
| `lqc/indel.py` | Insertion/deletion counts, length distribution, normalized location. | `Indel` |
| `lqc/mismatch.py` | Mismatch type counts and binned normalized location. | `Mismatch` |
| `lqc/splice.py` | Splice (intron) pair/site counts. | `Splice` |
| `lqc/stat.py` | Open BAM/FASTA, iterate one contig, build the five objects above. The multiprocessing worker. | `stat_element_from_bam_by_contig` |
| `lqc/utils.py` | Small helpers: file-type detection, tag detection, complement conversion, `write_readcs`. | `bam_or_sam`, `check_bam_with_cs_or_md`, `write_readcs` |
| `lqc/report_table.py` | Build pandas DataFrames for the `.txt` summary tables. | `create_*_table` |
| `lqc/report_figure.py` | matplotlib (Agg) figures. | `plot_*` |
| `lqc/report_html.py` | Substitute `{%token%}` placeholders in the template with computed tables. | `html_add_*` |
| `lqc/template/` | `template.html` + SVG logos shipped as package data. | `get_html_template`, `copy_logo` |
| `tests/` | pytest unit + integration tests (outside the package). | `test_cs.py`, `test_cli.py`, `test_stat.py`, ... |
| `tests/data/` | Committed small CS/CIGAR+MD fixture used by the unit tests. | `cs_test.test_data` |
| `scripts/` | Developer utilities (not shipped, not tested). | `generate_cs_test_data.py` |
| `tmp/data/` | Large real BAM + index for manual smoke runs (gitignored, not CI). | `ENCFF417VHJ.chr22.sorted.bam` |
| `lqc/cli.py` | CLI entry point; composes the above and defines every output filename. | — |

## The data objects

- `CS` — a parsed alignment, holding the element list and precomputed counts. Built via
  `from_cs_tag_string` or `from_cigar_string`; used by `stat.py` to feed the other four
  classes.
- `ReadStat`, `Indel`, `Mismatch`, `Splice` — per-contig accumulators with getters used
  directly by `report_table.py` and `report_figure.py`. Their public getters (e.g.
  `get_mean_length()`, `get_length_NL(50)`, `get_type_count()`) are the interface the
  report layer depends on.

## Two rules that keep the layering safe

- **Opening input files belongs only to `lqc/stat.py` and `lqc/utils.py`.**
  *why: keeps report modules pure and testable · when: adding new input handling ·
  expire: if an explicit I/O layer is introduced.*
- **A new output artifact is wired end-to-end in `lqc/cli.py` only** (compute → table/figure
  → `report_html` → filename in `o_files`). *why: `lqc/cli.py` is the single composition
  root · when: adding a table, figure, or report section · expire: if composition moves
  out of `lqc/cli.py`.*

## What is documented in code, not here

Per-method semantics (formulas, defaults, arguments) live in the docstrings and getters of
each class. Do not re-document them here — read the module when you touch it.