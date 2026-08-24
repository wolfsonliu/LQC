# Reporting & Output Format

> **When to read:** before touching `src/lqc/report_table.py`, `src/lqc/report_figure.py`,
> `src/lqc/report_html.py`, or `src/lqc/template/`. Full code (formats, columns, colors) lives in
> the modules; this doc records the cross-file contracts.

## Output directory layout

`src/lqc/cli.py` builds this structure under `-o OUT`:

```
OUT/
  LQC_report.html
  result.pickle            # only with --output-pickle
  read.cs                  # only with --output-cs
  table/
    read_stat.txt
    insertion.txt
    deletion.txt
    mismatch.txt
    splice.txt
    mapping.txt
    splice_all.txt
  fig/                     # many PNG + PDF pairs (incl. mapping_* figures), names fixed in o_files
```

All tables are **tab-separated** with `index=False`. *why: users/later analysis read them
· when: changing table writers · expire: if a versioned format is introduced.*

## Summary tables (columns)

- `read_stat.txt` — `create_readstat_table`: `label, read_count, total_base, aligned_base,
  read_length_mean/median/N50/L50, mean_{insertion,deletion,mismatch,intron}_per_read,
  mean_*_per_read_per_kb, {insertion,deletion,mismatch}_per_query_kb,
  {insertion,deletion,mismatch}_per_aligned_kb`. Per-contig rows plus a `'Total'` row.
  (`total_base` is the query-length denominator of the `*_per_query_kb` rates.)
- `insertion.txt` / `deletion.txt` — `create_indel_summary_table`: `label, total_count,
  total_length, mean_length, median_length`.
- `mismatch.txt` — `create_mismatch_normalized_read_location_table`: `label, bin` plus one
  column per mismatch type (in a fixed canonical order), binned into 10 equal read-location
  bins (`cuts = 0..1`), plus a per-bin `bin_total` column.
- `splice.txt` — `create_splice_table`: `label, gt-ag, gt-ag_pct, gc-ag, gc-ag_pct,
  at-ac, at-ac_pct, other, other_pct` (four canonical categories). The full per-pair matrix
  is preserved in `splice_all.txt`.
- `mapping.txt` — `create_mapping_table`: `label, read_count, query_base, aligned_base,
  aligned_fraction_mean/median, mapq_mean/median, reads_aligned_fraction_lt_0.9,
  reads_fully_aligned`.
- `splice_all.txt` — `create_splice_all_table`: `label` plus one column per observed splice
  pair (full matrix; column set derives from the aggregate).

**Contract:** column names are looked up by name downstream (e.g. `report_html.py` reads
`mean_mismatch_per_read_per_kb`, `mean_intron_per_read`). Renaming a column breaks the HTML
report. *why: no shared schema object exists yet · when: changing report_table.py ·
expire: if a shared schema/constants module is introduced.*

## Figures

`src/lqc/report_figure.py` renders each stat as: per-contig facets (one figure), a `'Total'`
figure, and one figure per contig (`.png` and `.pdf` for every file, via `savefig` in
`src/lqc/cli.py`). Mapping-specific figures are `plot_mapping_mapq_hist` (MAPQ history),
`plot_mapping_aligned_fraction_hist` (aligned-fraction histogram), and
`plot_mapping_aligned_vs_query` (aligned-vs-query-length hexbin scatter). Feature names in
`plot_readstat_bar` must exactly match the strings listed in `src/lqc/cli.py` (`"Read count"`,
`"Median read length"`, …). All figure code assumes the `Agg`
backend. *why: filenames and facet titles are derived from these strings · when: adding a
figure or feature · expire: if plotting is moved behind a config-driven registry.*

## HTML report contract

- `src/lqc/template/template.html` contains `{%token%}` placeholders; `src/lqc/report_html.py`
  (`html_add_*`) is the **only** code that substitutes them. Token groups: `readstat_*`,
  `mapping_*`, `mismatch_*`, `insertion_*`, `deletion_*`, `intron_*`, `splice_*`, plus the self-contained
  asset tokens `bootstrap_css`, `bootstrap_js`, `lqc_data`, `version`.
- Every placeholder in the template must be replaced for a given run; the `html_add_*`
  functions are chained in `src/lqc/cli.py` and must cover the template's full token set. Add a
  new token only together with the `html_add_*` that fills it.
  *why: unreplaced tokens would surface verbatim as `{%...%}` in the report · when: adding
  a report section · expire: if template rendering moves to a real templating engine.*
- The `'Total'` row is styled specially (`table-secondary`) and also sourced for the
  headline metric values — keep the aggregate row labeled exactly `'Total'`.
- **Self-contained output:** `src/lqc/cli.py` finalizes the report with
  `html_add_bootstrap` (inlines the vendored Bootstrap CSS/JS from `src/lqc/template/`,
  MIT-licensed), `html_add_data` (embeds the summary tables — readstat, mapping,
  insertion, deletion, mismatch, splice — as JSON in
  `<script type="application/json" id="lqc-data">`), and `inline_figures` (rewrites
  `src="fig/..."` to base64 `data:` URIs). The rendered `LQC_report.html` therefore
  references no `fig/` files and no remote resources and opens offline as a single file.

## Things already documented in code (do not duplicate)

- Exact figure layout/colors (`magentas` palette, `get_facet_row_col` grid) — in
  `src/lqc/report_figure.py`.
- Exact per-placeholder replacement strings and table markup — in `src/lqc/report_html.py` and
  `src/lqc/template/template.html`.