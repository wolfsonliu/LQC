# LQC HTML Report: Self-Contained, Offline-Safe Single File (Design)

Date: 2026-08-22 · Branch: `improvement` · Status: approved (Approach A, incl. embedded JSON)

## Problem Statement

`LQC_report.html` is not a standalone artifact today. Opening the single file by itself
shows a broken page because it depends on runtime-remote resources:

- Bootstrap CSS and the JS bundle are loaded from `cdn.jsdelivr.net`.
- Every figure is a static PNG/SVG referenced by `<img src="fig/...">`, so the HTML must
  live next to the `fig/` directory to render.
- Two `img.shields.io` remote badge images are loaded over the network.

Summary tables and headline metrics are already embedded in the file (Python renders them
into `<table>` markup and `{%token%}` values), so only styling, scripting, and figures are
missing from the file. The goal is a **single, self-contained HTML file** that opens
correctly offline and can be shared by copying one file — without changing the existing
PNG/PDF/TSV artifacts.

## Proposed Solution

Do not rewrite matplotlib or re-render charts in JavaScript. Inline every runtime resource
the report needs directly into `LQC_report.html`: the vendored Bootstrap CSS/JS, the
figures the page actually references (as base64 data URIs), and a machine-readable JSON
copy of the summary tables. Replace the remote badge images with local links. Leave the
`fig/` (PNG+PDF) and `table/` (TSV) outputs exactly as they are.

This is "Approach A" from the brainstorming session: it achieves a 100% self-contained,
offline-safe file with near-zero rewriting and no second (JS) copy of the plotting logic to
keep in sync.

## Detailed Design

### 1. Vendor Bootstrap 5.1.3 (remove CDN dependency)

- Add two vendor files under `lqc/template/`:
  - `lqc/template/bootstrap.min.css` (Bootstrap 5.1.3, MIT license; keep the license/version
    header comment intact).
  - `lqc/template/bootstrap.bundle.min.js` (same version/license; includes Popper).
- MIT is GPL-compatible, satisfying the license constraint (AGENTS.md #10).
- In `lqc/template/template.html`, replace the two CDN lines:

  - `<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css"
    ...>` → `{%bootstrap_css%}`
  - `<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"
    ...></script>` → `{%bootstrap_js%}`

- A new `report_html.py` function substitutes those tokens by wrapping the vendored file
  contents in `<style>...</style>` and `<script>...</script>` respectively. Visual output
  stays pixel-identical to today.

### 2. Inline referenced figures as base64 data URIs

- Add a render-time post-process step `inline_figures(html_string, fig_dir)` in
  `lqc/report_html.py`. It:
  - Scans for every `src="fig/<name>.png"` / `src="fig/<name>.svg"` in the HTML.
  - Reads the file from `fig_dir`, base64-encodes it, and rewrites the `src` to
    `data:image/png;base64,<...>` or `data:image/svg+xml;base64,<...>` (MIME chosen by
    extension; no PDFs are referenced by the HTML).
  - Only inlines the ~three dozen images the page actually uses; it does **not** touch the
    remaining per-contig PNG/PDF files under `fig/`.
- `cli.py` calls it as the **last** step of the HTML assembly chain, after `copy_logo` and
  after all `html_add_*` calls, because `copy_logo` materializes `fig/logo_white.svg` and
  `fig/card.svg` that the inliner must then embed.
- Expected size increase is about the PNG set the page shows (on the order of ~1.5 MB in
  the sample run); acceptable for the offline/single-file goal.

### 3. Replace remote badge images with local links

In `template.html`, remove the two remote `<img src="https://img.shields.io/...">` badges
(version and GitHub stars) and replace them with plain local markup that carries the same
information without network access:

- Version → a locally-styled `LQC v0.0.5` text (source the version string from
  `{%version%}` so it does not drift from `lqc.__version__`).
- Stars → a local `github.com/gxiaolab/LQC` text link.

Plain `<a href="https://...">` anchors (PyPI, GitHub) are kept; they are navigation
targets, not resource loads, and do not affect offline rendering.

### 4. Embed machine-readable JSON of the summary tables

- Add one `<script type="application/json" id="lqc-data">{...}</script>` block to the
  rendered output, serializing the five summary DataFrames — `read_stat`, `insertion`,
  `deletion`, `mismatch`, `splice` — as `records`-oriented JSON under top-level keys
  `readstat`, `insertion`, `deletion`, `mismatch`, `splice`.
- Escape it for safe inline placement in HTML: replace any `</` sequence (in particular
  `</script>`) so the injected JSON cannot terminate its own script tag.
- Volume is a few KB. This makes the HTML self-describing ("the key data lives in the file")
  and parseable by downstream scripts.

### 5. Wiring and code locations

- `lqc/report_html.py`: add `html_add_bootstrap(html_string)` (fills
  `{%bootstrap_css%}`, `{%bootstrap_js%}`, and `{%version%}` from `lqc.__version__`),
  `html_add_data(html_string, tables)`, and `inline_figures(html_string, fig_dir)`.
- `lqc/template/template.html`: token changes in §1 and §3.
- `lqc/template/bootstrap.min.css`, `lqc/template/bootstrap.bundle.min.js`: new vendor files.
- `lqc/cli.py`: in the report block, extend the existing `html_add_*` chain with the new
  calls in order — `html_add_bootstrap` → `html_add_data` → `inline_figures` — before
  writing `o_files['html']`.
- `lqc/__init__.py`: export the new public functions (and keep `__all__` in sync), following
  the existing export pattern used for the current `html_add_*` functions.
- `docs/reporting-and-output.md`: document the new tokens, the vendor files, the finalize
  pass, and the new contract "the rendered HTML references no external or `fig/` resources."

### 6. Symbol contract preserved

- The `{%token%}` substitution contract is maintained: `html_add_*` functions remain the
  only code that substitutes tokens, and any new token is introduced together with the
  function that fills it (AGENTS.md invariants in the reporting doc).
- Figure rendering (`report_figure.py`), table building (`report_table.py`), the column
  names, and the `'Total'` sentinel are untouched.

## User Experience

- **Offline / single file:** copying `LQC_report.html` alone, with the network disabled, opens
  a fully styled page with all five accordion sections, all tables, and all figures.
- **No behavioral change otherwise:** navigation, accordion toggling, styling, and figure
  content look the same as today. The only visible difference is that the two remote
  `img.shields.io` badges are replaced by local text/links (no network fetch). Charts remain
  static (no hover/zoom — interactivity is out of scope).
- **Data access:** tables remain human-readable in-page and now also machine-readable via the
  embedded `#lqc-data` JSON.
- **Edge cases:**
  - A figure file referenced by the template but missing from `fig/` should raise a clear
    error during `inline_figures` (fail loudly rather than emit a broken `data:` URI). This
    mirrors the "every placeholder must be replaced" safety posture.
  - Non-image referenced assets (none expected today) are left untouched.

## Implementation Notes / Trade-offs

- **Pro:** no second rendering path. The PNG/PDF and the embedded PNG share one matplotlib
  source of truth, so they can never drift.
- **Trade-off:** vendoring Bootstrap adds ~270 KB of third-party (MIT, GPL-compatible) files
  to the repository; base64 figures increase the HTML by roughly the size of the shown PNGs.
  This is the accepted cost for a true single-file, offline report.
- **Explicitly NOT being done:** no Plotly/ECharts/Chart.js; no JS re-rendering of charts; no
  removal of `fig/` or `table/` outputs; no change to plotting/table code.

## Success Criteria

1. `uv run lqc -b <cs-tagged.bam> -o out` produces `out/LQC_report.html` containing:
   - no `src="fig/` and no `src="https://` (nor any other remote resource URL) — `<a
     href="https://...">` anchors are allowed;
   - at least one `data:image/...;base64`; and inline `<style>`/`<script>` from the vendored
     Bootstrap.
2. Opening the report in a browser with the network disabled renders all five sections,
   tables, and figures with the same appearance as before.
3. `fig/` and `table/` artifacts are byte-identical to the previous version's outputs.
4. `uvx ruff check` passes; unit tests for `inline_figures` / JSON escaping / token coverage
   and a smoke assertion "no external or fig/ resource references" all pass.

## Dependencies

- Bootstrap 5.1.3 distribution files (MIT) to vendor. Nothing else new; this builds directly
  on the existing `report_html.py` string-substitution pipeline.

## Decisions (resolved, no placeholders)

- Approach: **A** (self-contained via inlining), chosen over B (full JS re-render) and C
  (hybrid JS + PNG).
- Keep all concurrent artifacts: **PNG, PDF, TSV stay**.
- Vendor Bootstrap into the repo and inline it: **yes**.
- Embed the summary-table JSON: **yes, in this change (v1)**.