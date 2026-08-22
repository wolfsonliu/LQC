# Self-Contained Offline LQC HTML Report — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `LQC_report.html` a single, self-contained file that renders fully offline by inlining vendored Bootstrap CSS/JS, embedding referenced figures as base64 data URIs, replacing remote badge images with local links, and embedding the summary tables as JSON — while leaving `fig/` (PNG+PDF) and `table/` (TSV) outputs byte-identical.

**Architecture:** No matplotlib or table changes. The change is confined to the HTML assembly pipeline: `template.html` swaps CDN `<link>`/`<script>` and remote badges for `{%...%}` tokens; `lqc/report_html.py` gains three functions (`html_add_bootstrap`, `html_add_data`, `inline_figures`) that fill those tokens and post-process `<img src="fig/...">` into base64; `lqc/cli.py` chains the new functions at the end of the existing `html_add_*` sequence.

**Tech Stack:** Python 3.11, Bootstrap 5.1.3 (vendored, MIT), stdlib `base64`/`json`/`os`/`re`, pandas DataFrame → `to_json(orient='records')`, pytest (dev group), ruff (via `uvx`).

---

## File Structure

| Action | Path | Responsibility |
| --- | --- | --- |
| Create | `lqc/template/bootstrap.min.css` | Vendored Bootstrap 5.1.3 CSS (MIT) |
| Create | `lqc/template/bootstrap.bundle.min.js` | Vendored Bootstrap 5.1.3 JS bundle incl. Popper (MIT) |
| Modify | `lqc/template/template.html` | Replace CDN link/script + remote badges with `{%bootstrap_css%}` / `{%bootstrap_js%}` / `{%lqc_data%}` / `{%version%}` tokens and local badges |
| Modify | `lqc/report_html.py` | Add `html_add_bootstrap`, `html_add_data`, `inline_figures` |
| Modify | `lqc/__init__.py` | Export the three new functions |
| Modify | `lqc/cli.py` | Import the new functions and chain them before writing HTML |
| Create | `tests/test_report_html_assets.py` | Unit tests for template tokens + new functions |
| Modify | `tests/test_cli.py` | Add self-containment end-to-end assertion |
| Modify | `docs/reporting-and-output.md` | Document new tokens, vendor files, and the self-contained contract |

Convention used in this plan: new unit tests import the function under test **locally inside the test function** (e.g. `from lqc.report_html import html_add_bootstrap`). This keeps each task's failing test runnable on its own without import errors from not-yet-written siblings.

---

## Task 1: Vendor Bootstrap 5.1.3 (MIT)

**Files:**
- Create: `lqc/template/bootstrap.min.css`
- Create: `lqc/template/bootstrap.bundle.min.js`

- [ ] **Step 1: Download the two Bootstrap 5.1.3 distribution files**

```bash
curl -sSLo lqc/template/bootstrap.min.css \
  https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css
curl -sSLo lqc/template/bootstrap.bundle.min.js \
  https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js
```

- [ ] **Step 2: Verify the files carry the expected version + MIT license header**

Run:

```bash
grep -o "Bootstrap v5.1.3" lqc/template/bootstrap.min.css | head -1
grep -o "Bootstrap v5.1.3" lqc/template/bootstrap.bundle.min.js | head -1
head -c 300 lqc/template/bootstrap.min.css
```

Expected: both `grep` commands print `Bootstrap v5.1.3`; the `head` shows `/*! * Bootstrap v5.1.3 ... Licensed under MIT (...) */`. (MIT is GPL-compatible — satisfies AGENTS.md constraint #10.)

- [ ] **Step 3: Commit**

```bash
git add lqc/template/bootstrap.min.css lqc/template/bootstrap.bundle.min.js
git commit -m "vendor: add Bootstrap 5.1.3 CSS and JS bundle (MIT)"
```

---

## Task 2: Add asset tokens + replace remote badges in the template

**Files:**
- Modify: `lqc/template/template.html`
- Create: `tests/test_report_html_assets.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_report_html_assets.py`:

```python
from lqc import get_html_template


def test_template_has_asset_tokens_and_no_remote_resources():
    html = get_html_template()
    assert '{%bootstrap_css%}' in html
    assert '{%bootstrap_js%}' in html
    assert '{%lqc_data%}' in html
    assert '{%version%}' in html
    assert 'cdn.jsdelivr.net' not in html
    assert 'img.shields.io' not in html
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
uv run pytest tests/test_report_html_assets.py -v
```

Expected: FAIL — the `{%bootstrap_css%}` / `{%bootstrap_js%}` / `{%lqc_data%}` / `{%version%}` assertions (and the `cdn.jsdelivr.net` one) fail because the template still has CDN references and no asset tokens.

- [ ] **Step 3: Replace the CDN link/script header with tokens**

In `lqc/template/template.html`, replace:

```html
<!-- Bootstrap CSS -->
<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3" crossorigin="anonymous">
<!-- JavaScript Bundle with Popper -->
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
<title>LQC report</title>
```

with:

```html
<!-- Bootstrap CSS (inlined at render time) -->
{%bootstrap_css%}
<!-- Bootstrap JavaScript bundle with Popper (inlined at render time) -->
{%bootstrap_js%}
<!-- Machine-readable summary tables -->
{%lqc_data%}
<title>LQC report</title>
```

- [ ] **Step 4: Replace the two remote badge images with local badges**

In the About section, replace:

```html
                <a href="https://pypi.org/project/lqc/">
                  <img src="https://img.shields.io/badge/version-v0.0.5-7a0177"/>
                </a>
                <a href="https://github.com/gxiaolab/LQC">
                  <img src="https://img.shields.io/github/stars/gxiaolab/LQC?style=social"/>
                </a>
```

with:

```html
                <a href="https://pypi.org/project/lqc/">
                  <span class="badge rounded-pill text-light" style="background-color:#7a0177;">LQC v{%version%}</span>
                </a>
                <a href="https://github.com/gxiaolab/LQC">
                  <span class="badge rounded-pill text-dark" style="background-color:#e9ecef;">github.com/gxiaolab/LQC</span>
                </a>
```

- [ ] **Step 5: Remove the third remote badge and update the now-stale guidance text**

Replace:

```html
                <li>
                  The correct view of this page requires network
                  connection to bootstrap
                  (<a href="https://getbootstrap.com/">
                    <img src="https://img.shields.io/badge/bootstrap-v5.1.3-blueviolet"/>
                  </a>).
                </li>
```

with:

```html
                <li>
                  This page is self-contained: all styling, scripts, and
                  figures are embedded, so it can be opened offline.
                </li>
```

- [ ] **Step 6: Run the test to verify it passes**

```bash
uv run pytest tests/test_report_html_assets.py -v
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add lqc/template/template.html tests/test_report_html_assets.py
git commit -m "feat: add self-contained asset tokens to HTML template"
```

---

## Task 3: Implement `html_add_bootstrap`

**Files:**
- Modify: `lqc/report_html.py`
- Modify: `tests/test_report_html_assets.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_report_html_assets.py`:

```python
def test_html_add_bootstrap_inlines_vendored_assets():
    from lqc.report_html import html_add_bootstrap

    html = html_add_bootstrap(get_html_template(), '0.0.5')
    assert '{%bootstrap_css%}' not in html
    assert '{%bootstrap_js%}' not in html
    assert '{%version%}' not in html
    assert '<style>' in html
    assert '</style>' in html
    assert 'Bootstrap v5.1.3' in html
    assert 'LQC v0.0.5' in html
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
uv run pytest tests/test_report_html_assets.py::test_html_add_bootstrap_inlines_vendored_assets -v
```

Expected: FAIL with `ModuleNotFoundError` / `ImportError` (the function does not exist yet).

- [ ] **Step 3: Add imports to `lqc/report_html.py`**

Replace the first line of `lqc/report_html.py`:

```python
import re
```

with:

```python
import os
import re
```

- [ ] **Step 4: Add the `html_add_bootstrap` function**

Append to the end of `lqc/report_html.py`:

```python
_BOOTSTRAP_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'template'
)


def _read_text(path):
    with open(path, 'r', encoding='utf-8') as fh:
        return fh.read()


def html_add_bootstrap(html_string, version):
    css = _read_text(os.path.join(_BOOTSTRAP_DIR, 'bootstrap.min.css'))
    # Escape the sequence that would terminate the inline <script> early.
    js = _read_text(
        os.path.join(_BOOTSTRAP_DIR, 'bootstrap.bundle.min.js')
    ).replace('</script>', '<\\/script>')

    new_html = re.sub(
        r'\{%bootstrap_css%\}',
        lambda _match: '<style>' + css + '</style>',
        html_string
    )
    new_html = re.sub(
        r'\{%bootstrap_js%\}',
        lambda _match: '<script>' + js + '</script>',
        new_html
    )
    new_html = re.sub(
        r'\{%version%\}',
        lambda _match: version,
        new_html
    )
    return new_html
```

> **Why callable replacements:** the vendored JS/CSS contains backslash sequences
> (`\s`, `\d`, …). A string `re.sub` replacement interprets those as escapes and raises
> `re.error: bad escape`. A callable replacement is inserted verbatim, which is the intent.

- [ ] **Step 5: Run the test to verify it passes**

```bash
uv run pytest tests/test_report_html_assets.py::test_html_add_bootstrap_inlines_vendored_assets -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add lqc/report_html.py tests/test_report_html_assets.py
git commit -m "feat: inline vendored Bootstrap assets into the HTML report"
```

---

## Task 4: Implement `html_add_data`

**Files:**
- Modify: `lqc/report_html.py`
- Modify: `tests/test_report_html_assets.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_report_html_assets.py`:

```python
def test_html_add_data_injects_parseable_json():
    import json

    import pandas as pd

    from lqc.report_html import html_add_data

    readstat = pd.DataFrame([{'label': 'Total', 'read_count': 2}])
    html = html_add_data(
        '<html>{%lqc_data%}</html>',
        {'readstat': readstat},
    )
    assert '{%lqc_data%}' not in html
    assert '<script type="application/json" id="lqc-data">' in html
    start = html.index('id="lqc-data">') + len('id="lqc-data">')
    end = html.index('</script>', start)
    payload = json.loads(html[start:end])
    assert payload == {'readstat': [{'label': 'Total', 'read_count': 2}]}


def test_html_add_data_escapes_script_terminator():
    import pandas as pd

    from lqc.report_html import html_add_data

    df = pd.DataFrame([{'label': '</script>'}])
    html = html_add_data('<html>{%lqc_data%}</html>', {'block': df})
    inner = html.split('id="lqc-data">')[1].split('</script>')[0]
    assert '</script>' not in inner
    assert '<\\/script>' in inner
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
uv run pytest tests/test_report_html_assets.py::test_html_add_data_injects_parseable_json tests/test_report_html_assets.py::test_html_add_data_escapes_script_terminator -v
```

Expected: FAIL with `ImportError` (the function does not exist yet).

- [ ] **Step 3: Add `import json`, then the `html_add_data` function**

First, in `lqc/report_html.py`, add `import json` before `import os` (the import block
becomes `import json` / `import os` / `import re`). Then append to the end of `lqc/report_html.py`:

```python
def html_add_data(html_string, tables):
    # to_json -> json.loads round-trip converts numpy scalars to native JSON
    # types so json.dumps never hits a non-serializable int64/float64.
    payload = {
        name: json.loads(table.to_json(orient='records'))
        for name, table in tables.items()
    }
    raw = json.dumps(payload).replace('</', '<\\/')
    script = (
        '<script type="application/json" id="lqc-data">'
        + raw
        + '</script>'
    )
    # Callable replacement: ``script`` can contain backslashes (the ``<\/``
    # escaping above), which a string re.sub replacement would mis-parse.
    return re.sub(r'\{%lqc_data%\}', lambda _match: script, html_string)
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
uv run pytest tests/test_report_html_assets.py::test_html_add_data_injects_parseable_json tests/test_report_html_assets.py::test_html_add_data_escapes_script_terminator -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add lqc/report_html.py tests/test_report_html_assets.py
git commit -m "feat: embed summary tables as JSON in the HTML report"
```

---

## Task 5: Implement `inline_figures`

**Files:**
- Modify: `lqc/report_html.py`
- Modify: `tests/test_report_html_assets.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_report_html_assets.py`:

```python
def test_inline_figures_replaces_png_and_svg(tmp_path):
    from lqc.report_html import inline_figures

    (tmp_path / 'a.png').write_bytes(b'\x89PNG\r\n\x1a\nfake')
    (tmp_path / 'b.svg').write_bytes(b'<svg></svg>')
    html = '<img src="fig/a.png" /><img src="fig/b.svg" />'
    out = inline_figures(html, str(tmp_path))
    assert 'src="fig/' not in out
    assert 'src="data:image/png;base64,' in out
    assert 'src="data:image/svg+xml;base64,' in out


def test_inline_figures_raises_on_missing_file(tmp_path):
    import pytest

    from lqc.report_html import inline_figures

    with pytest.raises(FileNotFoundError):
        inline_figures('<img src="fig/missing.png" />', str(tmp_path))


def test_inline_figures_raises_on_unknown_type(tmp_path):
    import pytest

    from lqc.report_html import inline_figures

    (tmp_path / 'x.foo').write_bytes(b'x')
    with pytest.raises(ValueError):
        inline_figures('<img src="fig/x.foo" />', str(tmp_path))
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
uv run pytest tests/test_report_html_assets.py::test_inline_figures_replaces_png_and_svg tests/test_report_html_assets.py::test_inline_figures_raises_on_missing_file tests/test_report_html_assets.py::test_inline_figures_raises_on_unknown_type -v
```

Expected: FAIL with `ImportError` (the function does not exist yet).

- [ ] **Step 3: Add `import base64`, then the `inline_figures` function**

First, in `lqc/report_html.py`, add `import base64` before `import json` (the import block
becomes `import base64` / `import json` / `import os` / `import re`). Then append to the end
of `lqc/report_html.py`:

```python
_MIME_BY_EXT = {
    '.png': 'image/png',
    '.svg': 'image/svg+xml',
    '.jpg': 'image/jpeg',
    '.jpeg': 'image/jpeg',
    '.gif': 'image/gif',
}


def inline_figures(html_string, fig_dir):
    def _to_data_uri(match):
        rel = match.group(1)
        ext = os.path.splitext(rel)[1].lower()
        if ext not in _MIME_BY_EXT:
            raise ValueError(f'Unsupported figure type: {rel}')
        path = os.path.join(fig_dir, rel)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f'Figure referenced by HTML is missing: {path}'
            )
        with open(path, 'rb') as fh:
            b64 = base64.b64encode(fh.read()).decode('ascii')
        return f'src="data:{_MIME_BY_EXT[ext]};base64,{b64}"'

    return re.sub(r'src="fig/([^"]+)"', _to_data_uri, html_string)
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
uv run pytest tests/test_report_html_assets.py::test_inline_figures_replaces_png_and_svg tests/test_report_html_assets.py::test_inline_figures_raises_on_missing_file tests/test_report_html_assets.py::test_inline_figures_raises_on_unknown_type -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add lqc/report_html.py tests/test_report_html_assets.py
git commit -m "feat: embed report figures as base64 data URIs"
```

---

## Task 6: Wire the new functions into the CLI and export them

**Files:**
- Modify: `lqc/__init__.py`
- Modify: `lqc/cli.py`
- Modify: `tests/test_cli.py`

- [ ] **Step 1: Write the failing end-to-end test**

Append to `tests/test_cli.py`:

```python
def test_main_report_is_self_contained(cs_bam, tmp_path):
    out = tmp_path / 'out'
    rc = main(['-b', cs_bam, '-o', str(out), '-c', 'chr1', '-t', '1'])
    assert rc == 0
    html = (out / 'LQC_report.html').read_text()
    assert 'src="fig/' not in html
    assert 'cdn.jsdelivr.net' not in html
    assert 'img.shields.io' not in html
    assert 'id="lqc-data"' in html
    assert 'data:image/png;base64,' in html
    assert 'data:image/svg+xml;base64,' in html
    assert '{%' not in html
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
uv run pytest tests/test_cli.py::test_main_report_is_self_contained -v
```

Expected: FAIL — the still-not-wired HTML contains `src="fig/`, `cdn.jsdelivr.net`, and `img.shields.io`, and the new functions are not exported/imported.

- [ ] **Step 3: Export the new functions from `lqc/__init__.py`**

In the `from lqc.report_html import (...)` block, change:

```python
from lqc.report_html import (
    html_add_deletion_table,
    html_add_insertion_table,
    html_add_mismatch_table,
    html_add_readstat_table,
    html_add_splice_table,
)
```

to:

```python
from lqc.report_html import (
    html_add_bootstrap,
    html_add_data,
    html_add_deletion_table,
    html_add_insertion_table,
    html_add_mismatch_table,
    html_add_readstat_table,
    html_add_splice_table,
    inline_figures,
)
```

And in `__all__`, add the three names in alphabetical position — change:

```python
    'get_html_template',
    'html_add_deletion_table',
```

to:

```python
    'get_html_template',
    'html_add_bootstrap',
    'html_add_data',
    'html_add_deletion_table',
```

and change:

```python
    'html_add_splice_table',
    'list_bam_contigs',
```

to:

```python
    'html_add_splice_table',
    'inline_figures',
    'list_bam_contigs',
```

- [ ] **Step 4: Import the new functions in `lqc/cli.py`**

In the `from lqc import (...)` block, change:

```python
    get_html_template,
    html_add_deletion_table,
```

to:

```python
    get_html_template,
    html_add_bootstrap,
    html_add_data,
    html_add_deletion_table,
```

and change:

```python
    html_add_splice_table,
    list_bam_contigs,
```

to:

```python
    html_add_splice_table,
    inline_figures,
    list_bam_contigs,
```

- [ ] **Step 5: Chain the new functions before writing the HTML**

In `lqc/cli.py`, locate the end of the report block (right after the `html_add_splice_table(...)` call and before `with open(o_files['html'], 'w')`):

```python
    new_html_string = html_add_splice_table(
        new_html_string, t_splice,
        t_readstat.loc[
            t_readstat['label'] == 'Total',
            'mean_intron_per_read'
        ].values[0]
    )

    with open(o_files['html'], 'w') as f:
        f.write(new_html_string)
```

Replace with:

```python
    new_html_string = html_add_splice_table(
        new_html_string, t_splice,
        t_readstat.loc[
            t_readstat['label'] == 'Total',
            'mean_intron_per_read'
        ].values[0]
    )

    new_html_string = html_add_bootstrap(new_html_string, VERSION)

    new_html_string = html_add_data(
        new_html_string,
        {
            'readstat': t_readstat,
            'insertion': t_insertion,
            'deletion': t_deletion,
            'mismatch': t_mismatch,
            'splice': t_splice,
        }
    )

    new_html_string = inline_figures(new_html_string, o_dirs['fig'])

    with open(o_files['html'], 'w') as f:
        f.write(new_html_string)
```

(`copy_logo(o_dirs['fig'])` already runs earlier in the same block, so `fig/logo_white.svg` and `fig/card.svg` exist before `inline_figures` reads them.)

- [ ] **Step 6: Run the test to verify it passes**

```bash
uv run pytest tests/test_cli.py::test_main_report_is_self_contained -v
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add lqc/__init__.py lqc/cli.py tests/test_cli.py
git commit -m "feat: wire self-contained HTML assembly into the CLI"
```

---

## Task 7: Update the reporting/output docs

**Files:**
- Modify: `docs/reporting-and-output.md`

- [ ] **Step 1: Update the "HTML report contract" section**

In `docs/reporting-and-output.md`, replace the current section:

```markdown
## HTML report contract

- `lqc/template/template.html` contains `{%token%}` placeholders; `lqc/report_html.py`
  (`html_add_*`) is the **only** code that substitutes them. Token groups: `readstat_*`,
  `mismatch_*`, `insertion_*`, `deletion_*`, `intron_*`, `splice_*`.
- Every placeholder in the template must be replaced for a given run; the `html_add_*`
  functions are chained in `lqc/cli.py` and must cover the template's full token set. Add a
  new token only together with the `html_add_*` that fills it.
  *why: unreplaced tokens would surface verbatim as `{%...%}` in the report · when: adding
  a report section · expire: if template rendering moves to a real templating engine.*
- The `'Total'` row is styled specially (`table-secondary`) and also sourced for the
  headline metric values — keep the aggregate row labeled exactly `'Total'`.
```

with:

```markdown
## HTML report contract

- `lqc/template/template.html` contains `{%token%}` placeholders; `lqc/report_html.py`
  (`html_add_*`) is the **only** code that substitutes them. Token groups: `readstat_*`,
  `mismatch_*`, `insertion_*`, `deletion_*`, `intron_*`, `splice_*`, plus the self-contained
  asset tokens `bootstrap_css`, `bootstrap_js`, `lqc_data`, `version`.
- Every placeholder in the template must be replaced for a given run; the `html_add_*`
  functions are chained in `lqc/cli.py` and must cover the template's full token set. Add a
  new token only together with the `html_add_*` that fills it.
  *why: unreplaced tokens would surface verbatim as `{%...%}` in the report · when: adding
  a report section · expire: if template rendering moves to a real templating engine.*
- The `'Total'` row is styled specially (`table-secondary`) and also sourced for the
  headline metric values — keep the aggregate row labeled exactly `'Total'`.
- **Self-contained output:** `lqc/cli.py` finalizes the report with
  `html_add_bootstrap` (inlines the vendored Bootstrap CSS/JS from `lqc/template/`,
  MIT-licensed), `html_add_data` (embeds the five summary tables as JSON in
  `<script type="application/json" id="lqc-data">`), and `inline_figures` (rewrites
  `src="fig/..."` to base64 `data:` URIs). The rendered `LQC_report.html` therefore
  references no `fig/` files and no remote resources and opens offline as a single file.
```

- [ ] **Step 2: Verify no stale wording remains**

Run:

```bash
grep -n "cdn.jsdelivr\|img.shields.io\|requires network" docs/reporting-and-output.md
```

Expected: no matches (exit code 1).

- [ ] **Step 3: Commit**

```bash
git add docs/reporting-and-output.md
git commit -m "docs: document self-contained HTML report tokens and contract"
```

---

## Task 8: Full verification

**Files:** none (verification only)

- [ ] **Step 1: Run the full test suite**

```bash
uv run pytest
```

Expected: all tests PASS (the new `test_report_html_assets.py`, the extended `test_cli.py`, and every pre-existing test).

- [ ] **Step 2: Lint the changed files (no new violations)**

The repo has **104 pre-existing `ruff check` violations** (LOG015, C408, RUF059, SIM115,
PERF102, FLY002, RUF015, B006, F841, SIM113 across `lqc/*.py`) that predate this feature;
whole-repo `uvx ruff check` will NOT be clean until a separate lint-cleanup task lands. Gate
this feature on introducing **no new** violations in the files it changes.

```bash
export UV_CACHE_DIR="$PWD/.uv-cache" UV_TOOL_DIR="$PWD/.uv-tools"
uvx ruff check lqc/report_html.py lqc/cli.py lqc/__init__.py
```

Expected: the reported violation count per file is **≤ the pre-change baseline**
(`report_html.py` ≤ 8, `cli.py` ≤ 40, `__init__.py` = 0), i.e. this change adds none.

- [ ] **Step 3: Manual offline smoke check**

```bash
uv run lqc -b tmp/data/ENCFF417VHJ.chr22.sorted.bam -o tmp/verify_selfcontained -c chr22 -t 1
python - <<'PY'
import pathlib
html = pathlib.Path('tmp/verify_selfcontained/LQC_report.html').read_text()
assert 'src="fig/' not in html
assert 'cdn.jsdelivr.net' not in html
assert 'img.shields.io' not in html
assert 'data:image/png;base64,' in html
assert 'id="lqc-data"' in html
print("self-contained OK:", len(html), "bytes")
PY
```

Expected: prints `self-contained OK: <bytes>`. (Skip if the real BAM at `tmp/data/ENCFF417VHJ.chr22.sorted.bam` is absent; the `cs_bam` end-to-end test already covers the assertion.)

- [ ] **Step 4: Confirm the design's success criteria**

Open `tmp/verify_selfcontained/LQC_report.html` in a browser with the network disabled and confirm all five accordion sections, tables, and figures render identically to the pre-change version (the only visible differences are the two local badges replacing the remote ones). Manually confirm `fig/` and `table/` were still produced unchanged.

---

## Self-Review

**1. Spec coverage** — each spec section maps to a task:

| Spec section | Task(s) |
| --- | --- |
| §1 Vendor Bootstrap 5.1.3 | Task 1 (vendor) + Task 3 (`html_add_bootstrap`) |
| §2 Inline referenced figures as base64 | Task 5 (`inline_figures`) + Task 6 (wiring) |
| §3 Replace remote badge images | Task 2 (template edits) |
| §4 Embed machine-readable JSON | Task 4 (`html_add_data`) + Task 6 (wiring) |
| §5 Wiring / code locations | Task 6 (`__init__.py`, `cli.py`) |
| §6 Symbol contract preserved | Tasks 2–6 keep `html_add_*` as the sole token-filler; no `report_figure.py`/`report_table.py` changes |
| Success criteria 1–4 | Task 6 e2e test + Task 8 (pytest/ruff + offline smoke) |

**2. Placeholder scan** — every code step contains complete code; commands include expected output; no TBD/TODO/"similar to Task N".

**3. Type/signature consistency** — `html_add_bootstrap(html_string, version)`, `html_add_data(html_string, tables)`, and `inline_figures(html_string, fig_dir)` are defined with those exact signatures in Tasks 3–5 and called with matching arguments in Task 6 (cli.py) and tested in the corresponding tasks. The `{%lqc_data%}` / `{%bootstrap_css%}` / `{%bootstrap_js%}` / `{%version%}` token names are identical across template edits (Task 2), function implementations (Tasks 3–4), and assertions.

**Note:** the spec also listed "no `src="https://`" as a success signal; the plan's tests assert the stronger, more precise `cdn.jsdelivr.net` / `img.shields.io` absence plus `src="fig/` absence, and the template edit in Task 2 removes the only three `img.shields.io` images, so no remote image reference remains. Plain `<a href="https://...">` anchors are intentionally preserved (they are navigation, not resource loads).