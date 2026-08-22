# Lint Cleanup + Rule-Set Freeze Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `uv run ruff check` pass (0 violations) and freeze the rule set + ruff version so the lint gate is deterministic across machines and future ruff releases.

**Architecture:** Two parts. (1) Freeze the tooling: pin `ruff==0.16.4` as a dev dependency, add an explicit `[tool.ruff.lint] select` list, set `required-version`, and switch the documented gate from `uvx ruff check` to `uv run ruff check`. (2) Fix all 104 existing violations in code (no `ignore` entries). The 104 violations are distributed across 12 files and 10 rule codes: LOG015 ×38, C408 ×31, RUF059 ×13, F841 ×5, FLY002 ×5, RUF015 ×5, B006 ×3, SIM115 ×2, PERF102 ×1, SIM113 ×1.

**Tech Stack:** Python 3.11, uv, ruff 0.16.4, pytest.

---

## Environment note (READ FIRST)

The DSH sandbox mounts `~/.cache/uv` and `~/.local/share/uv` read-only, so bare `uv`/`uvx` fail with
`Read-only file system`. Export these before **every** `uv`/`uvx` command in this plan:

```bash
export UV_CACHE_DIR=/tmp/lqc-uv-cache \
       UV_TOOL_DIR=/tmp/lqc-uv-tools \
       UV_TOOL_BIN_DIR=/tmp/lqc-uv-toolbin \
       UV_PYTHON_INSTALL_DIR=/tmp/lqc-uv-python
```

On a normal dev machine (or CI) this line is unnecessary — it is sandbox-specific, not project config.

**Important ruff caveat:** ruff 0.16.4 classifies *all* of the fixable violations in this repo as
"unsafe" fixes. `ruff check --fix` alone applies **zero** fixes here ("No fixes available (63 hidden
fixes…)"). Therefore every fix in this plan is an explicit manual edit — do not rely on `--fix` /
`--unsafe-fixes`.

**Commit policy:** Per request, this is a **single combined commit at the end** (Task 15) that bundles
the `.gitignore` update with all code/config/doc changes. Do not commit intermediate tasks.

---

## File map

| File | Violations | Fixes in task |
|---|---|---|
| `pyproject.toml`, `uv.lock`, `AGENTS.md`, `docs/development.md` | (config) | Task 1 |
| `lqc/cli.py` | 40 (LOG015 ×38, C408 ×1, PERF102 ×1) | Task 2 |
| `lqc/cs.py` | 15 (C408 ×8, RUF059 ×3, F841 ×3, RUF015 ×1) | Task 3 |
| `lqc/report_figure.py` | 14 (RUF059 ×10, C408 ×4) | Task 4 |
| `lqc/report_html.py` | 8 (C408 ×7, F841 ×1) | Task 5 |
| `lqc/utils.py` | 7 (F841 ×1, C408 ×1, SIM113 ×1, SIM115 ×1, FLY002 ×1, RUF015 ×2) | Task 6 |
| `lqc/indel.py` | 6 (C408 ×4, B006 ×1, FLY002 ×1) | Task 7 |
| `lqc/mismatch.py` | 6 (C408 ×3, B006 ×2, FLY002 ×1) | Task 8 |
| `lqc/splice.py` | 2 (C408 ×1, FLY002 ×1) | Task 9 |
| `lqc/stat.py` | 2 (RUF015 ×2) | Task 10 |
| `lqc/report_table.py` | 2 (C408 ×2) | Task 11 |
| `lqc/readstat.py` | 1 (FLY002 ×1) | Task 12 |
| `lqc/template/__init__.py` | 1 (SIM115 ×1) | Task 13 |
| (whole repo) | verification | Task 14 |
| `docs/specs/…`, `.gitignore`, all above | design doc + commit | Task 15 |

---

## Task 1: Freeze ruff version + rule set

**Files:**
- Modify: `pyproject.toml:37-44`
- Modify: `AGENTS.md` (Quick Start + constraint #11)
- Modify: `docs/development.md:13,47-48`
- Regenerate: `uv.lock`

- [ ] **Step 1.1: Add ruff dev-dependency + lint config to `pyproject.toml`**

Change the `[tool.ruff]` section (currently only `target-version`) to:

```toml
[tool.ruff]
target-version = "py311"
required-version = "0.16.4"

[tool.ruff.lint]
select = ["E4", "E7", "E9", "F", "I", "B", "C4", "SIM", "RUF", "PERF", "LOG", "FLY"]
```

Change the dependency group:

```toml
[dependency-groups]
dev = ["pytest>=8.0", "ruff==0.16.4"]
```

- [ ] **Step 1.2: Regenerate the lockfile and sync the venv**

Run:

```bash
export UV_CACHE_DIR=/tmp/lqc-uv-cache UV_TOOL_DIR=/tmp/lqc-uv-tools UV_TOOL_BIN_DIR=/tmp/lqc-uv-toolbin UV_PYTHON_INSTALL_DIR=/tmp/lqc-uv-python
uv lock
uv sync
```

Expected: `uv.lock` gains a `ruff` entry; `.venv/bin/ruff` now exists.

Verify:

```bash
uv run ruff --version
```

Expected: `ruff 0.16.4`.

- [ ] **Step 1.3: Update the documented lint command**

In `AGENTS.md`, change the Quick Start bullet and constraint #11:

- `Lint: \`uvx ruff check\`  (append \`--fix\` to auto-fix what it can)` → `Lint: \`uv run ruff check\``
- `11. **Code must pass \`uvx ruff check\`.**` → `11. **Code must pass \`uv run ruff check\`.**`

In `docs/development.md`:

- Line 13: `uvx ruff check                             # lint (add --fix to auto-fix)` → `uv run ruff check                            # lint`
- Update the "Style" paragraph (line 47) that says `adopt \`uvx ruff check\` as the gate` to say `\`uv run ruff check\``.

- [ ] **Step 1.4: Confirm the config resolves and the violation count is unchanged (still 104)**

```bash
uv run ruff check --statistics
```

Expected: still `Found 104 errors` (config freeze does not change the count; the code fixes land next).

---

## Task 2: `lqc/cli.py` — LOG015 ×38, C408 ×1, PERF102 ×1

**Files:**
- Modify: `lqc/cli.py:52-53` (add logger), `:78`, `:87-89`, and every `logging.<call>` site listed below

- [ ] **Step 2.1: Add a module logger**

Insert after `matplotlib.use('Agg')` (line 52), before `def stat_bam`:

```python
matplotlib.use('Agg')

logger = logging.getLogger(__name__)
```

- [ ] **Step 2.2: Route the 38 root-logger calls through `logger`**

Run this one-line substitution (matches only `logging.debug/info/error/warning(`, so it leaves
`logging.basicConfig(`, `logging.INFO`, and `logging.DEBUG` untouched):

```bash
sed -i -E 's/logging\.(debug|info|error|warning)\(/logger.\1(/g' lqc/cli.py
```

Verify — the count of remaining root-logger calls must be 0 (only `basicConfig`/`INFO`/`DEBUG` may mention `logging.`):

```bash
grep -cE 'logging\.(debug|info|error|warning)\(' lqc/cli.py   # expect 0
```

- [ ] **Step 2.3: C408 (line 78)**

Change:

```python
    stat_list = list()
```

to:

```python
    stat_list = []
```

- [ ] **Step 2.4: PERF102 (lines 87-89)**

Change:

```python
def build_directories(dir_dict):
    for a, b in dir_dict.items():
        os.makedirs(b, exist_ok = True)
```

to:

```python
def build_directories(dir_dict):
    for b in dir_dict.values():
        os.makedirs(b, exist_ok = True)
```

- [ ] **Step 2.5: Verify this file is clean**

```bash
uv run ruff check lqc/cli.py
```

Expected: no warnings ("All checks passed!").

---

## Task 3: `lqc/cs.py` — C408 ×8, RUF059 ×3, F841 ×3, RUF015 ×1

**Files:**
- Modify: `lqc/cs.py:29,55,78,125,148,168,170-171,182-185,204-206,207,265,505,543`

- [ ] **Step 3.1: C408 `list()` → `[]` (8 sites)**

Replace each exact line:

```python
    cslist = list()              # line 29
    cigarlist = list()           # line 55
    mdlist = list()              # line 78
    new_cigar = list()           # line 168
            new_items = list()   # line 207
    cigar_pos = list()      # sequence position by cigar   # line 265
        read_location = list()   # line 505
            new_cs = list()      # line 543
```

with the same line using `[]` instead of `list()` (keep indentation and the trailing comment on line 265).

- [ ] **Step 3.2: RUF059 — unused unpacked variables**

Line 125:

```python
        low, high, cigar_mark, cigar_value = item
```

→

```python
        _, _, cigar_mark, cigar_value = item
```

Line 148:

```python
        md_value, md_type, md_len = item
```

→

```python
        _, md_type, md_len = item
```

- [ ] **Step 3.3: F841 — remove dead assignments**

Delete lines 170-171 (the `cigar_start` / `cigar_end` assignments; the loop only ever uses `cigar[i][6]`/`[7]`):

```python
        cigar_start = cigar[i][4]
        cigar_end = cigar[i][5]
```

Delete lines 204-206 (the `sorted(...)` result is never referenced — the loop below iterates `inside_mismatches`):

```python
            inside_mismatch_ref_pos = sorted(
                inside_mismatches, key = lambda a: a[3]
            )
```

- [ ] **Step 3.4: RUF015 (lines 182-185)**

Change:

```python
                deletion_item = [
                    a for a in deletions
                    if a[3] == md_start
                ][0]
```

to:

```python
                deletion_item = next(
                    a for a in deletions
                    if a[3] == md_start
                )
```

- [ ] **Step 3.5: Verify this file is clean**

```bash
uv run ruff check lqc/cs.py
```

Expected: "All checks passed!".

---

## Task 4: `lqc/report_figure.py` — RUF059 ×10, C408 ×4

**Files:**
- Modify: `lqc/report_figure.py:172,240,313,372,448,548,597,653,656,662,711,760,763,777`

- [ ] **Step 4.1: RUF059 — rename unused `axes` to `_` (ONLY these 10 lines)**

Change `fig, axes = plt.subplots(` → `fig, _ = plt.subplots(` at exactly these line numbers:
**172, 240, 313, 372, 448, 548, 597, 662, 711, 777**.

> ⚠️ Do **not** touch lines 114 and 505 — there `axes` is actually used (`axes.bar(...)`).

- [ ] **Step 4.2: C408 `list()` → `[]` (4 sites)**

```python
    mislist = list()   # line 653  ->  mislist = []
        alist = list()  # line 656  ->  alist = []
    splist = list()    # line 760  ->  splist = []
        alist = list()  # line 763  ->  alist = []
```

- [ ] **Step 4.3: Verify**

```bash
uv run ruff check lqc/report_figure.py
```

Expected: "All checks passed!".

---

## Task 5: `lqc/report_html.py` — C408 ×7, F841 ×1

**Files:**
- Modify: `lqc/report_html.py:8,95-98,101,141,146,171,230,291`

- [ ] **Step 5.1: F841 — delete the unused `bins` list (lines 95-98)**

```python
    bins = ['[0.0,0.1)', '[0.1,0.2)', '[0.2,0.3)',
            '[0.3,0.4)', '[0.4,0.5)', '[0.5,0.6)',
            '[0.6,0.7)', '[0.7,0.8)', '[0.8,0.9)',
            '[0.9,1.0]']
```

- [ ] **Step 5.2: C408 `list()` → `[]` (7 sites)**

Replace `rowstring_list = list()` at lines **8, 101, 171, 230, 291**, and:

```python
    mistype_list1 = list()   # line 141  ->  mistype_list1 = []
    mistype_list2 = list()   # line 146  ->  mistype_list2 = []
```

- [ ] **Step 5.3: Verify**

```bash
uv run ruff check lqc/report_html.py
```

Expected: "All checks passed!".

---

## Task 6: `lqc/utils.py` — F841 + C408 (line 47), SIM113, SIM115, FLY002, RUF015 ×2

**Files:**
- Modify: `lqc/utils.py:41-76` (`check_bam_with_cs_or_md`) and `:79-149` (`write_readcs`)

- [ ] **Step 6.1: Fix `check_bam_with_cs_or_md` (F841 + C408 + SIM113)**

Replace the loop header block:

```python
    bam = pysam.AlignmentFile(bam_file, file_read)
    i = 0
    read_cs_md = list()
    bam_type = None
    for read in bam:
        i += 1
        if i >= 10:
            break
        else:
```

with (removes dead `read_cs_md`, uses `enumerate(start=1)` to preserve the exact "9 reads, break on the 10th" behavior):

```python
    bam = pysam.AlignmentFile(bam_file, file_read)
    bam_type = None
    for i, read in enumerate(bam, start = 1):
        if i >= 10:
            break
        else:
```

The rest of the function body (lines 53-76) is unchanged.

- [ ] **Step 6.2: Rewrite `write_readcs` (SIM115 + FLY002 + RUF015 ×2)**

Replace the entire `write_readcs` function (lines 79-149) with:

```python
def write_readcs(bam_file,
                 genome_file,
                 output_file,
                 method = 'cs'):
    assert method in ['cs', 'MD', 'both'],\
        "method should be either: cs, MD, both."
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)

    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)
    else:
        pass

    with open(output_file, 'w') as output:
        output.write(
            'read_name\tcontig\tlow\thigh\tcs_mark\tcs_value\n'
        )
        for read in bam:
            strand = '-' if read.is_reverse else '+'
            if method == 'cs':
                # there're cs tags in the bam file
                cs_string = next(
                    a[1] for a in read.tags
                    if a[0] == 'cs'
                )
                cs = CS.from_cs_tag_string(
                    cs_tag_string = cs_string,
                    contig = read.reference_name,
                    start_pos = read.reference_start,
                    strand = strand
                )
            else:
                # there's no cs tag in the bam file, and there're MD tags.
                read_seq = read.query_sequence
                ref_seq = genome.fetch(
                    read.reference_name,
                    read.reference_start,
                    read.reference_end
                )
                md_string = next(
                    a[1] for a in read.tags
                    if a[0] == 'MD'
                )
                cs = CS.from_cigar_string(
                    cigar_string = read.cigarstring,
                    md_string = md_string,
                    read_seq = read_seq,
                    ref_seq = ref_seq,
                    contig = read.reference_name,
                    start_pos = read.reference_start,
                    strand = strand
                )
            # read cs
            cs_list = cs.get_contig_position()
            output.writelines('\t'.join(
                        [read.query_name,
                         read.reference_name] +
                        [f'{a}' for a in line]
                    ) + '\n' for line in cs_list)
    bam.close()
    if method not in ['cs', 'both']:
        genome.close()
    else:
        pass
```

> Changes vs original: `open(...)` wrapped in `with` (no trailing `output.close()`), header tab-join
> replaced by the literal string, and the two `[…][0]` list slices replaced with `next(…)`.

- [ ] **Step 6.3: Verify**

```bash
uv run ruff check lqc/utils.py
```

Expected: "All checks passed!".

---

## Task 7: `lqc/indel.py` — C408 ×4, B006 ×1, FLY002 ×1

**Files:**
- Modify: `lqc/indel.py:22-24,133-134,163,169`

- [ ] **Step 7.1: C408 `list()`/`dict()` → `[]`/`{}` (4 sites)**

```python
        self._indel_strings = list()   # line 22  ->  self._indel_strings = []
        self._indel_index = dict()     # line 23  ->  self._indel_index = {}
        self._indels = list()          # line 24  ->  self._indels = []
        other_new_idx = list()         # line 169 ->  other_new_idx = []
```

- [ ] **Step 7.2: B006 (lines 133-134)**

Change:

```python
    def get_location_bin_count(self,
                               cuts = [0, 0.25, 0.5, 0.75, 1]):
        edges, hist = self.get_location_histogram(cuts = cuts)
```

to:

```python
    def get_location_bin_count(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        edges, hist = self.get_location_histogram(cuts = cuts)
```

- [ ] **Step 7.3: FLY002 (line 163)**

Change:

```python
            ' '.join([self.label, other.label])
```

to:

```python
            f'{self.label} {other.label}'
```

- [ ] **Step 7.4: Verify**

```bash
uv run ruff check lqc/indel.py
```

Expected: "All checks passed!".

---

## Task 8: `lqc/mismatch.py` — C408 ×3, B006 ×2, FLY002 ×1

**Files:**
- Modify: `lqc/mismatch.py:21-22,78-79,111-112,123,149`

- [ ] **Step 8.1: C408 (3 sites)**

```python
        self._mismatch_types = list()   # line 21 ->  self._mismatch_types = []
        self._mismatch_index = dict()   # line 22 ->  self._mismatch_index = {}
        locations = list()              # line 123 -> locations = []
```

- [ ] **Step 8.2: B006 — `get_location_bin_count` (lines 78-79)**

Change:

```python
    def get_location_bin_count(self,
                               cuts = [0, 0.25, 0.5, 0.75, 1]):
        bin_count = Counter()
```

to:

```python
    def get_location_bin_count(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        bin_count = Counter()
```

- [ ] **Step 8.3: B006 — `get_location_bin_count_by_type` (lines 111-112)**

Change:

```python
    def get_location_bin_count_by_type(self,
                                       cuts = [0, 0.25, 0.5, 0.75, 1]):
        type_bin_count_dict = defaultdict(
```

to:

```python
    def get_location_bin_count_by_type(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        type_bin_count_dict = defaultdict(
```

- [ ] **Step 8.4: FLY002 (line 149)**

Change:

```python
        newMis = type(self)(' '.join([self.label, other.label]))
```

to:

```python
        newMis = type(self)(f'{self.label} {other.label}')
```

- [ ] **Step 8.5: Verify**

```bash
uv run ruff check lqc/mismatch.py
```

Expected: "All checks passed!".

---

## Task 9: `lqc/splice.py` — C408 ×1, FLY002 ×1

**Files:**
- Modify: `lqc/splice.py:44,89`

- [ ] **Step 9.1: C408 (line 44)**

```python
        new_dict = dict()
```
→
```python
        new_dict = {}
```

- [ ] **Step 9.2: FLY002 (line 89)**

```python
            ' '.join([self.label, other.label])
```
→
```python
            f'{self.label} {other.label}'
```

- [ ] **Step 9.3: Verify**

```bash
uv run ruff check lqc/splice.py
```

Expected: "All checks passed!".

---

## Task 10: `lqc/stat.py` — RUF015 ×2

**Files:**
- Modify: `lqc/stat.py:38-41,54-57`

- [ ] **Step 10.1: RUF015 (cs tag)**

Change:

```python
            cs_string = [
                a[1] for a in read.tags
                if a[0] == 'cs'
            ][0]
```

to:

```python
            cs_string = next(
                a[1] for a in read.tags
                if a[0] == 'cs'
            )
```

- [ ] **Step 10.2: RUF015 (MD tag)**

Change:

```python
            md_string = [
                a[1] for a in read.tags
                if a[0] == 'MD'
            ][0]
```

to:

```python
            md_string = next(
                a[1] for a in read.tags
                if a[0] == 'MD'
            )
```

- [ ] **Step 10.3: Verify**

```bash
uv run ruff check lqc/stat.py
```

Expected: "All checks passed!".

---

## Task 11: `lqc/report_table.py` — C408 ×2

**Files:**
- Modify: `lqc/report_table.py:24,90`

- [ ] **Step 11.1: C408**

```python
    row_list = list()    # line 24 ->  row_list = []
    data_list = list()   # line 90 ->  data_list = []
```

- [ ] **Step 11.2: Verify**

```bash
uv run ruff check lqc/report_table.py
```

Expected: "All checks passed!".

---

## Task 12: `lqc/readstat.py` — FLY002 ×1

**Files:**
- Modify: `lqc/readstat.py:251`

- [ ] **Step 12.1: FLY002**

```python
            ' '.join([self.label, other.label])
```
→
```python
            f'{self.label} {other.label}'
```

- [ ] **Step 12.2: Verify**

```bash
uv run ruff check lqc/readstat.py
```

Expected: "All checks passed!".

---

## Task 13: `lqc/template/__init__.py` — SIM115 ×1

**Files:**
- Modify: `lqc/template/__init__.py:8-17`

- [ ] **Step 13.1: SIM115 — use a context manager**

Change:

```python
def get_html_template():
    html = open(
        os.path.join(file_dir, "template.html"),
        "r"
    )
    html_string = "\n".join(
        [read.strip() for read in html]
    )
    html.close()
    return html_string
```

to:

```python
def get_html_template():
    with open(
        os.path.join(file_dir, "template.html"),
        "r"
    ) as html:
        html_string = "\n".join(
            [read.strip() for read in html]
        )
    return html_string
```

- [ ] **Step 13.2: Verify**

```bash
uv run ruff check lqc/template/__init__.py
```

Expected: "All checks passed!".

---

## Task 14: Full verification

- [ ] **Step 14.1: Whole-repo lint must be clean**

```bash
uv run ruff check
```

Expected: `All checks passed!` (0 errors; previously 104).

- [ ] **Step 14.2: Tests must still pass**

```bash
uv run pytest -q
```

Expected: `86 passed`. (Baseline before changes was 86 passed.)

---

## Task 15: Design doc + single combined commit

**Files:**
- Create: `docs/specs/2026-08-22-lint-rule-freeze-design.md`
- Include in commit: `.gitignore`, `pyproject.toml`, `uv.lock`, all 12 modified `lqc/*.py`
  files, `AGENTS.md`, `docs/development.md`, `docs/plans/2026-08-22-lint-rule-freeze.md`, and the new design doc.

- [ ] **Step 15.1: Write the design doc**

Create `docs/specs/2026-08-22-lint-rule-freeze-design.md`:

```markdown
# Lint Rule-Set Freeze Design

## Problem Statement
`uvx ruff check` pulls the latest ruff (0.16.4), whose default rule set is broad and
version-dependent. Because the repo pinned nothing beyond `target-version`, it reported 104
violations and AGENTS.md constraint #11 was unmet. Any future ruff release could silently change
the rule set and re-break the gate.

## Proposed Solution
Freeze the tool and the rules: pin `ruff==0.16.4` as a dev dependency, add an explicit
`[tool.ruff.lint] select` list, set `required-version = "0.16.4"`, and switch the documented gate
to `uv run ruff check`. Then fix all 104 existing violations in code so the gate is clean with no
`ignore` entries.

## Detailed Design
- **Rule set** (curated, a subset of ruff 0.16.4's default): `E4, E7, E9, F, I, B, C4, SIM, RUF,
  PERF, LOG, FLY`. No `ignore` entries.
- **Version pin**: `ruff==0.16.4` in `[dependency-groups] dev` + `uv.lock`, guarded by
  `required-version = "0.16.4"`.
- **Code fixes** (104 total): LOG015 → module logger in `cli.py`; C408 → `[]`/`{}` literals;
  RUF059 → rename unused unpacked vars to `_`; F841 → remove verified-dead assignments; FLY002 →
  f-strings/literals; RUF015 → `next(...)`; B006 → `None` sentinel for mutable `cuts` defaults;
  SIM115 → `with open(...)`; PERF102 → `.values()`; SIM113 → `enumerate(..., start=1)`.

## Success Criteria
`uv run ruff check` reports 0 errors and `uv run pytest` reports 86 passed, committed as a single
change.

## Out of Scope
`ruff format`, line-length (`E501`), docstring (`D`), and type-annotation rules; CI integration.
```

- [ ] **Step 15.2: Confirm `.gitignore` is staged as-is**

```bash
git status --short
```

Expected: `.gitignore` listed as `modified` alongside the other changed files (it is the user's
own update — stage it verbatim, do not edit its contents).

- [ ] **Step 15.3: Single commit**

```bash
git add -A
git commit -m "chore: freeze ruff ruleset/version and resolve all 104 lint violations"
```

- [ ] **Step 15.4: Post-commit sanity**

```bash
git status --short
uv run ruff check
```

Expected: `git status` clean after the commit; `ruff check` → `All checks passed!`.

---

## Scope Expansion (user-approved): 104 → 153

The original `select` written in Task 1 used group prefixes, which enable **7 rule types the ruff
0.16.4 *default* excludes**. The user chose to keep the broad select and fix these too, raising the
target from 104 to **153** violations. The additions, per file (merged into the tasks below):

- **`lqc/cli.py`** — PERF401 (line 82-84 `get_stat_list` loop → list comprehension, together with C408 at 78); B007 on line 88 is the same line as PERF102 (`dir_dict.items()` → `dir_dict.values()`).
- **`lqc/cs.py`** — B905 ×3 (`zip(cs_mark, cs_value)`, `zip(cigar_mark, cigar_num)`, `zip(read_loc, relative_pos)` → add `strict = True`); RUF005 ×5 (`item + [...]` → `[*item, ...]` at 134/154/203, `cigar[i][0:3] + [x]` → `[*cigar[i][0:3], x]` at 187, `[a,b] + b[2:]` → `[a, b, *b[2:]]` at 549); B007 (loop var → `_` at 267 and 488/508 → `for _, _, c, d in ...`).
- **`lqc/indel.py`** — B007 (`for iidx, ilen, loc in self._indels:` → `iidx, _, _`); E721 (`type(other) == type(self)` → `isinstance(other, type(self))`).
- **`lqc/mismatch.py`** — B905 ×2 (`zip(self._mismatch_types, new_types)` → `strict = True`); E721 (`type(other) == type(self)` → `isinstance`).
- **`lqc/readstat.py`** — E721 (`type(other) == type(self)` → `isinstance`).
- **`lqc/splice.py`** — E721 (`type(other) == type(self)` → `isinstance`).
- **`lqc/stat.py`** — B007 ×3 (`for a, b, c, d in cs.get_insertions/deletions/mismatches(...)` → `for a, _, _, d in ...`); SIM108 ×3 (strand if/else → ternary).
- **`lqc/report_figure.py`** — C416 (`accumulate([a for a in lengths])` → `accumulate(lengths)`); PERF401 ×2 (`alist.append(count)` inner loop → list comprehension); SIM108 (`ratio` if/else → ternary).
- **`lqc/report_html.py`** — B007 ×5 (`for ri, row in ...iterrows():` → `for _, row in ...`); PERF401 ×2 (`mistype_list1/2` append loops → comprehensions).
- **`lqc/report_table.py`** — PERF401 ×2 (`data_list` append loops → nested comprehension + `extend`); RUF005 ×2 (`['label','bin'] + mistypes` → `['label','bin',*mistypes]`; `['label'] + sptypes` → `['label',*sptypes]`).
- **`tests/conftest.py`** — PERF401 ×1 (`records.append(...)` loop → list comprehension) — **new Task 13b**.

`zip(..., strict = True)` is behavior-preserving here: every zipped pair is derived from the same
source (same `cs`/cigar string, or `read_loc`/`relative_pos` both built from `get_relative_position()`),
so lengths are always equal. `isinstance(other, type(self))` is the standard `__add__` type guard and
accepts subclasses (the intended semantics).

## Self-review notes

- **Spec coverage:** all 104 violations (10 rule codes / 12 files) map to Tasks 2-13; drift
  prevention (`select`, `required-version`, dev-dep pin, `uv run` gate) → Task 1; verification →
  Task 14; design doc + combined commit → Task 15.
- **Type/name consistency:** the `None`-sentinel names (`cuts`), the `next(...)` generator form, and
  `logger`/`_` renames are consistent across every file that uses them.
- **Behavior preservation:** `enumerate(..., start=1)` reproduces the original 9-read scan;
  `next(...)` over `[…][0]` keeps the required-tag path (no default supplied, matching the original
  "tag must exist" assumption); B006 `None` sentinels rebuild the exact same `[0, 0.25, 0.5, 0.75, 1]`
  default on each call.