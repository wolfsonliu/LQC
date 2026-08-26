# M3 — report-phase polish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: subagent-driven-development or executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Resolve the two §4.4 items — (1) the per-contig figure cost policy (decided by measurement now that the full run completes post-M2) and (2) the `convert_reverse_complement` micro-fix — without touching the non-negotiable byte-parity bar.

**Measurement (post-M2 full-genome run, `-t 4`, `tmp/m2/full`):** wall 8:26.86 total — stat 2:38, reduce/`'Total'` 0:28, **figures 5:16 (62%)**, html ~0. `generate_multiple_figs` renders `1 aggregate + 1 Total + 24 per-contig` figures for each of ~15 plot types (≈390 figures × PNG+PDF ≈ 780 files), sequentially in the main thread.

**Figure-policy decision (Task 1):** **Option 1 — keep current behavior.** Rationale: the 15 plot functions in `src/lqc/report_figure.py` mix `plt.*`/`pyplot` module-level state (`plt.subplots`, `plt.tight_layout`, `plt.close`), so thread-safety would require a full object-API refactor of all 15; process-parallelism would require pickling the (now-numpy, so cheap — but still large: `l_readstat`+`l_mismatch`+… across 24 contigs) stat objects into workers, plus filename/`savefig` marshaling. Both are disproportionate to a lower-priority polish item, and either risks the byte-identical `fig/*.png` bar. The run now **completes in 8:27** (previously OOM-killed), so 5:16 of figures is acceptable. This decision is recorded here and is the substantive M3 deliverable — no code change.

**Micro-fix (Task 2):** replace `convert_reverse_complement`'s `convert_complement(string)[::-1]` (translate, then extra function call, then reverse) with a direct `string.translate(_COMPLEMENT_TABLE)[::-1]`, i.e. one translate with the reverse slice, dropping the redundant `convert_complement` indirection. Byte-identical. `convert_complement`/`convert_reverse_complement` stay pure helpers.

---

## File Structure

- **Modify `src/lqc/utils.py`**: inline the translate in `convert_reverse_complement`.
- **Modify `tests/test_utils.py`**: add a test pinning `convert_reverse_complement` equivalence and case handling.

---

## Task 1: Figure policy — keep (no code change)

- [ ] **Step 1:** Confirm the measurement above is in `tmp/m2/full.time.txt`/`tmp/m2/full.run.log`. No code changes.

This task's only output is this decision record. (No test, no commit for this task alone; it ships together with Task 2's commit.)

---

## Task 2: `convert_reverse_complement` micro-fix

**Files:**
- Modify: `src/lqc/utils.py`
- Test: `tests/test_utils.py` (create if absent)

- [ ] **Step 1: Add a characterization test**

Add to `tests/test_utils.py` (create the file if it does not exist):

```python
from lqc.utils import convert_reverse_complement


def test_convert_reverse_complement():
    assert convert_reverse_complement('ag') == 'ct'
    assert convert_reverse_complement('atcg') == 'cgat'
    assert convert_reverse_complement('n-') == '-n'
    assert convert_reverse_complement('') == ''


def test_convert_reverse_complement_uppercase_unchanged():
    # uppercase is not in _COMPLEMENT_TABLE, so it passes through unchanged
    assert convert_reverse_complement('AG') == 'GA'
```

- [ ] **Step 2: run to verify they pass against current code** — `.venv/bin/python -m pytest tests/test_utils.py -q` (characterization; current impl already satisfies them).

- [ ] **Step 3: implement**

In `src/lqc/utils.py`, replace:

```python
def convert_reverse_complement(string):
    return convert_complement(string)[::-1]
```

with:

```python
def convert_reverse_complement(string):
    return string.translate(_COMPLEMENT_TABLE)[::-1]
```

- [ ] **Step 4: run tests + ruff** — `.venv/bin/python -m pytest tests/test_utils.py -q` and `.venv/bin/python -m ruff check src/lqc/utils.py tests/test_utils.py` → all pass.

- [ ] **Step 5: full suite** — `.venv/bin/python -m pytest -q` (expect **142 passed** = 140 + 2 new), `.venv/bin/python -m ruff check src tests` clean.

- [ ] **Step 6: commit**

```bash
git add src/lqc/utils.py tests/test_utils.py
git commit -m "perf: inline translate in convert_reverse_complement"
```

---

## Task 3: Verification (byte-parity + full re-run)

- [ ] **Step 1: chr22 byte-parity** — re-run `.venv/bin/lqc -b tmp/full_bam/ENCFF417VHJ.sorted.bam -c chr22 -o tmp/m3/chr22 -t 1 --output-cs` and `cmp` all `table/*.txt`, `read.cs`, `LQC_report.html`, `fig/*.png` against `tmp/verify/chr22` → all match (the reverse-complement path is exercised by reverse-strand reads on chr22).

- [ ] **Step 2: full-genome re-run** — `bash tmp/run_lqc_full.sh tmp/full_bam/ENCFF417VHJ.sorted.bam tmp/m3/full 4` → exit 0, peak RSS recorded, wall recorded.

- [ ] **Step 3: record results** in this plan.

---

## Self-Review

- **Spec coverage:** §4.4 — figure policy resolved (keep, measured, documented) and `convert_reverse_complement` micro-fix (one translate + reverse slice, no redundant call; stays a pure helper). Covered.
- **Placeholder scan:** no TBD/TODO; complete code in every step.
- **Key risk:** none — the micro-fix is byte-identical (chr22 PNG + table byte-parity guards it), and the figure decision introduces no code change.