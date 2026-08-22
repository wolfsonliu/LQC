# Lint Rule-Set Freeze Design

## Problem Statement

`uvx ruff check` pulls the latest ruff (0.16.4), whose default rule set is broad and
version-dependent. Because the repo pinned nothing beyond `target-version`, it reported 104
violations (under ruff's default select) and AGENTS.md constraint #11 was unmet. Any future ruff
release could silently change the rule set and re-break the gate.

## Proposed Solution

Freeze the tool and the rules: pin `ruff==0.16.4` as a dev dependency, set
`required-version = "0.16.4"`, and switch the documented gate from `uvx ruff check` to
`uv run ruff check`. Use an explicit `[tool.ruff.lint] select` list, then fix all resulting
violations in code so the gate is clean with no `ignore` entries.

## Detailed Design

- **Rule set** (explicit group-prefix select): `E4, E7, E9, F, I, B, C4, SIM, RUF, PERF, LOG, FLY`.
  This is broader than ruff 0.16.4's default (which omits some rules inside B/C4/SIM/RUF/PERF/E7),
  so it surfaced 49 additional violations beyond the initial 104 — **153 total**. The user chose to
  keep this broader set and fix all 153.
- **Version pin**: `ruff==0.16.4` in `[dependency-groups] dev` + `uv.lock`, guarded by
  `required-version = "0.16.4"`.
- **Code fixes (153 total)**:
  - LOG015 → module logger in `cli.py`; C408 → `[]`/`{}` literals; RUF059 → `_` for unused
    unpacked variables; F841 → remove dead assignments; FLY002 → f-strings/literals; RUF015 →
    `next(...)`; B006 → `None` sentinel; SIM115 → `with open(...)`; PERF102 → `.values()`;
    SIM113 → `enumerate(..., start=1)`.
  - Broader-select additions: B007 (unused loop vars → `_`), PERF401 (append-loops →
    comprehensions/`extend`), RUF005 (list concat → `*` unpacking), B905 (`zip(..., strict=True)`),
    SIM108 (if/else → ternary), E721 (`type(x)==type(y)` → `isinstance`), C416 (identity
    comprehension → iterable).

## Success Criteria

`uv run ruff check` reports 0 errors and `uv run pytest` reports 86 passed, committed as a single
change.

## Out of Scope

`ruff format`, line-length (`E501`), docstring (`D`), and type-annotation rules; CI integration.