"""Keep each figure notebook's filename in sync with the panels it writes.

A notebook is named after the published panels it generates, so the file itself
says what it is for::

    Fig4K_S12-S14__Interaction_Doxo1.ipynb

The name is derived, never typed by hand: this script reads every ``fu.savefig``
panel name out of a notebook, resolves it through ``PAPER_PANELS`` in
``_figutils.py`` (the single source of truth for published names), and builds the
expected filename from the result. Panels mapped to ``None`` are not published and
so do not appear in the name; a notebook that writes no figure at all (the data
preparation notebook) keeps whatever it is called.

Usage::

    python _panel_names.py           # check every notebook, report mismatches
    python _panel_names.py --fix     # rename the mismatched notebooks (git mv)
    python _panel_names.py --table   # notebook -> panels table, for the legend docs

Run it after adding, removing or renumbering a panel. ``PAPER_PANELS`` stays the
only place a published name is edited; this script propagates the consequence to
the filenames.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGUTILS = HERE / "_figutils.py"

# Notebooks that write no published panel keep their own name.
NO_PANEL_KEEP = {"00_Prepare_data.ipynb"}

# The separator between the panel spec and the topic part of a filename.
SEP = "__"


# ---------------------------------------------------------------------------
# PAPER_PANELS, read without importing _figutils (which pulls in scanpy)
# ---------------------------------------------------------------------------
def load_paper_panels(path: Path = FIGUTILS) -> dict:
    """Return the PAPER_PANELS dict from _figutils.py source, no import needed."""
    src = path.read_text()
    start = src.index("PAPER_PANELS = {")
    block = src[start:]
    block = block[: block.index("\n}") + 2]
    ns: dict = {}
    exec(block, ns)                       # a literal dict; nothing else runs
    return ns["PAPER_PANELS"]


# ---------------------------------------------------------------------------
# What a notebook writes
# ---------------------------------------------------------------------------
_STR_RE = re.compile(r"""f?(['"])((?:[^'"\\\n]|\\.)*?)\1""")
_SAVEFIG_RE = re.compile(r"""fu\.savefig\(\s*f?(['"])((?:[^'"\\\n]|\\.)*?)\1""")
_PLACEHOLDER_RE = re.compile(r"\{[^{}]*\}")


def notebook_code(nb_path: Path) -> str:
    """All code-cell source of a notebook, concatenated."""
    cells = json.loads(nb_path.read_text()).get("cells", [])
    return "\n".join(
        "".join(c.get("source", [])) for c in cells if c.get("cell_type") == "code"
    )


def _as_pattern(literal: str) -> re.Pattern:
    """Turn a (possibly f-string) literal into a regex matching the panel names it can take.

    ``Fig2e1_doxo1_vs_pseudotime_{mkey}_{day}`` becomes a pattern matching every
    ``PAPER_PANELS`` key of that shape, so the loop variables need not be known.
    """
    parts = _PLACEHOLDER_RE.split(literal)
    return re.compile(".+?".join(re.escape(p) for p in parts) + r"\Z")


def written_panels(nb_path: Path, panels: dict) -> tuple[list[str], list[str]]:
    """Panels a notebook writes.

    Returns ``(published, unmapped)``: ``published`` are the names the files are
    actually written under (PAPER_PANELS values, ``None`` entries dropped);
    ``unmapped`` are savefig names with no PAPER_PANELS entry, which savefig writes
    verbatim. A savefig name built from a variable is still found, because every
    string literal in the notebook is matched against the known panel keys.
    """
    code = notebook_code(nb_path)
    keys: set[str] = set()
    unmapped: list[str] = []

    # Literals passed straight to savefig are panels by definition.
    for _, lit in _SAVEFIG_RE.findall(code):
        if _PLACEHOLDER_RE.search(lit):
            continue                      # resolved below, against the known keys
        if lit in panels:
            keys.add(lit)
        else:
            unmapped.append(lit)

    # Any literal in the notebook that matches a known panel key counts too. This
    # is what catches a name passed through a variable (the Interaction notebook
    # hands savefig an `fname` argument) without having to trace the variable.
    # A literal that opens with a placeholder is skipped: f"{v:,}" would compile to
    # the pattern ".+?" and match every panel in the table.
    for _, lit in _STR_RE.findall(code):
        if not _PLACEHOLDER_RE.split(lit)[0]:
            continue
        pat = _as_pattern(lit)
        keys.update(k for k in panels if pat.match(k))

    published = sorted({panels[k] for k in keys if panels.get(k) is not None})
    return published, sorted(set(unmapped))


# ---------------------------------------------------------------------------
# Filename construction
# ---------------------------------------------------------------------------
_MAIN4_RE = re.compile(r"\AFig4([A-Z])_")
_MAIN5_RE = re.compile(r"\AFig5([A-Z])_")      # lettered panel, e.g. Fig5F_
_MAIN5_ANY_RE = re.compile(r"\AFig5")           # any Figure 5 panel, lettered or not
_SUPP_RE = re.compile(r"\AFigS(\d+)_")


def _collapse(nums: list[int]) -> str:
    """[1,2,3,5,9,10] -> 'S01-S03_S05_S09-S10' (consecutive runs collapsed)."""
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"S{nums[i]:02d}" if i == j else f"S{nums[i]:02d}-S{nums[j]:02d}")
        i = j + 1
    return "_".join(out)


def panel_spec(published: list[str]) -> str:
    """The panel part of a filename, e.g. 'Fig4L_Fig5_S15'."""
    letters = sorted({m.group(1) for p in published if (m := _MAIN4_RE.match(p))})
    letters5 = sorted({m.group(1) for p in published if (m := _MAIN5_RE.match(p))})
    has5 = any(_MAIN5_ANY_RE.match(p) for p in published)
    supp = sorted({int(m.group(1)) for p in published if (m := _SUPP_RE.match(p))})

    tokens = []
    if letters:
        tokens.append("Fig4" + "".join(letters))
    if letters5:
        tokens.append("Fig5" + "".join(letters5))
    elif has5:
        tokens.append("Fig5")
    if supp:
        tokens.append(_collapse(supp))
    return "_".join(tokens)


def topic_of(nb_path: Path) -> str:
    """The descriptive half of a notebook name, independent of the panel spec."""
    stem = nb_path.stem
    if SEP in stem:
        return stem.split(SEP, 1)[1]
    return re.sub(r"\AFig\d+_", "", stem)      # legacy 'Fig2_Pseudotime' naming


def expected_name(nb_path: Path, published: list[str]) -> str:
    """The filename this notebook should carry, given what it writes."""
    if not published:
        return nb_path.name
    return f"{panel_spec(published)}{SEP}{topic_of(nb_path)}.ipynb"


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _pretty(published: list[str]) -> str:
    """'4K, S12, S13, S14' for human-readable output."""
    out = []
    for p in published:
        if m := _MAIN4_RE.match(p):
            out.append("4" + m.group(1))
        elif m := _MAIN5_RE.match(p):
            out.append("5" + m.group(1))
        elif _MAIN5_ANY_RE.match(p):
            out.append("5")
        elif m := _SUPP_RE.match(p):
            out.append(f"S{int(m.group(1)):02d}")
        else:
            out.append(p)
    return ", ".join(out)


def survey(root: Path = HERE):
    """[(path, published, unmapped, expected_name), ...] for every notebook."""
    panels = load_paper_panels()
    rows = []
    for nb in sorted(root.glob("*.ipynb")):
        published, unmapped = written_panels(nb, panels)
        rows.append((nb, published, unmapped, expected_name(nb, published)))
    return rows


def cmd_check(rows) -> int:
    bad = 0
    for nb, published, unmapped, want in rows:
        if not published:
            keep = nb.name in NO_PANEL_KEEP
            print(f"{'OK  ' if keep else 'NOTE'} {nb.name}")
            print(f"       writes no published panel"
                  + ("" if keep else "; name left alone, add it to NO_PANEL_KEEP if intended"))
            continue
        if nb.name == want:
            print(f"OK   {nb.name}")
            print(f"       panels: {_pretty(published)}")
        else:
            bad += 1
            claimed = nb.stem.split(SEP, 1)[0] if SEP in nb.stem else "(no panel spec)"
            print(f"MISMATCH  {nb.name}")
            print(f"    name claims : {claimed}")
            print(f"    actual write: {_pretty(published)}  ->  {panel_spec(published)}")
            print(f"    expected    : {want}")
        if unmapped:
            print(f"       not in PAPER_PANELS, written verbatim: {', '.join(unmapped)}")
    print()
    print(f"{len(rows)} notebooks, {bad} mismatched")
    return 1 if bad else 0


def cmd_fix(rows) -> int:
    renamed = 0
    for nb, published, _unmapped, want in rows:
        if not published or nb.name == want:
            continue
        dest = nb.with_name(want)
        if dest.exists():
            print(f"SKIP {nb.name}: {want} already exists")
            continue
        r = subprocess.run(["git", "mv", nb.name, want], cwd=nb.parent,
                           capture_output=True, text=True)
        if r.returncode != 0:                      # not tracked by git, or no repo
            nb.rename(dest)
            print(f"renamed (plain) {nb.name}  ->  {want}")
        else:
            print(f"renamed (git mv) {nb.name}  ->  {want}")
        renamed += 1
    print(f"\n{renamed} notebooks renamed")
    return 0


def cmd_table(rows) -> int:
    print("| notebook | main panels | supplementary |")
    print("|---|---|---|")
    for nb, published, _unmapped, want in rows:
        name = want if published else nb.name
        main = [p for p in published if _MAIN4_RE.match(p) or _MAIN5_ANY_RE.match(p)]
        supp = [p for p in published if _SUPP_RE.match(p)]
        print(f"| `{name}` | {_pretty(main) or 'none'} | {_pretty(supp) or 'none'} |")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--fix", action="store_true", help="rename mismatched notebooks")
    g.add_argument("--table", action="store_true", help="print the notebook -> panel table")
    args = ap.parse_args()

    rows = survey()
    if args.table:
        return cmd_table(rows)
    if args.fix:
        return cmd_fix(rows)
    return cmd_check(rows)


if __name__ == "__main__":
    sys.exit(main())
