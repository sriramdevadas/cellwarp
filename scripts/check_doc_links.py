#!/usr/bin/env python3
"""Check that every internal markdown anchor link resolves to a real heading.

Run:  python scripts/check_doc_links.py
Exit: 0 if every anchor resolves, 1 if any does not (and the dead ones are printed).

Why this exists. A dead internal link is a defect nobody notices and every reader hits, and
nothing else in this repository checks documentation at all. It has already caught one:
dispatch 77 found a link introduced by dispatch 76.

CAVEAT 1 -- THIS REIMPLEMENTS GITHUB'S SLUG ALGORITHM AND CAN DRIFT FROM IT.
There is no published specification; the rules below are what GitHub's renderer is observed to
do, and they are the whole of what this file implements:

    1. take the heading text after the leading #'s
    2. strip inline-code backticks and bold/italic markers (`, **, *)
    3. lowercase
    4. delete every character that is not a word character, whitespace or a hyphen
       -- note \\w keeps [A-Za-z0-9_], so UNDERSCORES SURVIVE
    5. collapse whitespace runs to single hyphens

If GitHub changes, or if a heading uses something these rules do not cover (emoji, HTML tags,
duplicate headings -- which GitHub disambiguates with a -1, -2 suffix that this does NOT
implement), this checker can report a live link dead or a dead link live. Treat a disagreement
between this and the rendered page as a bug HERE first.

    THIS HAS BEEN WRONG ONCE, on 2026-08-26 (dispatch 85). Step 4 above also stripped
    underscores, because the same regex that removes emphasis markers (`_`) was applied to the
    whole heading rather than only to the markers. It reported
    `#the-supplementary-tables-run_allsh-does-not-finish-and-the-one-it-now-does` as DEAD when
    the link was correct and the checker was not. The bug had been invisible for 17 anchors
    because none of them had contained an underscore. If you edit the slug rules, add a case to
    SELF_TEST below rather than trusting a green run.

CAVEAT 2 -- SCOPE IS THE TWO READMEs, NOT THE REPOSITORY.
FILES below lists what is checked. Other tracked markdown carries internal links -- CROSSWALK.md,
DEPOSIT_MANIFEST.md, SCOPE.md and the docs/ tree among them -- and none of it is checked here.
That is a deliberate limit, not an oversight: these two files are the reader path, they are the
ones CI executes, and they are the ones that have broken. Widening FILES is safe and cheap if
that stops being true.
"""

from __future__ import annotations

import os
import re
import sys

FILES = ["README.md", "reproduce/README.md"]

_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*$")
_EMPHASIS = re.compile(r"`|\*\*|\*")          # NOT underscore -- see CAVEAT 1
_NON_SLUG = re.compile(r"[^\w\s-]")           # \w keeps A-Za-z0-9_
_WHITESPACE = re.compile(r"\s+")
_LINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")


def slugify(heading_text: str) -> str:
    """Turn heading text into the anchor GitHub is observed to generate."""
    t = _EMPHASIS.sub("", heading_text)
    t = t.lower()
    t = _NON_SLUG.sub("", t)
    return _WHITESPACE.sub("-", t.strip())


def headings(path: str) -> list[str]:
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = _HEADING.match(line)
            if m:
                out.append(slugify(m.group(2)))
    return out


# Cases that pin the rules above. The underscore case is the dispatch-85 bug; it must stay.
SELF_TEST = [
    ("The supplementary tables `run_all.sh` does not finish, and the one it now does",
     "the-supplementary-tables-run_allsh-does-not-finish-and-the-one-it-now-does"),
    ("Install prerequisites first", "install-prerequisites-first"),
    ("Setup (requires Python 3.12)", "setup-requires-python-312"),
    ("**Bold** and *italic* and `code`", "bold-and-italic-and-code"),
]


def self_test() -> bool:
    ok = True
    for text, expected in SELF_TEST:
        got = slugify(text)
        if got != expected:
            print(f"  SELF-TEST FAILED: {text!r}\n    expected {expected!r}\n    got      {got!r}")
            ok = False
    return ok


def main() -> int:
    if not self_test():
        print("slug rules are wrong; fix them before trusting the link results")
        return 1

    cache: dict[str, list[str] | None] = {}
    checked = dead = 0

    for f in FILES:
        if not os.path.exists(f):
            print(f"  MISSING FILE  {f}")
            dead += 1
            continue
        with open(f, encoding="utf-8") as fh:
            text = fh.read()
        for m in _LINK.finditer(text):
            link = m.group(1)
            if "#" not in link:
                continue
            target_rel, anchor = link.split("#", 1)
            # Resolve the anchor against the file the link TARGETS, not the file it sits in.
            target = f if target_rel == "" else os.path.normpath(
                os.path.join(os.path.dirname(f), target_rel))
            checked += 1
            if target not in cache:
                cache[target] = headings(target) if os.path.exists(target) else None
            if cache[target] is None:
                print(f"  DEAD FILE    {f}: {link}")
                dead += 1
            elif anchor not in cache[target]:
                print(f"  DEAD ANCHOR  {f}: #{anchor}")
                print(f"               not a heading in {target}")
                dead += 1

    print(f"  {checked} anchors checked in {len(FILES)} files, {dead} dead")
    return 1 if dead else 0


if __name__ == "__main__":
    sys.exit(main())
