#!/usr/bin/env python3
"""
Citation & Reference Order Audit for Cell Systems manuscript.

Parses Unicode superscript citations, traces first appearances,
verifies sequential numbering, and flags issues.

Biology: Cell Systems requires numbered superscript citations in order
of first appearance. This script audits compliance.

Math: Sequential ordering verification — each ref N must first appear
before any ref N+1 in the manuscript body.
"""

import re
import sys
from pathlib import Path
from collections import OrderedDict

# ─── Configuration ───────────────────────────────────────────────────────────

MANUSCRIPT = Path(__file__).resolve().parent.parent / "docs/submission/manuscript_combined.txt"
REPORT_OUT = Path(__file__).resolve().parent.parent / "docs/submission/citation_audit_report.txt"

# Unicode superscript digit mapping
SUP_DIGITS = {
    '\u2070': '0',  # ⁰
    '\u00b9': '1',  # ¹
    '\u00b2': '2',  # ²
    '\u00b3': '3',  # ³
    '\u2074': '4',  # ⁴
    '\u2075': '5',  # ⁵
    '\u2076': '6',  # ⁶
    '\u2077': '7',  # ⁷
    '\u2078': '8',  # ⁸
    '\u2079': '9',  # ⁹
}
SUP_MINUS = '\u207b'  # ⁻

# All superscript characters (digits + minus)
SUP_CHARS = set(SUP_DIGITS.keys()) | {SUP_MINUS}

# Section processing order (Cell Systems mandated)
SECTION_ORDER = [
    "TITLE PAGE",
    "SUMMARY",
    "INTRODUCTION",
    "RESULTS",
    "DISCUSSION",
    "ACKNOWLEDGMENTS",
    "AUTHOR CONTRIBUTIONS",
    "DECLARATION OF INTERESTS",
    "FIGURE LEGENDS",
    "STAR METHODS",
    "SUPPLEMENTAL ITEM LEGENDS",
]

# ─── Helpers ─────────────────────────────────────────────────────────────────


def decode_superscript(s: str) -> str:
    """Convert a string of Unicode superscript characters to ASCII digits/minus."""
    result = []
    for ch in s:
        if ch in SUP_DIGITS:
            result.append(SUP_DIGITS[ch])
        elif ch == SUP_MINUS:
            result.append('-')
        else:
            result.append(ch)
    return ''.join(result)


def parse_superscript_group(sup_str: str) -> list[int]:
    """
    Parse a superscript group into a list of reference numbers.

    Handles:
    - Single numbers: '12' -> [12]
    - Ranges: '2-11' -> [2,3,4,...,11]
    - Comma-separated: '1,2' -> [1,2]
    - Combinations: '1,3-5' -> [1,3,4,5]

    Returns empty list if unparseable.
    """
    decoded = decode_superscript(sup_str)

    # Split by comma first
    refs = []
    parts = decoded.split(',')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            # Range like "2-11"
            range_parts = part.split('-')
            if len(range_parts) == 2:
                try:
                    start = int(range_parts[0])
                    end = int(range_parts[1])
                    if 1 <= start <= 200 and 1 <= end <= 200 and start <= end:
                        refs.extend(range(start, end + 1))
                    else:
                        return []  # Invalid range
                except ValueError:
                    return []
            else:
                return []  # Multiple dashes
        else:
            try:
                num = int(part)
                if 1 <= num <= 200:
                    refs.append(num)
                else:
                    return []
            except ValueError:
                return []
    return refs


def is_superscript_char(ch: str) -> bool:
    """Check if a character is a Unicode superscript digit or minus."""
    return ch in SUP_CHARS


def extract_context(line: str, start: int, end: int, width: int = 15) -> str:
    """Extract ±width characters around the superscript in the line."""
    ctx_start = max(0, start - width)
    ctx_end = min(len(line), end + width)
    ctx = line[ctx_start:ctx_end]
    return ctx.replace('\n', ' ').replace('\t', ' ')


# ─── Section Parsing ────────────────────────────────────────────────────────


def parse_sections(lines: list[str]) -> dict[str, list[tuple[int, str]]]:
    """
    Parse manuscript into sections based on ======== dividers.

    Returns dict mapping section name -> list of (line_number, line_text).
    Line numbers are 1-based.
    """
    sections = OrderedDict()
    current_section = None
    divider_pattern = re.compile(r'^={8,}$')

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()

        if divider_pattern.match(line):
            # Next non-empty line is the section name
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                candidate = lines[j].strip()
                # Check if next line after candidate is also a divider
                k = j + 1
                while k < len(lines) and not lines[k].strip():
                    k += 1
                if k < len(lines) and divider_pattern.match(lines[k].strip()):
                    current_section = candidate
                    sections[current_section] = []
                    i = k + 1
                    continue
            i += 1
            continue

        if current_section is not None:
            sections[current_section].append((i + 1, lines[i]))  # 1-based line num

        i += 1

    return sections


# ─── Citation Extraction ────────────────────────────────────────────────────


def find_citations_in_line(line: str, line_num: int, section: str) -> list[dict]:
    """
    Find all superscript citation groups in a line.

    Returns list of dicts with keys:
    - line_num: int
    - section: str
    - refs: list[int]
    - superscript_str: str (raw Unicode)
    - context: str
    - is_exponent: bool
    - ambiguous: bool
    - ambiguity_note: str
    """
    results = []

    i = 0
    while i < len(line):
        if is_superscript_char(line[i]):
            # Found start of a superscript group
            start = i
            while i < len(line) and (is_superscript_char(line[i]) or line[i] == ','):
                # Allow commas within superscript groups for comma-separated citations
                # But only if the comma is between superscript chars
                if line[i] == ',':
                    # Check if next non-space char is also superscript
                    next_i = i + 1
                    while next_i < len(line) and line[next_i] == ' ':
                        next_i += 1
                    if next_i < len(line) and is_superscript_char(line[next_i]):
                        i += 1
                        continue
                    else:
                        break
                i += 1
            end = i
            sup_str = line[start:end]

            # Determine if this is an exponent or citation
            # Rule: if preceded by a pure number (only digits, decimal points,
            # ×, etc.), it's a mathematical exponent. But if the preceding
            # token contains letters (e.g., "L1000"), the digits are part of
            # a proper name and the superscript is a citation.
            is_exponent = False
            preceding_char = line[start - 1] if start > 0 else ''

            if preceding_char.isdigit():
                # Look back to see if this is a pure number or an alphanumeric name
                j = start - 1
                while j >= 0 and (line[j].isdigit() or line[j] in '.×,'):
                    j -= 1
                # Check if the character before the numeric part is a letter
                # (indicating a name like "L1000" or "H3K27ac")
                if j >= 0 and line[j].isalpha():
                    is_exponent = False  # It's a name ending in digits → citation
                else:
                    is_exponent = True   # Pure number → exponent

            # Decode the superscript
            decoded = decode_superscript(sup_str)

            # Additional exponent check: purely negative like ⁻⁶ after a digit
            # This is already handled by the preceding digit check

            # Check for ambiguity
            ambiguous = False
            ambiguity_note = ""

            if not is_exponent:
                refs = parse_superscript_group(sup_str)

                # Check for multi-digit ambiguity
                # e.g., ¹² could be ref 12 or refs 1,2
                # Only flag if BOTH interpretations produce valid ref numbers (1-32)
                MAX_REF = 32
                if refs and len(decoded) >= 2 and '-' not in decoded and ',' not in decoded:
                    single_digits = [int(d) for d in decoded if d.isdigit()]
                    multi_digit = int(decoded) if decoded.isdigit() else None
                    if (multi_digit and multi_digit != single_digits[0] and len(single_digits) > 1
                            and all(1 <= d <= MAX_REF for d in single_digits)):
                        ambiguous = True
                        ambiguity_note = (
                            f"'{sup_str}' decoded as '{decoded}' → interpreted as ref {multi_digit}, "
                            f"could also be refs {','.join(str(d) for d in single_digits)}"
                        )

                context = extract_context(line, start, end)
                results.append({
                    'line_num': line_num,
                    'section': section,
                    'refs': refs,
                    'superscript_str': sup_str,
                    'decoded': decoded,
                    'context': context,
                    'is_exponent': False,
                    'ambiguous': ambiguous,
                    'ambiguity_note': ambiguity_note,
                })
            else:
                context = extract_context(line, start, end)
                results.append({
                    'line_num': line_num,
                    'section': section,
                    'refs': [],
                    'superscript_str': sup_str,
                    'decoded': decoded,
                    'context': context,
                    'is_exponent': True,
                    'ambiguous': False,
                    'ambiguity_note': "",
                })
        else:
            i += 1

    return results


# ─── Reference List Parsing ─────────────────────────────────────────────────


def parse_references(ref_lines: list[tuple[int, str]]) -> list[dict]:
    """
    Parse the REFERENCES section into individual reference entries.

    References are separated by blank lines, as sequential paragraphs.
    Returns list of dicts with keys: position, first_author_year, full_text, first_100.
    """
    refs = []
    current_ref = []
    current_start_line = None

    for line_num, line_text in ref_lines:
        stripped = line_text.strip()

        # Skip the "References" header
        if stripped.lower() == 'references':
            continue

        if stripped == '':
            if current_ref:
                full_text = ' '.join(current_ref)
                # Extract first author + year
                first_author_year = extract_first_author_year(full_text)
                refs.append({
                    'position': len(refs) + 1,
                    'first_author_year': first_author_year,
                    'full_text': full_text,
                    'first_100': full_text[:100],
                    'start_line': current_start_line,
                })
                current_ref = []
                current_start_line = None
        else:
            if not current_ref:
                current_start_line = line_num
            current_ref.append(stripped)

    # Don't forget last ref
    if current_ref:
        full_text = ' '.join(current_ref)
        first_author_year = extract_first_author_year(full_text)
        refs.append({
            'position': len(refs) + 1,
            'first_author_year': first_author_year,
            'full_text': full_text,
            'first_100': full_text[:100],
            'start_line': current_start_line,
        })

    return refs


def extract_first_author_year(text: str) -> str:
    """Extract 'LastName et al. (YYYY)' or 'LastName (YYYY)' from a reference."""
    # Try to match patterns like "Author, ... (YYYY)." or "Author et al. (YYYY)."
    # First get the first author last name
    first_author = text.split(',')[0].strip()
    # Remove "The " prefix for consortium refs
    if first_author.startswith("The "):
        first_author = first_author[4:]

    # Find year in parentheses
    year_match = re.search(r'\((\d{4})\)', text)
    if year_match:
        year = year_match.group(1)
    else:
        # Try to find a standalone year
        year_match = re.search(r'(\d{4})', text)
        year = year_match.group(1) if year_match else '????'

    # Check if "et al." is present before the year
    if 'et al.' in text[:text.find(year) if year != '????' else len(text)]:
        return f"{first_author} et al. ({year})"
    else:
        return f"{first_author} ({year})"


# ─── Main Audit ─────────────────────────────────────────────────────────────


def main():
    # Read manuscript
    text = MANUSCRIPT.read_text(encoding='utf-8')
    lines = text.split('\n')

    # Parse sections
    sections = parse_sections(lines)

    print(f"Found {len(sections)} sections: {', '.join(sections.keys())}")

    # ─── Step 1 & 2: Extract all citations section by section ────────────

    all_citations = []  # All citation occurrences (non-exponent)
    all_exponents = []  # All excluded exponents
    all_ambiguous = []  # All ambiguous parses

    for section_name in SECTION_ORDER:
        if section_name not in sections:
            print(f"  WARNING: Section '{section_name}' not found in manuscript")
            continue

        section_lines = sections[section_name]
        for line_num, line_text in section_lines:
            cites = find_citations_in_line(line_text, line_num, section_name)
            for c in cites:
                if c['is_exponent']:
                    all_exponents.append(c)
                else:
                    if c['refs']:  # Only include if we got valid ref numbers
                        all_citations.append(c)
                    if c['ambiguous']:
                        all_ambiguous.append(c)

    # ─── Step 3: Trace first appearance ──────────────────────────────────

    first_appearance = {}  # ref_num -> citation dict (first occurrence)
    appearance_order = []  # ref numbers in order of first appearance

    for cite in all_citations:
        for ref in cite['refs']:
            if ref not in first_appearance:
                first_appearance[ref] = cite
                appearance_order.append(ref)

    # ─── Step 4: Compute new numbering ───────────────────────────────────

    # appearance_order gives us the order refs appear
    old_to_new = {}
    for new_num, old_num in enumerate(appearance_order, 1):
        old_to_new[old_num] = new_num

    # ─── Step 5: Parse reference list ────────────────────────────────────

    ref_section_lines = sections.get("REFERENCES", [])
    references = parse_references(ref_section_lines)

    # ─── Step 7: Flag issues ─────────────────────────────────────────────

    all_cited_refs = set()
    for cite in all_citations:
        all_cited_refs.update(cite['refs'])

    ref_positions = set(r['position'] for r in references)

    # Uncited references
    uncited = [r for r in references if r['position'] not in all_cited_refs]

    # Orphan citations
    orphans = [r for r in all_cited_refs if r not in ref_positions]

    # Order violations
    order_violations = []
    for i in range(len(appearance_order)):
        old_num = appearance_order[i]
        new_num = i + 1
        if old_num != new_num:
            order_violations.append((old_num, new_num))

    # Check if renumbering needed
    renumbering_needed = len(order_violations) > 0

    # Supplemental-only refs: cited only in SUPPLEMENTAL ITEM LEGENDS
    supp_only_refs = set()
    non_supp_refs = set()
    for cite in all_citations:
        if cite['section'] == "SUPPLEMENTAL ITEM LEGENDS":
            supp_only_refs.update(cite['refs'])
        else:
            non_supp_refs.update(cite['refs'])
    supp_only_refs = supp_only_refs - non_supp_refs

    # Duplicate references check
    ref_texts_lower = {}
    duplicates = []
    for r in references:
        key = r['full_text'][:80].lower()
        if key in ref_texts_lower:
            duplicates.append((ref_texts_lower[key], r['position']))
        else:
            ref_texts_lower[key] = r['position']

    # ─── Build Report ────────────────────────────────────────────────────

    report = []
    report.append("=" * 72)
    report.append("CITATION & REFERENCE ORDER AUDIT REPORT")
    report.append(f"Manuscript: {MANUSCRIPT}")
    report.append(f"Generated: scripts/citation_audit.py")
    report.append("=" * 72)
    report.append("")

    # ─── Summary ─────────────────────────────────────────────────────────
    report.append("=== CITATION AUDIT SUMMARY ===")
    report.append(f"Total references in list: {len(references)}")
    report.append(f"Total unique refs cited in body: {len(all_cited_refs)}")
    report.append(f"Total citation occurrences: {len(all_citations)}")
    report.append(f"Renumbering needed: {'YES' if renumbering_needed else 'NO'}")

    if order_violations:
        violations_str = ', '.join(f"ref {old}→should be {new}" for old, new in order_violations[:10])
        if len(order_violations) > 10:
            violations_str += f" ... and {len(order_violations) - 10} more"
        report.append(f"Order violations found: [{violations_str}]")
    else:
        report.append("Order violations found: none")

    report.append(f"Uncited references: {[r['position'] for r in uncited] if uncited else 'none'}")
    report.append(f"Orphan citations: {sorted(orphans) if orphans else 'none'}")
    report.append(f"Ambiguous parses: {len(all_ambiguous)} found" if all_ambiguous else "Ambiguous parses: none")
    report.append(f"Supplemental-only refs: {sorted(supp_only_refs) if supp_only_refs else 'none'}")
    report.append(f"Flagged exponents excluded: {len(all_exponents)}")
    report.append(f"Duplicate references: {duplicates if duplicates else 'none'}")
    report.append("")

    if not renumbering_needed:
        report.append("*** NO RENUMBERING NEEDED — all {} references are in correct "
                       "first-citation order. ***".format(len(references)))
        report.append("")

    # ─── Step 3 Table: First Appearance ──────────────────────────────────
    report.append("")
    report.append("=" * 72)
    report.append("STEP 3: FIRST APPEARANCE OF EACH REFERENCE")
    report.append("=" * 72)
    report.append("")
    report.append(f"{'Old Ref #':<10} | {'First Appears In':<25} | {'Line #':<7} | Context")
    report.append("-" * 10 + "-+-" + "-" * 25 + "-+-" + "-" * 7 + "-+-" + "-" * 40)

    for ref_num in appearance_order:
        cite = first_appearance[ref_num]
        ctx = cite['context'].replace('\n', ' ')[:50]
        report.append(f"{ref_num:<10} | {cite['section']:<25} | {cite['line_num']:<7} | {ctx}")

    # ─── Step 4 Table: Numbering Mapping ─────────────────────────────────
    report.append("")
    report.append("=" * 72)
    report.append("STEP 4: NUMBERING MAPPING (OLD → NEW)")
    report.append("=" * 72)
    report.append("")
    report.append(f"{'Old #':<6} → {'New #':<6} | {'First Author, Year':<40} | Change?")
    report.append("-" * 6 + "---" + "-" * 6 + "-+-" + "-" * 40 + "-+-" + "-" * 7)

    for i, old_num in enumerate(appearance_order):
        new_num = i + 1
        # Find the reference entry
        ref_entry = None
        for r in references:
            if r['position'] == old_num:
                ref_entry = r
                break
        author_year = ref_entry['first_author_year'] if ref_entry else "NOT IN REF LIST"
        changed = "YES" if old_num != new_num else "NO"
        report.append(f"{old_num:<6} → {new_num:<6} | {author_year:<40} | {changed}")

    # Also list any cited refs not in the reference list
    for ref_num in sorted(all_cited_refs - set(appearance_order)):
        report.append(f"{ref_num:<6} → {'??':<6} | {'CITED BUT NOT IN APPEARANCE ORDER':<40} | ERROR")

    # ─── Step 5 Table: Reference List ────────────────────────────────────
    report.append("")
    report.append("=" * 72)
    report.append("STEP 5: REFERENCE LIST ENTRIES")
    report.append("=" * 72)
    report.append("")
    report.append(f"{'Pos':<5} | {'First Author, Year':<40} | {'First 100 chars'}")
    report.append("-" * 5 + "-+-" + "-" * 40 + "-+-" + "-" * 60)

    for r in references:
        cited = "✓" if r['position'] in all_cited_refs else "✗ UNCITED"
        report.append(f"{r['position']:<5} | {r['first_author_year']:<40} | {r['first_100'][:60]}...")
        if r['position'] not in all_cited_refs:
            report.append(f"       *** WARNING: Reference {r['position']} is NEVER CITED in the manuscript body ***")

    report.append("")
    report.append(f"Total references found: {len(references)}")
    report.append(f"Expected: 32")
    report.append(f"Match: {'YES' if len(references) == 32 else 'NO — MISMATCH'}")

    # ─── Step 6 Table: Full Citation Inventory ───────────────────────────
    report.append("")
    report.append("=" * 72)
    report.append("STEP 6: FULL CITATION INVENTORY")
    report.append("=" * 72)
    report.append("")
    report.append(f"{'Line #':<7} | {'Section':<25} | {'Ref(s)':<12} | {'Superscript':<15} | Context")
    report.append("-" * 7 + "-+-" + "-" * 25 + "-+-" + "-" * 12 + "-+-" + "-" * 15 + "-+-" + "-" * 40)

    for cite in all_citations:
        refs_str = ','.join(str(r) for r in cite['refs'])
        if len(cite['refs']) > 3:
            refs_str = f"{cite['refs'][0]}-{cite['refs'][-1]}"
        ctx = cite['context'].replace('\n', ' ')[:40]
        sup = cite['superscript_str'][:15]
        report.append(
            f"{cite['line_num']:<7} | {cite['section']:<25} | [{refs_str}]{'':>{11-len(refs_str)}} | "
            f"{sup:<15} | {ctx}"
        )

    report.append("")
    report.append(f"Total citation occurrences: {len(all_citations)}")

    # ─── Step 7: Issues ──────────────────────────────────────────────────
    report.append("")
    report.append("=" * 72)
    report.append("STEP 7: ISSUE FLAGS")
    report.append("=" * 72)
    report.append("")

    # 1. Uncited references
    report.append("--- 1. Uncited References ---")
    if uncited:
        for r in uncited:
            report.append(f"  Ref {r['position']}: {r['first_author_year']} — NEVER CITED")
    else:
        report.append("  None — all references are cited at least once.")
    report.append("")

    # 2. Orphan citations
    report.append("--- 2. Orphan Citations ---")
    if orphans:
        for o in sorted(orphans):
            report.append(f"  Citation ref {o} — NO MATCHING REFERENCE ENTRY")
    else:
        report.append("  None — all citations have matching reference entries.")
    report.append("")

    # 3. Ambiguous parses
    report.append("--- 3. Ambiguous Parses ---")
    if all_ambiguous:
        # Group by ref number for readability
        ambig_by_ref = {}
        for a in all_ambiguous:
            refs = a['refs']
            key = refs[0] if refs else 0
            if key not in ambig_by_ref:
                ambig_by_ref[key] = []
            ambig_by_ref[key].append(a)

        report.append(f"  NOTE: {len(all_ambiguous)} occurrences of multi-digit Unicode superscripts")
        report.append(f"  are technically ambiguous (e.g., ¹² could be ref 12 or refs 1,2).")
        report.append(f"  All are interpreted as single multi-digit references, which is correct")
        report.append(f"  for a manuscript with 32 references. If separate refs were intended,")
        report.append(f"  commas would be used (e.g., ¹,²). Grouped by reference number below:")
        report.append("")

        for ref_num in sorted(ambig_by_ref.keys()):
            entries = ambig_by_ref[ref_num]
            report.append(f"  Ref {ref_num}: {len(entries)} occurrences — {entries[0]['ambiguity_note']}")
            # Show first 2 contexts
            for e in entries[:2]:
                report.append(f"    e.g., Line {e['line_num']} ({e['section']}): \"{e['context']}\"")
            if len(entries) > 2:
                report.append(f"    ... and {len(entries) - 2} more occurrences")
    else:
        report.append("  None.")
    report.append("")

    # 4. Exponent/citation exclusions
    report.append("--- 4. Exponent Exclusions (digit-preceding rule applied) ---")
    for e in all_exponents:
        report.append(f"  Line {e['line_num']} ({e['section']}): '{e['superscript_str']}' "
                       f"decoded as '{e['decoded']}' — EXCLUDED (preceded by digit)")
        report.append(f"    Context: \"{e['context']}\"")
    report.append(f"  Total excluded: {len(all_exponents)}")
    report.append("")

    # 5. Order violations
    report.append("--- 5. Order Violations ---")
    if order_violations:
        for old, new in order_violations:
            ref_entry = None
            for r in references:
                if r['position'] == old:
                    ref_entry = r
                    break
            author = ref_entry['first_author_year'] if ref_entry else "??"
            first_cite = first_appearance.get(old)
            section = first_cite['section'] if first_cite else "??"
            line = first_cite['line_num'] if first_cite else "??"
            report.append(f"  Ref {old} (={author}) first cited at line {line} ({section})")
            report.append(f"    → Should be renumbered to ref {new}")
    else:
        report.append("  None — all references appear in correct sequential order.")
    report.append("")

    # 6. Duplicate references
    report.append("--- 6. Duplicate References ---")
    if duplicates:
        for d1, d2 in duplicates:
            report.append(f"  Refs {d1} and {d2} appear to be duplicates")
    else:
        report.append("  None.")
    report.append("")

    # 7. Supplemental-only refs
    report.append("--- 7. Supplemental-Only References ---")
    if supp_only_refs:
        for s in sorted(supp_only_refs):
            ref_entry = None
            for r in references:
                if r['position'] == s:
                    ref_entry = r
                    break
            author = ref_entry['first_author_year'] if ref_entry else "??"
            report.append(f"  Ref {s} ({author}) — cited ONLY in SUPPLEMENTAL ITEM LEGENDS")
            report.append(f"    Per Cell Systems rules, consider moving to supplemental numbering [S1], [S2]...")
    else:
        report.append("  None — all cited references appear in at least one non-supplemental section.")
    report.append("")

    # ─── Cross-check: verify each ref 1-32 is cited ─────────────────────
    report.append("--- Cross-check: Reference coverage ---")
    for i in range(1, len(references) + 1):
        if i in all_cited_refs:
            report.append(f"  Ref {i}: ✓ cited")
        else:
            report.append(f"  Ref {i}: ✗ NOT CITED")
    report.append("")

    # ─── Write report ────────────────────────────────────────────────────
    report_text = '\n'.join(report)
    REPORT_OUT.write_text(report_text, encoding='utf-8')

    # ─── Print summary to stdout ─────────────────────────────────────────
    print()
    print("=== CITATION AUDIT SUMMARY ===")
    print(f"Total references in list: {len(references)}")
    print(f"Total unique refs cited in body: {len(all_cited_refs)}")
    print(f"Total citation occurrences: {len(all_citations)}")
    print(f"Renumbering needed: {'YES' if renumbering_needed else 'NO'}")

    if order_violations:
        violations_str = ', '.join(f"ref {old}→{new}" for old, new in order_violations[:10])
        if len(order_violations) > 10:
            violations_str += f" ... +{len(order_violations) - 10} more"
        print(f"Order violations found: [{violations_str}]")
    else:
        print("Order violations found: none")

    print(f"Uncited references: {[r['position'] for r in uncited] if uncited else 'none'}")
    print(f"Orphan citations: {sorted(orphans) if orphans else 'none'}")
    print(f"Ambiguous parses: {len(all_ambiguous)} found" if all_ambiguous else "Ambiguous parses: none")
    print(f"Supplemental-only refs: {sorted(supp_only_refs) if supp_only_refs else 'none'}")
    print(f"Flagged exponents excluded: {len(all_exponents)}")
    print(f"Duplicate references: {duplicates if duplicates else 'none'}")
    print()
    print(f"Full report saved to: {REPORT_OUT}")


if __name__ == '__main__':
    main()
