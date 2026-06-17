#!/usr/bin/env python3
"""
Citation Renumbering & Reference List Reordering for Cell Systems manuscript.

Uses the same parsing logic as citation_audit.py. Renumbers in-text
superscript citations and reorders the reference list to match
Cell Systems first-citation-order requirements.

Biology: Cell Systems requires numbered superscript citations in order
of first appearance.

Math: Bijective mapping from old ref numbers to new ref numbers,
applied via position-based replacement (right-to-left per line)
to avoid collision.
"""

import re
import sys
from pathlib import Path
from collections import OrderedDict

# ─── Configuration ───────────────────────────────────────────────────────────

MANUSCRIPT = Path(__file__).resolve().parent.parent / "docs/submission/manuscript_combined.txt"
REPORT_OUT = Path(__file__).resolve().parent.parent / "docs/submission/citation_audit_report.txt"

# Unicode superscript mappings
SUP_DIGITS = {
    '\u2070': '0', '\u00b9': '1', '\u00b2': '2', '\u00b3': '3',
    '\u2074': '4', '\u2075': '5', '\u2076': '6', '\u2077': '7',
    '\u2078': '8', '\u2079': '9',
}
DIGIT_TO_SUPER = {
    '0': '\u2070', '1': '\u00b9', '2': '\u00b2', '3': '\u00b3',
    '4': '\u2074', '5': '\u2075', '6': '\u2076', '7': '\u2077',
    '8': '\u2078', '9': '\u2079',
}
SUP_MINUS = '\u207b'
SUP_CHARS = set(SUP_DIGITS.keys()) | {SUP_MINUS}

# Renumbering map (old → new), derived from the citation audit
OLD_TO_NEW = {
    16: 19, 17: 20, 18: 21, 19: 23,
    20: 17, 21: 16, 22: 18, 23: 22,
    27: 29, 28: 27, 29: 28,
}
# Fill identity for unchanged refs
for i in range(1, 33):
    if i not in OLD_TO_NEW:
        OLD_TO_NEW[i] = i

# Inverse: new position → old position (for reordering reference list)
NEW_TO_OLD = {v: k for k, v in OLD_TO_NEW.items()}

# Sections to process for citation renumbering (everything except REFERENCES)
SECTIONS_TO_RENUMBER = [
    "TITLE PAGE", "SUMMARY", "INTRODUCTION", "RESULTS", "DISCUSSION",
    "ACKNOWLEDGMENTS", "AUTHOR CONTRIBUTIONS", "DECLARATION OF INTERESTS",
    "FIGURE LEGENDS", "STAR METHODS", "SUPPLEMENTAL ITEM LEGENDS",
]

# ─── Parsing functions (reused from citation_audit.py) ───────────────────────


def is_superscript_char(ch):
    return ch in SUP_CHARS


def decode_superscript(s):
    result = []
    for ch in s:
        if ch in SUP_DIGITS:
            result.append(SUP_DIGITS[ch])
        elif ch == SUP_MINUS:
            result.append('-')
        else:
            result.append(ch)
    return ''.join(result)


def num_to_super(n):
    """Convert an integer to Unicode superscript string."""
    return ''.join(DIGIT_TO_SUPER[d] for d in str(n))


def find_superscript_groups(line):
    """
    Find all superscript groups in a line with their character positions.

    Returns list of dicts: {start, end, sup_str, decoded, is_exponent}
    """
    groups = []
    i = 0
    while i < len(line):
        if is_superscript_char(line[i]):
            start = i
            while i < len(line) and (is_superscript_char(line[i]) or line[i] == ','):
                if line[i] == ',':
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

            # Exponent detection: preceded by pure number (not alphanumeric name)
            is_exponent = False
            preceding_char = line[start - 1] if start > 0 else ''
            if preceding_char.isdigit():
                j = start - 1
                while j >= 0 and (line[j].isdigit() or line[j] in '.×,'):
                    j -= 1
                if j >= 0 and line[j].isalpha():
                    is_exponent = False  # Name ending in digits (e.g., L1000)
                else:
                    is_exponent = True   # Pure number (e.g., 10)

            decoded = decode_superscript(sup_str)
            groups.append({
                'start': start,
                'end': end,
                'sup_str': sup_str,
                'decoded': decoded,
                'is_exponent': is_exponent,
            })
        else:
            i += 1
    return groups


def remap_superscript_group(decoded):
    """
    Given a decoded superscript string (e.g., '21', '2-11'),
    remap ref numbers using OLD_TO_NEW and return the new superscript string.

    Returns (new_sup_str, changed) where changed is True if any ref was remapped.
    """
    changed = False

    if '-' in decoded:
        # Range like "2-11"
        parts = decoded.split('-')
        if len(parts) == 2:
            try:
                range_start = int(parts[0])
                range_end = int(parts[1])
            except ValueError:
                return None, False  # Can't parse — leave unchanged

            new_start = OLD_TO_NEW.get(range_start, range_start)
            new_end = OLD_TO_NEW.get(range_end, range_end)
            changed = (new_start != range_start) or (new_end != range_end)

            if changed:
                # Verify the range is still contiguous after remapping
                old_refs = list(range(range_start, range_end + 1))
                new_refs = [OLD_TO_NEW.get(r, r) for r in old_refs]
                new_refs_sorted = sorted(new_refs)
                if (new_refs_sorted == list(range(new_refs_sorted[0],
                                                  new_refs_sorted[-1] + 1))):
                    # Still a contiguous range
                    new_sup = (num_to_super(new_refs_sorted[0]) + SUP_MINUS +
                               num_to_super(new_refs_sorted[-1]))
                else:
                    # No longer contiguous — expand to comma-separated
                    new_sup = ','.join(num_to_super(r) for r in sorted(new_refs))
                return new_sup, changed
            else:
                return None, False  # No change needed
        else:
            return None, False
    elif ',' in decoded:
        # Comma-separated like "1,2,3"
        parts = decoded.split(',')
        new_parts = []
        for p in parts:
            p = p.strip()
            if p:
                try:
                    old_num = int(p)
                    new_num = OLD_TO_NEW.get(old_num, old_num)
                    if new_num != old_num:
                        changed = True
                    new_parts.append(num_to_super(new_num))
                except ValueError:
                    return None, False
        if changed:
            return ','.join(new_parts), True
        return None, False
    else:
        # Single number like "21"
        try:
            old_num = int(decoded)
        except ValueError:
            return None, False
        new_num = OLD_TO_NEW.get(old_num, old_num)
        if new_num != old_num:
            return num_to_super(new_num), True
        return None, False


def find_references_section_boundary(lines):
    """
    Find the line index (0-based) where the REFERENCES section header starts.

    Returns the index of the first ======== divider line of the REFERENCES block.
    """
    divider_re = re.compile(r'^={8,}$')
    for i in range(len(lines)):
        if divider_re.match(lines[i].strip()):
            # Check if next non-blank line is "REFERENCES"
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and lines[j].strip() == "REFERENCES":
                return i
    return None


def parse_reference_entries_raw(lines, ref_section_start):
    """
    Parse reference entries from the REFERENCES section, preserving
    original line formatting.

    Returns:
    - header_lines: lines from section start to first entry (dividers, blanks, "References")
    - entries: list of (list_of_lines) for each reference entry
    """
    # Find the end of the header (dividers + "REFERENCES" + divider + blanks + "References" + blank)
    divider_re = re.compile(r'^={8,}$')

    i = ref_section_start
    header_lines = []

    # First divider
    header_lines.append(lines[i])
    i += 1

    # "REFERENCES" and possible blanks until second divider
    while i < len(lines):
        header_lines.append(lines[i])
        if divider_re.match(lines[i].strip()):
            i += 1
            break
        i += 1

    # Blank lines after second divider
    while i < len(lines) and not lines[i].strip():
        header_lines.append(lines[i])
        i += 1

    # "References" label line
    if i < len(lines) and lines[i].strip().lower() == 'references':
        header_lines.append(lines[i])
        i += 1

    # Blank line after "References"
    while i < len(lines) and not lines[i].strip():
        header_lines.append(lines[i])
        i += 1

    # Now parse entries separated by blank lines
    entries = []
    current_entry = []

    while i < len(lines):
        line = lines[i]
        if line.strip() == '':
            if current_entry:
                entries.append(current_entry)
                current_entry = []
            # Keep blank line as separator (will be re-added during reconstruction)
        else:
            current_entry.append(line)
        i += 1

    if current_entry:
        entries.append(current_entry)

    return header_lines, entries


# ─── Main ────────────────────────────────────────────────────────────────────


def main():
    # Read original manuscript
    text = MANUSCRIPT.read_text(encoding='utf-8')
    lines = text.split('\n')
    print(f"Read {len(lines)} lines from manuscript")

    # ─── Step 1: Renumber in-text citations ──────────────────────────────

    # Find REFERENCES section boundary
    ref_boundary = find_references_section_boundary(lines)
    if ref_boundary is None:
        print("ERROR: Could not find REFERENCES section boundary")
        sys.exit(1)
    print(f"REFERENCES section starts at line {ref_boundary + 1} (0-indexed: {ref_boundary})")

    # Process all lines before REFERENCES for citation renumbering
    modified_lines = list(lines)  # Copy
    citations_modified = 0
    lines_modified = set()

    for line_idx in range(ref_boundary):
        line = modified_lines[line_idx]
        groups = find_superscript_groups(line)

        # Filter to non-exponent citations that need remapping
        replacements = []  # (start, end, new_sup_str)
        for g in groups:
            if g['is_exponent']:
                continue
            new_sup, changed = remap_superscript_group(g['decoded'])
            if changed and new_sup is not None:
                replacements.append((g['start'], g['end'], new_sup))

        if replacements:
            # Sort by start position descending (right to left) to preserve positions
            replacements.sort(key=lambda r: r[0], reverse=True)
            for start, end, new_sup in replacements:
                line = line[:start] + new_sup + line[end:]
                citations_modified += 1
            modified_lines[line_idx] = line
            lines_modified.add(line_idx + 1)  # 1-based for reporting

    print(f"Step 1: {citations_modified} citation occurrences modified across "
          f"{len(lines_modified)} lines")

    # ─── Step 2: Reorder reference list ──────────────────────────────────

    header_lines, old_entries = parse_reference_entries_raw(lines, ref_boundary)
    print(f"Step 2: Parsed {len(old_entries)} reference entries")

    if len(old_entries) != 32:
        print(f"WARNING: Expected 32 reference entries, found {len(old_entries)}")

    # Build new order: for new position N (1-indexed), find old position
    new_entries = []
    for new_pos in range(1, 33):
        old_pos = NEW_TO_OLD[new_pos]
        old_idx = old_pos - 1  # 0-indexed
        if old_idx < len(old_entries):
            new_entries.append(old_entries[old_idx])
        else:
            print(f"ERROR: Old position {old_pos} out of range")
            sys.exit(1)

    # Reconstruct the REFERENCES section
    ref_section_lines = list(header_lines)
    for i, entry in enumerate(new_entries):
        for entry_line in entry:
            ref_section_lines.append(entry_line)
        # Add blank line between entries (not after last)
        if i < len(new_entries) - 1:
            ref_section_lines.append('')

    # ─── Part C: Reassemble and write ────────────────────────────────────

    # Take modified lines up to REFERENCES boundary + reordered REFERENCES section
    final_lines = modified_lines[:ref_boundary] + ref_section_lines

    # Preserve trailing newline if original had one
    final_text = '\n'.join(final_lines)
    if text.endswith('\n') and not final_text.endswith('\n'):
        final_text += '\n'

    MANUSCRIPT.write_text(final_text, encoding='utf-8')
    print(f"Part C: Wrote updated manuscript ({len(final_lines)} lines)")

    # ─── Part D: Verification ────────────────────────────────────────────

    print("\n" + "=" * 60)
    print("VERIFICATION")
    print("=" * 60)

    # Re-read the just-written file and run the audit logic
    verify_text = MANUSCRIPT.read_text(encoding='utf-8')
    verify_lines = verify_text.split('\n')

    # Find sections in the new file
    verify_ref_boundary = find_references_section_boundary(verify_lines)

    # Extract all citations from non-REFERENCES sections
    verify_citations = []
    for line_idx in range(verify_ref_boundary):
        line = verify_lines[line_idx]
        groups = find_superscript_groups(line)
        for g in groups:
            if g['is_exponent']:
                continue
            decoded = g['decoded']
            # Parse refs from decoded
            refs = parse_refs_from_decoded(decoded)
            if refs:
                verify_citations.append({
                    'line_num': line_idx + 1,
                    'refs': refs,
                    'decoded': decoded,
                    'sup_str': g['sup_str'],
                    'start': g['start'],
                    'end': g['end'],
                })

    # Trace first appearances
    verify_first = {}
    verify_order = []
    for cite in verify_citations:
        for ref in cite['refs']:
            if ref not in verify_first:
                verify_first[ref] = cite
                verify_order.append(ref)

    # Check 1: First-citation order
    order_correct = True
    order_issues = []
    for i, ref_num in enumerate(verify_order):
        expected = i + 1
        if ref_num != expected:
            order_correct = False
            order_issues.append(f"Position {expected}: found ref {ref_num}")

    # Check 2: Citation count preserved
    new_count = len(verify_citations)
    count_preserved = (new_count == 158)  # Expected from audit

    # Check 3: All refs cited
    verify_cited = set()
    for c in verify_citations:
        verify_cited.update(c['refs'])
    all_cited = verify_cited == set(range(1, 33))

    # Check 4: No orphans
    verify_header, verify_entries = parse_reference_entries_raw(
        verify_lines, verify_ref_boundary)
    ref_count = len(verify_entries)
    no_orphans = max(verify_cited) <= ref_count if verify_cited else True
    orphan_refs = [r for r in verify_cited if r > ref_count]

    # Check 5: Reference list count
    ref_list_ok = ref_count == 32

    # Check 6: Spot checks
    spot_checks = [
        (16, "Krzanowski"),
        (17, "Tabula Sapiens"),
        (21, "CellHint"),
        (23, "Andrews"),
        (29, "Storey"),
    ]
    spot_results = []
    for new_ref, expected_context in spot_checks:
        # Find first citation of this ref number and check surrounding text
        found = False
        for c in verify_citations:
            if new_ref in c['refs']:
                line = verify_lines[c['line_num'] - 1]
                # Check if expected context word appears near the citation
                start = max(0, c['start'] - 50)
                end = min(len(line), c['end'] + 50)
                neighborhood = line[start:end]
                if expected_context.lower() in neighborhood.lower():
                    spot_results.append((new_ref, expected_context, True,
                                        neighborhood[:60]))
                    found = True
                else:
                    spot_results.append((new_ref, expected_context, False,
                                        neighborhood[:60]))
                    found = True
                break
        if not found:
            spot_results.append((new_ref, expected_context, False,
                                 "NOT FOUND IN BODY"))

    all_spots_pass = all(r[2] for r in spot_results)

    # Print verification
    all_pass = (order_correct and count_preserved and all_cited
                and not orphan_refs and ref_list_ok and all_spots_pass)

    print(f"\n=== RENUMBERING SUMMARY ===")
    print(f"Citations modified: {citations_modified} occurrences across "
          f"{len(lines_modified)} lines")
    print(f"References reordered: 11 of 32")
    print(f"Verification: {'PASS' if all_pass else 'FAIL'}")
    print(f"  - First-citation order: {'PASS' if order_correct else 'FAIL'}")
    if not order_correct:
        for iss in order_issues[:5]:
            print(f"      {iss}")
        if len(order_issues) > 5:
            print(f"      ... and {len(order_issues) - 5} more")
    print(f"  - Citation count preserved: "
          f"{'PASS' if count_preserved else 'FAIL'} "
          f"(old: 158, new: {new_count})")
    print(f"  - All refs cited: {'PASS' if all_cited else 'FAIL'}"
          f" (cited: {sorted(verify_cited)})")
    print(f"  - No orphans: {'PASS' if not orphan_refs else 'FAIL'}"
          f"{' (' + str(orphan_refs) + ')' if orphan_refs else ''}")
    print(f"  - Ref list count: {'PASS' if ref_list_ok else 'FAIL'} ({ref_count})")
    print(f"  - Spot checks: {sum(1 for r in spot_results if r[2])}/5 PASS")
    for ref, ctx, passed, neighborhood in spot_results:
        status = "PASS" if passed else "FAIL"
        print(f"      Ref {ref} (expect '{ctx}'): {status}")
        if not passed:
            print(f"        Found: \"{neighborhood}\"")

    # ─── Append verification to audit report ─────────────────────────────

    verification_text = []
    verification_text.append("\n\n" + "=" * 72)
    verification_text.append("POST-RENUMBERING VERIFICATION (appended by citation_renumber.py)")
    verification_text.append("=" * 72)
    verification_text.append("")
    verification_text.append(f"Citations modified: {citations_modified} occurrences across "
                             f"{len(lines_modified)} lines")
    verification_text.append(f"References reordered: 11 of 32")
    verification_text.append(f"Overall: {'PASS' if all_pass else 'FAIL'}")
    verification_text.append(f"  First-citation order: {'PASS' if order_correct else 'FAIL'}")
    if not order_correct:
        for iss in order_issues:
            verification_text.append(f"    {iss}")
    verification_text.append(f"  Citation count: old=158, new={new_count} "
                             f"({'MATCH' if count_preserved else 'MISMATCH'})")
    verification_text.append(f"  All 32 refs cited: {'YES' if all_cited else 'NO'}")
    verification_text.append(f"  Orphan citations: {orphan_refs if orphan_refs else 'none'}")
    verification_text.append(f"  Reference list entries: {ref_count}")
    verification_text.append(f"  Spot checks:")
    for ref, ctx, passed, neighborhood in spot_results:
        verification_text.append(
            f"    Ref {ref} ('{ctx}'): {'PASS' if passed else 'FAIL'}")

    verification_text.append("")
    verification_text.append("New reference order (new pos → old pos → author):")
    for new_pos in range(1, 33):
        old_pos = NEW_TO_OLD[new_pos]
        changed = "  ← MOVED" if old_pos != new_pos else ""
        if old_pos - 1 < len(old_entries):
            first_line = old_entries[old_pos - 1][0].strip()[:60]
        else:
            first_line = "??"
        verification_text.append(f"  {new_pos:>2} ← old {old_pos:>2}: {first_line}...{changed}")

    # Write verification to audit report
    with open(REPORT_OUT, 'a', encoding='utf-8') as f:
        f.write('\n'.join(verification_text))

    print(f"\nVerification appended to: {REPORT_OUT}")

    if not all_pass:
        print("\n*** VERIFICATION FAILED — review issues above ***")
        sys.exit(1)


def parse_refs_from_decoded(decoded):
    """Parse decoded superscript string into list of ref numbers."""
    refs = []
    parts = decoded.split(',')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            range_parts = part.split('-')
            if len(range_parts) == 2:
                try:
                    start = int(range_parts[0])
                    end = int(range_parts[1])
                    if 1 <= start <= 200 and 1 <= end <= 200 and start <= end:
                        refs.extend(range(start, end + 1))
                except ValueError:
                    pass
        else:
            try:
                num = int(part)
                if 1 <= num <= 200:
                    refs.append(num)
            except ValueError:
                pass
    return refs


if __name__ == '__main__':
    main()
