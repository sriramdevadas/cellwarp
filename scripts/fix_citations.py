#!/usr/bin/env python3
"""
CellWarp manuscript: add missing citation superscripts, renumber references
to first-appearance order, fix ref 13 journal name.
"""
import re
import os
from pathlib import Path

BASE = str(Path(__file__).resolve().parent.parent / "docs")

# ---- Superscript helpers ----
S = {'0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'}
R = {v: k for k, v in S.items()}

def to_sup(n):
    return ''.join(S[d] for d in str(n))

def from_sup(s):
    return int(''.join(R[c] for c in s))

# ---- Old-to-new reference mapping ----
# New ordering by first appearance after superscript additions:
#   Intro: 1(Liang),2(Zhong),3(SAMap←22),4(SATURN←3),5(RIMA←4),6(Icebear←5),
#          7(CAMEX←6),8(SpeciesOT←7),9(CAJAL←8),10(Huynh←9),11(Qin←10),
#          12(Ma←11),13(MartinezAbadias←12)
#   Results: 14(Sun←13),15(PanSci←14),16(CellHint←15),17(Andrews←16),
#            18(TabulaSapiens←26),19(Krzanowski←21),20(TabulaMuris←27),
#            21(L1000←23),22(RIRA←17),23(Qu←18)
#   Discussion: 24(Patel←19),25(CELLxGENE←24 in DataAvail)
#   Methods: 26(Cardini←20),27(Eisenberg←25)
OLD_TO_NEW = {
    1:1, 2:2, 3:4, 4:5, 5:6, 6:7, 7:8, 8:9, 9:10, 10:11, 11:12, 12:13,
    13:14, 14:15, 15:16, 16:17, 17:22, 18:23, 19:24, 20:26, 21:19,
    22:3, 23:21, 24:25, 25:27, 26:18, 27:20
}
NEW_TO_OLD = {v: k for k, v in OLD_TO_NEW.items()}

# ---- Files to process ----
CITE_FILES = [
    'introduction.txt', 'results.txt', 'discussion.txt', 'methods.txt',
    'figure_legends.txt', 'supplementary_exploratory.txt',
    'supplementary_table_S3_caption.txt'
]

# ================================================================
# PHASE 1: Add missing superscripts (current ref numbers)
# ================================================================

def add_citations_to_file(text):
    """Add missing superscript citations everywhere a cited work is named."""
    # Longer names first to avoid partial matches
    patterns = [
        ('Tabula Muris Senis', 27),
        ('Tabula Sapiens',     26),
        ('Tabula Mouse',       27),
        ('CZ CELLxGENE Census',24),
        ('Liang-Wagner',        1),
        ('Liang and Wagner',    1),
        ('Eisenberg and Levanon',25),
        ('Patel and Yanai',    19),
        ('Andrews et al.',     16),
        ('Sun et al.',         13),
        ('Qu et al.',          18),
        ('CellHint',           15),
        ('Krzanowski',         21),
        ('SAMap',              22),
        ('PanSci',             14),
        ('RIRA',               17),
        ('L1000',              23),
    ]
    sup_chars = '⁰¹²³⁴⁵⁶⁷⁸⁹'

    for name, ref in patterns:
        sup = to_sup(ref)
        escaped = re.escape(name)
        # Match name NOT immediately followed by any superscript digit
        # Also NOT immediately followed by a hyphen (e.g. RIRA-only)
        pat = re.compile(f'({escaped})(?![{sup_chars}\\-])')

        def _make_repl(nm, sp, txt):
            def repl(m):
                pos = m.start()
                ls = txt.rfind('\n', 0, pos) + 1
                le = txt.find('\n', pos)
                if le == -1:
                    le = len(txt)
                line = txt[ls:le]

                # --- Exclusion rules ---
                # 1) Table rows (lines with ' | ')
                if ' | ' in line:
                    return m.group(0)

                # 2) Section headings: short lines, no sentence-ending punct
                stripped = line.strip()
                if (stripped and len(stripped) < 80
                        and stripped[-1] not in '.),]!?:;'
                        and not stripped.startswith('(')
                        and not stripped[0:1].isdigit()):
                    return m.group(0)

                # 3) Inside "(STAR Methods, …)" cross-references
                paren_s = txt.rfind('(', ls, pos)
                if paren_s != -1:
                    paren_e = txt.find(')', pos)
                    if paren_e != -1 and paren_e <= le:
                        if 'STAR Methods' in txt[paren_s:paren_e + 1]:
                            return m.group(0)

                # 4) Software package table descriptions
                if nm == 'SAMap' and 'cross-species validation' in line:
                    col = pos - ls
                    # Only skip if this is in the software table (samap | 1.0.14)
                    if 'samap' in line.lower() and '|' in line:
                        return m.group(0)

                return m.group(0) + sp
            return repl

        text = pat.sub(_make_repl(name, sup, text), text)

    return text


# Phase 1.5: Remove trailing redundant superscripts that are now
# superseded by the one placed right after the author/tool name.
CLEANUP = [
    # methods: Krzanowski trailing ²¹
    ('methods.txt', '(perfect alignment) ²¹.', '(perfect alignment).'),
    # methods: SAMap pipeline trailing ²²
    ('methods.txt', 'neighborhood graph)²².', 'neighborhood graph).'),
    # methods: CELLxGENE trailing ²⁴ on cancer line
    ('methods.txt', '(version 2025-11-08)²⁴,', '(version 2025-11-08),'),
]

def cleanup_doubled_superscripts(text):
    """Remove doubled superscripts like (RIRA¹⁷)¹⁷ → (RIRA¹⁷)."""
    sup_chars = '⁰¹²³⁴⁵⁶⁷⁸⁹'
    # Pattern: superscript_seq + ')' + same_superscript_seq
    # e.g. ¹⁷)¹⁷  →  ¹⁷)
    pat = re.compile(r'([' + sup_chars + r']+)\)(\1)')
    text = pat.sub(lambda m: m.group(1) + ')', text)
    # Also handle: superscript_seq + ') ' + same_superscript_seq (with space)
    pat2 = re.compile(r'([' + sup_chars + r']+)\) (\1)')
    text = pat2.sub(lambda m: m.group(1) + ')', text)
    return text


# ================================================================
# PHASE 2: Renumber all superscript citations
# ================================================================

def renumber_citations(text):
    """Single-pass replacement of old superscript numbers → new numbers."""
    original = text

    def repl(m):
        s = m.group(0)
        pos = m.start()
        # Don't touch negative exponents (preceded by ⁻)
        if pos > 0 and original[pos - 1] == '⁻':
            return s
        try:
            num = from_sup(s)
        except (ValueError, KeyError):
            return s
        if num in OLD_TO_NEW:
            return to_sup(OLD_TO_NEW[num])
        return s

    return re.sub(r'[⁰¹²³⁴⁵⁶⁷⁸⁹]+', repl, text)


# ================================================================
# PHASE 3: Rewrite references.txt in new order + fix journal name
# ================================================================

def rewrite_references():
    fpath = os.path.join(BASE, 'references.txt')
    with open(fpath, 'r') as f:
        content = f.read()

    # Parse into {old_num: text_after_number}
    refs = {}
    for line in content.strip().split('\n'):
        line = line.strip()
        if not line or line.lower() == 'references':
            continue
        m = re.match(r'^(\d+)\.\s+(.+)$', line)
        if m:
            refs[int(m.group(1))] = m.group(2)

    # TASK 4: Fix ref 13 journal name  Innovation → Innovation (Camb.)
    if 13 in refs:
        refs[13] = refs[13].replace(' Innovation 4,', ' Innovation (Camb.) 4,')

    # Build new list
    lines = ['References\n']
    for new_num in range(1, 28):
        old_num = NEW_TO_OLD[new_num]
        lines.append(f'\n{new_num}. {refs[old_num]}\n')

    with open(fpath, 'w') as f:
        f.write(''.join(lines))
    print('  Rewrote references.txt (new order + journal fix)')


# ================================================================
# MAIN
# ================================================================

def main():
    # --- Phase 1 ---
    print('Phase 1: Adding missing superscripts …')
    for fname in CITE_FILES:
        fpath = os.path.join(BASE, fname)
        with open(fpath, 'r') as f:
            orig = f.read()
        text = add_citations_to_file(orig)
        if text != orig:
            with open(fpath, 'w') as f:
                f.write(text)
            # count additions
            n = sum(1 for a, b in zip(orig, text) if a != b)
            print(f'  {fname}: changed')
        else:
            print(f'  {fname}: no changes')

    # --- Phase 1.5 ---
    print('\nPhase 1.5: Removing redundant trailing superscripts …')
    for fname, old, new in CLEANUP:
        fpath = os.path.join(BASE, fname)
        with open(fpath, 'r') as f:
            text = f.read()
        if old in text:
            text = text.replace(old, new, 1)   # replace first occurrence only
            with open(fpath, 'w') as f:
                f.write(text)
            print(f'  {fname}: cleaned "{old[:50]}…"')
        else:
            print(f'  {fname}: not found — "{old[:50]}…"')

    # --- Phase 2 ---
    print('\nPhase 1.6: Removing doubled superscripts …')
    for fname in CITE_FILES:
        fpath = os.path.join(BASE, fname)
        with open(fpath, 'r') as f:
            orig = f.read()
        text = cleanup_doubled_superscripts(orig)
        if text != orig:
            with open(fpath, 'w') as f:
                f.write(text)
            print(f'  {fname}: cleaned doubles')

    print('\nPhase 2: Renumbering citations …')
    for fname in CITE_FILES:
        fpath = os.path.join(BASE, fname)
        with open(fpath, 'r') as f:
            orig = f.read()
        text = renumber_citations(orig)
        if text != orig:
            with open(fpath, 'w') as f:
                f.write(text)
            print(f'  {fname}: renumbered')
        else:
            print(f'  {fname}: no changes')

    # --- Phase 3 ---
    print('\nPhase 3: Rewriting references.txt …')
    rewrite_references()

    print('\nAll done.')


if __name__ == '__main__':
    main()
