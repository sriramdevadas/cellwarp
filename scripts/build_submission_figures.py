#!/usr/bin/env python3
"""
Build submission supplementary figures S1-S3.

Task A: Fix S2 panel labels (A–E individual labels, Arial Bold 9pt).
Task B: Assemble composite S1 and build S3.

Output directory: figures/submission/supplementary/
  figS1_pipeline_validation.pdf          — composite of S1+S7 (panels A–F)
  figS2_parameter_protocol_sensitivity.pdf — fixed labels (panels A–E)
  figS3_bootstrap_rankings.pdf            — 4-panel composite (panels A–B)

This script writes only those three. The renumbered S4/S5 keepers
(figS4_matched_scale_control, figS5_markernull) land in the same output
directory but come from other producers. The former S4/S5 cellhint and
SAMap figures were cut from the submission and are no longer written here.
"""

from pathlib import Path

import fitz
from PIL import Image
import numpy as np
import io
import os
import sys

SUPP_DIR = str(Path(__file__).resolve().parent.parent / "figures/supplementary")
OUT_DIR = str(Path(__file__).resolve().parent.parent / "figures/submission/supplementary")
ARIAL_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

os.makedirs(OUT_DIR, exist_ok=True)


def add_label(page, text, x, y_baseline, size=9.0):
    """Add a label in Arial Bold to a PDF page."""
    font = fitz.Font(fontfile=ARIAL_BOLD)
    tw = fitz.TextWriter(page.rect)
    tw.append(fitz.Point(x, y_baseline), text, font=font, fontsize=size)
    tw.write_text(page)


# ═══════════════════════════════════════════════════════════════════
# TASK A — Figure S2: Fix panel labels
# ═══════════════════════════════════════════════════════════════════

print("=" * 60)
print("TASK A: Figure S2 — Build from individual panel PNGs")
print("=" * 60)

# Build row rasters from individual panel PNGs
PANEL_DIR = str(Path(__file__).resolve().parent.parent / "figures/panels")
panel_pairs = [
    (os.path.join(PANEL_DIR, "figs3a_pca_ranking.png"),
     os.path.join(PANEL_DIR, "figs3b_pca_obsnull.png")),
    (os.path.join(PANEL_DIR, "figs7a_smartseq2.png"),
     os.path.join(PANEL_DIR, "figs7b_smartseq2_fraction.png")),
]

modified_rasters = []
for idx, (left_path, right_path) in enumerate(panel_pairs):
    left = Image.open(left_path)
    right = Image.open(right_path)
    # Scale right to match left height
    rh = int(right.height * left.height / right.height)
    rw = int(right.width * left.height / right.height)
    right = right.resize((rw, rh), Image.LANCZOS)
    # Combine side by side with small gap
    gap = 20
    combined_w = left.width + gap + right.width
    combined = Image.new("RGB", (combined_w, left.height), (255, 255, 255))
    combined.paste(left, (0, 0))
    combined.paste(right, (left.width + gap, 0))
    out_path = f"/tmp/s2_raster_{idx}_clean.png"
    combined.save(out_path, dpi=(300, 300))
    w, h = combined.size
    modified_rasters.append((out_path, w, h))
    print(f"  Row {idx}: {w}×{h}")

# Layout constants
LEFT_MARGIN = 3.6
CONTENT_WIDTH = 382.2
RIGHT_EDGE = LEFT_MARGIN + CONTENT_WIDTH
LABEL_SIZE = 9.0
LABEL_Y_PAD = 3.4
ROW_GAP = 6.0

r0_w, r0_h = modified_rasters[0][1], modified_rasters[0][2]
r1_w, r1_h = modified_rasters[1][1], modified_rasters[1][2]
r0_disp_h = r0_h * (CONTENT_WIDTH / r0_w)
r1_disp_h = r1_h * (CONTENT_WIDTH / r1_w)

neg_doc = fitz.open(os.path.join(SUPP_DIR, "negative_control_distributions.pdf"))
neg_w = neg_doc[0].rect.width
neg_h = neg_doc[0].rect.height
e_scale = CONTENT_WIDTH / neg_w
e_disp_h = neg_h * e_scale

# Y positions
y = 0.0
lab_ab_y = y + LABEL_SIZE - 1
img0_y = y + LABEL_SIZE + LABEL_Y_PAD
img0_end = img0_y + r0_disp_h

lab_cd_y = img0_end + ROW_GAP + LABEL_SIZE - 1
img1_y = img0_end + ROW_GAP + LABEL_SIZE + LABEL_Y_PAD
img1_end = img1_y + r1_disp_h

lab_e_y = img1_end + ROW_GAP + LABEL_SIZE - 1
img_e_y = img1_end + ROW_GAP + LABEL_SIZE + LABEL_Y_PAD
img_e_end = img_e_y + e_disp_h
page_h = img_e_end + 4.0

s2_doc = fitz.open()
s2_page = s2_doc.new_page(width=RIGHT_EDGE, height=page_h)

s2_page.insert_image(fitz.Rect(LEFT_MARGIN, img0_y, RIGHT_EDGE, img0_end),
                     filename=modified_rasters[0][0])
s2_page.insert_image(fitz.Rect(LEFT_MARGIN, img1_y, RIGHT_EDGE, img1_end),
                     filename=modified_rasters[1][0])
s2_page.show_pdf_page(fitz.Rect(LEFT_MARGIN, img_e_y,
                                LEFT_MARGIN + CONTENT_WIDTH, img_e_end),
                      neg_doc, 0)

RIGHT_LABEL_X = LEFT_MARGIN + CONTENT_WIDTH / 2 + 5
for label, x, y_bl in [
    ("A", LEFT_MARGIN, lab_ab_y),
    ("B", RIGHT_LABEL_X, lab_ab_y),
    ("C", LEFT_MARGIN, lab_cd_y),
    ("D", RIGHT_LABEL_X, lab_cd_y),
    ("E", LEFT_MARGIN, lab_e_y),
]:
    add_label(s2_page, label, x, y_bl)

s2_out = os.path.join(OUT_DIR, "figS2_parameter_protocol_sensitivity.pdf")
s2_doc.save(s2_out, garbage=4, deflate=True)
s2_doc.close()
neg_doc.close()
print(f"  Saved: {s2_out} ({os.path.getsize(s2_out):,} bytes)")


# ═══════════════════════════════════════════════════════════════════
# TASK B — Composite S1 and renaming
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("TASK B: Composite S1 + renaming")
print("=" * 60)

# ── S1: Composite from S1_polished (A-B) + S7_polished (C-F) ─────

print("\n--- Building composite S1 ---")

s1_doc = fitz.open(os.path.join(SUPP_DIR, "figS1_independent_pca.pdf"))
s7_doc = fitz.open(os.path.join(SUPP_DIR, "figS7_simulation_study_polished.pdf"))

s1_page = s1_doc[0]
s7_page = s7_doc[0]

s1_w, s1_h = s1_page.rect.width, s1_page.rect.height
s7_w, s7_h = s7_page.rect.width, s7_page.rect.height

print(f"  S1 source: {s1_w:.1f} × {s1_h:.1f} pts")
print(f"  S7 source: {s7_w:.1f} × {s7_h:.1f} pts")

# Step 1: Modify S7 in-place — redact old "X.  Title" labels
# Use search_for to get exact bounding boxes
s7_labels_to_remove = [
    "A.  Detection power",
    "B.  Ranking recovery",
    "C.  Ranking stability",
    "D.  Null calibration",
]

for label_text in s7_labels_to_remove:
    instances = s7_page.search_for(label_text)
    for rect in instances:
        # Expand rect slightly to ensure full coverage
        expanded = fitz.Rect(rect.x0 - 2, rect.y0 - 2, rect.x1 + 2, rect.y1 + 2)
        s7_page.add_redact_annot(expanded, fill=(1, 1, 1))
        print(f"  Redacting '{label_text}' at ({rect.x0:.1f}, {rect.y0:.1f}, {rect.x1:.1f}, {rect.y1:.1f})")

s7_page.apply_redactions()

# Step 2: Add new labels to S7 in its native coordinate space.
# When placed in composite at scale s7_scale, these labels will be scaled down.
# To achieve 9pt in the composite, we need to place them at 9/s7_scale pt in S7.
s7_scale = s1_w / s7_w  # 0.809

# Target composite label positions (must match S1's A, B positions):
#   Left labels: composite x = 3.6  → S7 x = 3.6 / s7_scale = 4.45
#   Right labels: composite x = 219.6 → S7 x = 219.6 / s7_scale = 271.4
# Target font size in composite = 9pt → S7 font size = 9 / s7_scale = 11.1pt
# Target composite y for top row: match S1's label baseline positions
#   S1 label A is at y=14.6, label B at y=2.6
#   For S7 top row, we want composite y ≈ s1_h + gap + 2.6 (like S1's label B)
#   → S7 y = 2.6 / s7_scale = 3.2 (near the S7 page top edge, which is where old labels were)
# For S7 bottom row (panels E, F):
#   Old "C." was at y=222.6, so we use the same y position

s7_font_size = LABEL_SIZE / s7_scale  # ~11.1pt — will render as 9pt after scaling
s7_left_x = LEFT_MARGIN / s7_scale
s7_right_x = RIGHT_LABEL_X / s7_scale

# The S7 labels' original y positions: top row at y=1.9, bottom row at y=222.6
# These are baseline positions. The old labels were "X.  Title" at these baselines.
# Use the same y-baselines for the new single-letter labels.
s7_new_labels = [
    ("C", s7_left_x, 1.9 + LABEL_SIZE),     # top-left
    ("D", s7_right_x, 1.9 + LABEL_SIZE),     # top-right
    ("E", s7_left_x, 222.6 + LABEL_SIZE),    # bottom-left
    ("F", s7_right_x, 222.6 + LABEL_SIZE),   # bottom-right
]

for label, x, y_bl in s7_new_labels:
    add_label(s7_page, label, x, y_bl, size=s7_font_size)
    print(f"  Added S7 label '{label}' at ({x:.1f}, {y_bl:.1f}) size={s7_font_size:.1f}")

# Step 3: Build composite page
gap = 8.0
s7_disp_h = s7_h * s7_scale
composite_h = s1_h + gap + s7_disp_h

print(f"  Composite page: {s1_w:.1f} × {composite_h:.1f} pts")

comp_doc = fitz.open()
comp_page = comp_doc.new_page(width=s1_w, height=composite_h)

# Insert S1 content (panels A, B — labels already correct in Arial-BoldMT)
comp_page.show_pdf_page(fitz.Rect(0, 0, s1_w, s1_h), s1_doc, 0)

# Insert modified S7 content (panels C, D, E, F — with new labels)
s7_rect = fitz.Rect(0, s1_h + gap, s1_w, s1_h + gap + s7_disp_h)
comp_page.show_pdf_page(s7_rect, s7_doc, 0)

comp_out = os.path.join(OUT_DIR, "figS1_pipeline_validation.pdf")
comp_doc.save(comp_out, garbage=4, deflate=True)
comp_doc.close()
s1_doc.close()
s7_doc.close()
print(f"  Saved: {comp_out} ({os.path.getsize(comp_out):,} bytes)")


# ── S3: Direct copy ──────────────────────────────────────────────

print("\n--- S3: build 4-panel composite via composite_figS3.py ---")
import subprocess
composite_script = os.path.join(os.path.dirname(__file__), "composite_figS3.py")
subprocess.run([sys.executable, composite_script], check=True)
s3_out = os.path.join(OUT_DIR, "figS3_bootstrap_rankings.pdf")
print(f"  Saved: {s3_out} ({os.path.getsize(s3_out):,} bytes)")


# ── S4/S5: not produced here ─────────────────────────────────────
# The renumbered keepers figS4_matched_scale_control.{pdf,png} and
# figS5_markernull.pdf also live in OUT_DIR, but other producers write them
# (figS4 from scripts/49_build_figS7_matched_scale.py). The cut cellhint and
# SAMap figures are no longer copied into OUT_DIR at all.


# ═══════════════════════════════════════════════════════════════════
# VERIFICATION
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("VERIFICATION")
print("=" * 60)

expected = {
    "figS1_pipeline_validation.pdf": {"panels": 6, "labels": "ABCDEF"},
    "figS2_parameter_protocol_sensitivity.pdf": {"panels": 5, "labels": "ABCDE"},
    "figS3_bootstrap_rankings.pdf": {"panels": 2, "labels": "AB"},
}

all_ok = True
for fname, spec in expected.items():
    fpath = os.path.join(OUT_DIR, fname)
    if not os.path.exists(fpath):
        print(f"\n  MISSING: {fname}")
        all_ok = False
        continue

    doc = fitz.open(fpath)
    page = doc[0]
    fsize = os.path.getsize(fpath)

    # Find panel labels
    labels = []
    all_fonts = set()
    blocks = page.get_text("dict")["blocks"]
    for b in blocks:
        if b["type"] == 0:
            for line in b["lines"]:
                for span in line["spans"]:
                    all_fonts.add(span["font"])
                    t = span["text"].strip()
                    if (len(t) == 1 and t.isalpha() and t.isupper()
                            and span["size"] >= 8.0
                            and ("Bold" in span["font"] or "bold" in span["font"].lower())):
                        labels.append(t)

    found_labels = "".join(sorted(labels))
    ok = found_labels == spec["labels"]
    if not ok:
        all_ok = False

    print(f"\n  {fname}")
    print(f"    Page: {page.rect.width:.1f} × {page.rect.height:.1f} pts")
    exp_str = spec["labels"]
    print(f"    Labels: [{', '.join(labels)}] — {'OK' if ok else 'MISMATCH (expected ' + exp_str + ')'}")
    print(f"    Fonts: {all_fonts}")
    print(f"    Size: {fsize:,} bytes")

    # Check for embedded fonts
    font_list = page.get_fonts()
    embedded = all(f[3] != "" for f in font_list if f[3] != "")
    print(f"    Font embedding: {len(font_list)} fonts on page")

    doc.close()

# File inventory. This script produces a SUBSET of what lands in OUT_DIR (the
# renumbered S4/S5 keepers come from other producers), so the check is
# presence-of-expected, not a total file count for the shared directory.
print("\n--- Final inventory ---")
files = sorted(os.listdir(OUT_DIR))
print(f"  Total files in {OUT_DIR}: {len(files)}")
for f in files:
    print(f"    {f}  ({os.path.getsize(os.path.join(OUT_DIR, f)):,} bytes)")

missing = [f for f in expected if f not in files]
if missing:
    print(f"\n  MISSING expected figures: {', '.join(missing)}")

if all_ok and not missing:
    print("\n  ALL CHECKS PASSED")
else:
    print("\n  SOME CHECKS FAILED — review above")

print("\nDone.")
