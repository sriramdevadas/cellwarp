#!/usr/bin/env python3
"""
Task A: Fix Figure S2 panel labels.

Current state: Two raster images with overlay labels "A–B" and "C–D" as ranges,
plus embedded "A" and "B" sub-labels within each raster.

Target: Individual panel labels A, B, C, D, E matching manuscript legend:
  (A) k-sensitivity scatter plot (top-left)
  (B) k-sensitivity line plots (top-right)
  (C) Full vs 10x-only ranking comparison (bottom-left)
  (D) Protocol confound (bottom-right)
  (E) Expanded negative control distributions (new panel from separate file)

Approach:
  1. Extract raster images from existing PDF
  2. Mask embedded sub-labels with Pillow (white rectangles)
  3. Build new PDF with proper layout, including panel E
  4. Add clean A–E labels in Arial Bold 9pt
"""

import fitz
from PIL import Image
import numpy as np
import io
import os
from pathlib import Path

SUPP_DIR = str(Path(__file__).resolve().parent.parent / "figures/supplementary")
OUT_DIR = str(Path(__file__).resolve().parent.parent / "figures/submission/supplementary")

# ── Step 1: Extract and modify raster images ──────────────────────────

src_doc = fitz.open(os.path.join(SUPP_DIR, "figS2_parameter_protocol_sensitivity_polished.pdf"))
src_page = src_doc[0]

# Get image placement info from the page
# Top image: bbox [3.6, 6.0, 385.8, 126.8]  → 1593×504 px
# Bottom image: bbox [3.6, 229.7, 385.8, 376.4] → 1593×612 px
images_info = src_page.get_images(full=True)

modified_rasters = []
for idx, img_info in enumerate(images_info):
    xref = img_info[0]
    pix = fitz.Pixmap(src_doc, xref)
    if pix.n > 4:
        pix = fitz.Pixmap(fitz.csRGB, pix)
    img_bytes = pix.tobytes("png")
    img = Image.open(io.BytesIO(img_bytes))
    arr = np.array(img)
    h, w = arr.shape[:2]
    print(f"Raster {idx}: {w}×{h}")

    # Mask embedded "A" label in upper-left
    # Raster 0: dark at rows 15-43, cols 15-78
    # Raster 1: dark at rows 16-44, cols 15-42
    # Use generous rect to cover the bold letter cleanly
    arr[8:52, 8:55] = 255  # white out left "A" label region

    # Mask embedded "B" label in upper-right half
    # Raster 0: dark at rows 15-59, cols 894-995 (but 894-917 is the B itself)
    # Raster 1: dark at rows 15-43, cols 894-917
    # Be conservative — only mask the "B" letter area, not axis content
    arr[8:52, 886:928] = 255  # white out right "B" label region

    modified = Image.fromarray(arr)
    out_path = f"/tmp/s2_raster_{idx}_clean.png"
    modified.save(out_path, dpi=(300, 300))
    modified_rasters.append((out_path, w, h))
    print(f"  Saved cleaned raster to {out_path}")

src_doc.close()

# ── Step 2: Build new PDF ─────────────────────────────────────────────

# Layout constants
LEFT_MARGIN = 3.6
CONTENT_WIDTH = 382.2  # 3.6 to 385.8
RIGHT_EDGE = LEFT_MARGIN + CONTENT_WIDTH
LABEL_FONT_SIZE = 9.0
LABEL_Y_OFFSET = 3.4   # baseline above image top
ROW_GAP = 6.0           # gap between rows (label area + breathing room)

# Calculate image display dimensions (preserve aspect ratio at content width)
raster0_w, raster0_h = modified_rasters[0][1], modified_rasters[0][2]
raster1_w, raster1_h = modified_rasters[1][1], modified_rasters[1][2]

display_width = CONTENT_WIDTH
r0_display_h = raster0_h * (display_width / raster0_w)  # ~120.9 pts
r1_display_h = raster1_h * (display_width / raster1_w)  # ~146.8 pts

# Panel E: scale negative control figure
neg_doc = fitz.open(os.path.join(SUPP_DIR, "negative_control_distributions.pdf"))
neg_page = neg_doc[0]
neg_w, neg_h = neg_page.rect.width, neg_page.rect.height  # 568.6 × 423.0

# Scale E to content width
e_scale = display_width / neg_w
e_display_w = display_width
e_display_h = neg_h * e_scale

print(f"\nLayout calculations:")
print(f"  Raster 0 display: {display_width:.1f} × {r0_display_h:.1f}")
print(f"  Raster 1 display: {display_width:.1f} × {r1_display_h:.1f}")
print(f"  Panel E display: {e_display_w:.1f} × {e_display_h:.1f} (scale {e_scale:.3f})")

# Y positions
y = 0.0

# Row 1: Labels A, B + raster 0
label_ab_y = y + LABEL_FONT_SIZE - 1   # baseline for label text
img0_y = y + LABEL_FONT_SIZE + LABEL_Y_OFFSET
img0_bottom = img0_y + r0_display_h

# Row 2: Labels C, D + raster 1
label_cd_y = img0_bottom + ROW_GAP + LABEL_FONT_SIZE - 1
img1_y = img0_bottom + ROW_GAP + LABEL_FONT_SIZE + LABEL_Y_OFFSET
img1_bottom = img1_y + r1_display_h

# Row 3: Label E + negative control
label_e_y = img1_bottom + ROW_GAP + LABEL_FONT_SIZE - 1
img_e_y = img1_bottom + ROW_GAP + LABEL_FONT_SIZE + LABEL_Y_OFFSET
img_e_bottom = img_e_y + e_display_h

page_width = RIGHT_EDGE
page_height = img_e_bottom + 4.0  # small bottom margin

print(f"\n  Page size: {page_width:.1f} × {page_height:.1f} pts")
print(f"  Row 1 (A,B): label baseline y={label_ab_y:.1f}, image y={img0_y:.1f}–{img0_bottom:.1f}")
print(f"  Row 2 (C,D): label baseline y={label_cd_y:.1f}, image y={img1_y:.1f}–{img1_bottom:.1f}")
print(f"  Row 3 (E):   label baseline y={label_e_y:.1f}, image y={img_e_y:.1f}–{img_e_bottom:.1f}")

# Create new PDF
new_doc = fitz.open()
new_page = new_doc.new_page(width=page_width, height=page_height)

# Insert raster 0 (panels A, B)
r0_rect = fitz.Rect(LEFT_MARGIN, img0_y, RIGHT_EDGE, img0_bottom)
new_page.insert_image(r0_rect, filename=modified_rasters[0][0])

# Insert raster 1 (panels C, D)
r1_rect = fitz.Rect(LEFT_MARGIN, img1_y, RIGHT_EDGE, img1_bottom)
new_page.insert_image(r1_rect, filename=modified_rasters[1][0])

# Insert panel E (negative control distributions) as vector PDF
e_rect = fitz.Rect(LEFT_MARGIN, img_e_y, LEFT_MARGIN + e_display_w, img_e_bottom)
new_page.show_pdf_page(e_rect, neg_doc, 0)

# ── Step 3: Add panel labels ─────────────────────────────────────────

# Label positions: left labels at x=3.6, right labels at x≈219.6 (right-half start)
# Right panel starts at approximately pixel 880 out of 1593 → PDF x ≈ 3.6 + 880/1593*382.2 ≈ 214.8
# But for consistency with other figures, use x=219.6 (matching S1, S3, S5)
LEFT_LABEL_X = LEFT_MARGIN
RIGHT_LABEL_X = 219.6

# Use Arial Bold — fitz font name
font = fitz.Font("helv")  # Helvetica is the PDF standard; we'll use fontname for insert_text

label_specs = [
    ("A", LEFT_LABEL_X, label_ab_y),
    ("B", RIGHT_LABEL_X, label_ab_y),
    ("C", LEFT_LABEL_X, label_cd_y),
    ("D", RIGHT_LABEL_X, label_cd_y),
    ("E", LEFT_LABEL_X, label_e_y),
]

for label, x, y_baseline in label_specs:
    # insert_text uses the baseline y coordinate
    new_page.insert_text(
        fitz.Point(x, y_baseline),
        label,
        fontname="helv",       # Helvetica (maps to Arial in most renderers)
        fontsize=LABEL_FONT_SIZE,
        color=(0, 0, 0),
    )
    print(f"  Label '{label}' at ({x:.1f}, {y_baseline:.1f})")

# ── Step 4: Make labels bold ──────────────────────────────────────────
# PyMuPDF's insert_text with "helv" gives regular weight. For bold we need "hebo"
# Let's redo with the bold variant.

# Actually, let's rebuild — remove the regular-weight labels and add bold ones.
# Simpler: start fresh page.

new_doc.close()

# Rebuild with bold labels
new_doc = fitz.open()
new_page = new_doc.new_page(width=page_width, height=page_height)

# Insert images
new_page.insert_image(r0_rect, filename=modified_rasters[0][0])
new_page.insert_image(r1_rect, filename=modified_rasters[1][0])
new_page.show_pdf_page(e_rect, neg_doc, 0)

# Add bold labels using "hebo" (Helvetica-Bold)
for label, x, y_baseline in label_specs:
    new_page.insert_text(
        fitz.Point(x, y_baseline),
        label,
        fontname="hebo",       # Helvetica-Bold (≈ Arial Bold)
        fontsize=LABEL_FONT_SIZE,
        color=(0, 0, 0),
    )

# Save
out_path = os.path.join(OUT_DIR, "figS2_parameter_protocol_sensitivity.pdf")
new_doc.save(out_path, garbage=4, deflate=True)
new_doc.close()
neg_doc.close()

# ── Step 5: Verify ───────────────────────────────────────────────────

verify_doc = fitz.open(out_path)
vpage = verify_doc[0]
print(f"\n=== Verification: {out_path} ===")
print(f"  Page size: {vpage.rect.width:.1f} × {vpage.rect.height:.1f} pts")
blocks = vpage.get_text("dict")["blocks"]
for b in blocks:
    if b["type"] == 0:
        for line in b["lines"]:
            for span in line["spans"]:
                t = span["text"].strip()
                if t and len(t) <= 2:
                    print(f"  Label: '{t}' @ ({span['bbox'][0]:.1f}, {span['bbox'][1]:.1f}) "
                          f"font={span['font']} size={span['size']:.1f} bold={'Bold' in span['font'] or 'bold' in span['font'].lower()}")
images = vpage.get_images(full=True)
print(f"  Images: {len(images)}")
print(f"  File size: {os.path.getsize(out_path):,} bytes")

# Check fonts
fonts = set()
for b in blocks:
    if b["type"] == 0:
        for line in b["lines"]:
            for span in line["spans"]:
                fonts.add(span["font"])
print(f"  Fonts used: {fonts}")

verify_doc.close()
print("\nTask A complete.")
