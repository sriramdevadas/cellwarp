#!/usr/bin/env python3
"""Fig 2C — basal-ganglia three-pair Layer-2 replication.
Reads the deposited per-pair layer2_results JSONs and renders the compression
(pre vs post-rotation Krzanowski S at k=5) with permutation-null markers, under
both weightings. PLOS spec: Arial, 8-12 pt, RGB, <=19.05 cm wide."""
import csv, io, json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.font_manager as fm

# PLOS: Arial 8-12pt. Fall back to a sans if Arial absent.
for cand in ("Arial", "Helvetica", "DejaVu Sans"):
    if any(cand in f.name for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = cand; break
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "savefig.dpi": 600})

BG = pathlib.Path(__file__).resolve().parent / "bg_results"   # vendored, repo-relative

def marker_fraction(pair):
    """n/N cell types whose rank-1 CPC1 driver is a canonical identity marker
    (per-gene-standardized scheme B), derived from the vendored driver CSV."""
    with open(BG/f"layer2_cpc1_drivers_{pair}_W2_schemeB.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    numer = sum(1 for r in rows if r["rank1_class"] == "canonical identity marker")
    return f"{numer}/{len(rows)}"

PAIRS = [(k, lbl, marker_fraction(k)) for (k, lbl) in
         [("Human_Macaque","Human–Macaque"),
          ("Human_Marmoset","Human–Marmoset"),
          ("Macaque_Marmoset","Macaque–Marmoset")]]
SCHEMES = [("W0_unscaled","Unscaled"),("W2_schemeB","Per-gene standardized")]
PRE, POST = "#4a7fb5", "#e08214"   # CVD-validated pair (ΔE 22 protan)

def load(pair, scheme):
    d = json.load(open(BG/f"layer2_results_{pair}.json"))[scheme]["layer2"]["k5"]
    return d["S_pre"], d["S_post"], d["null_mean_pre"], d["null_mean_post"]

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), sharey=True)
x = range(len(PAIRS)); w = 0.34
for ax, (skey, sname) in zip(axes, SCHEMES):
    for i,(pkey,plabel,mk) in enumerate(PAIRS):
        s_pre,s_post,n_pre,n_post = load(pkey, skey)
        b1=ax.bar(i-w/2, s_pre, w, color=PRE, edgecolor="white", linewidth=0.5, zorder=3)
        b2=ax.bar(i+w/2, s_post, w, color=POST, edgecolor="white", linewidth=0.5, zorder=3)
        # permutation-null markers (open diamonds) — every bar clears its null
        ax.plot(i-w/2, n_pre, marker="D", ms=4, mfc="white", mec="#444", mew=0.8, zorder=4)
        ax.plot(i+w/2, n_post, marker="D", ms=4, mfc="white", mec="#444", mew=0.8, zorder=4)
        # direct value labels (relief for contrast; legibility)
        ax.text(i-w/2, s_pre+0.012, f"{s_pre:.2f}", ha="center", va="bottom", fontsize=6.5, color="#222")
        ax.text(i+w/2, s_post+0.012, f"{s_post:.2f}", ha="center", va="bottom", fontsize=6.5, color="#222")
        # compression bracket
        ax.annotate("", xy=(i+w/2, s_post+0.055), xytext=(i-w/2, s_pre+0.055),
                    arrowprops=dict(arrowstyle="->", color="#666", lw=0.7))
    ax.set_xticks(list(x)); ax.set_xticklabels([p[1] for p in PAIRS], fontsize=7.5, rotation=12)
    ax.set_title(sname, fontsize=8.5, pad=6)
    ax.set_ylim(0, 0.92); ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(length=2.5)
    if skey=="W2_schemeB":
        for i,(pkey,plabel,mk) in enumerate(PAIRS):
            ax.text(i, 0.875, mk, ha="center", va="bottom", fontsize=7.0, color=POST, style="italic")
        # figure-level so the two footnotes stack as two centred lines: this one
        # sat on almost the same baseline as the all-pairs note below, offset to
        # its right, and the pair read as one crowded line.
        fig.text(0.5, 0.078, "italic n/N = canonical identity-marker rank-1 CPC1 drivers (per-gene standardized)",
                 ha="center", fontsize=6.0, color=POST, style="italic")
axes[0].set_ylabel("Krzanowski subspace similarity  S (k = 5)", fontsize=8)
# legend + shared note
leg = [Patch(fc=PRE, label="Pre-rotation S"),
       Patch(fc=POST, label="Post-rotation S (centroid-optimal)"),
       plt.Line2D([0],[0], marker="D", mfc="white", mec="#444", ms=4, ls="none", label="Permutation null (mean)")]
axes[1].legend(handles=leg, fontsize=6.6, frameon=False, loc="upper right", handlelength=1.2)
# The note now says what this panel plots. Both axes draw k = 5 and nothing else
# (`load()` reads `["layer2"]["k5"]`), so "at every k tested" asserted a range the
# reader cannot check here; S14 Table is where every k is tabulated, and it covers
# k in {1, 3, 5} for all three pairs under both weightings. "All three primate
# pairs" replaces "All pairs" because the three drawn here are the primate set,
# not every pair the paper reports.
fig.text(0.5, 0.031, "All three primate pairs, both weightings, at k = 5: post-rotation below pre-rotation, "
                     "each above its null (all p < 10$^{-4}$). S14 Table gives every k.",
         ha="center", fontsize=6.6, color="#333")
fig.subplots_adjust(left=0.085, right=0.985, top=0.90, bottom=0.24, wspace=0.08)

out = pathlib.Path(__file__).resolve().parent   # repo-relative (script's own dir)
# Suppress the embedded creation timestamp so re-running this script in one
# environment reproduces its outputs byte-for-byte. Without it the PDF's
# /CreationDate moved on every run and the file showed as modified with no
# change to the drawing. Passed to the PNG writer too: matplotlib accepts the
# key there and writes no tEXt chunk for a None value, so the two writers are
# spelled the same way rather than one carrying an unexplained exception. This
# matches build_main_figures.py's save(), the writer the five composites go
# through; the TIFF below is written from a buffer and carries no timestamp.
fig.savefig(out/"Fig2C_bg_replication.pdf", metadata={"CreationDate": None})
fig.savefig(out/"Fig2C_bg_replication.png", dpi=300, metadata={"CreationDate": None})
# PLOS-spec TIFF (RGB, no alpha, LZW). matplotlib emits RGBA; flatten to RGB.
try:
    from PIL import Image
    _buf = io.BytesIO()
    fig.savefig(_buf, format="tiff", dpi=300, pil_kwargs={"compression":"tiff_lzw"})
    _buf.seek(0)
    Image.open(_buf).convert("RGB").save(
        out/"Fig2C_bg_replication.tiff", format="TIFF",
        compression="tiff_lzw", dpi=(300, 300))
    print("wrote PDF + PNG + TIFF(RGB, LZW)")
except Exception as e:
    print("PDF+PNG written; TIFF skipped:", e)
print("font:", plt.rcParams["font.family"])
