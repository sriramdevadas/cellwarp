"""
Cell Systems figure styling utilities for CellWarp publication figures.

Applies journal formatting: Arial/Helvetica fonts, 6-8pt text,
colorblind-safe palettes, proper dimensions (85/114/174 mm),
and uppercase panel labels.

Biology: These figures depict cross-species geometric morphometrics
in transcriptomic space — Procrustes analysis of cell type centroids.

Math: Procrustes distance, permutation null distributions,
Spearman correlations, enrichment statistics.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Dimensions (mm → inches)
# ---------------------------------------------------------------------------
MM_PER_INCH = 25.4
COL1_MM = 85        # single column
COL15_MM = 114      # 1.5 column
COL2_MM = 174       # full width (2-column)

COL1 = COL1_MM / MM_PER_INCH
COL15 = COL15_MM / MM_PER_INCH
COL2 = COL2_MM / MM_PER_INCH

DPI = 300

# ---------------------------------------------------------------------------
# Fonts
# ---------------------------------------------------------------------------
FONT_FAMILY = 'Arial'
FONT_SIZE_LABEL = 8       # axis labels
FONT_SIZE_TICK = 7        # tick labels
FONT_SIZE_TITLE = 8       # panel titles
FONT_SIZE_ANNOT = 7       # annotations
FONT_SIZE_PANEL = 10      # panel letter (A, B, C...)
FONT_SIZE_LEGEND = 7      # legend text

# ---------------------------------------------------------------------------
# Colorblind-safe palette (no red-green)
# Using a blue-orange-purple scheme
# ---------------------------------------------------------------------------
# Primary palette
C_BLUE = '#3574B0'       # primary data
C_ORANGE = '#E8820E'     # comparison / highlight
C_PURPLE = '#7B2D8E'     # tertiary
C_TEAL = '#00897B'       # positive controls
C_GRAY = '#757575'       # null / background
C_LIGHTGRAY = '#BDBDBD'  # subtle elements
C_DARKGRAY = '#424242'   # text / borders
C_BLACK = '#212121'

# Lineage category colors
LINEAGE_COLORS = {
    'Immune': '#3574B0',
    'Epithelial': '#E8820E',
    'Stromal': '#7B2D8E',
    'Endothelial': '#00897B',
    'Metabolic': '#C62828',  # dark red (fine with blue/orange, not green)
}

# Sequential palette for bar charts (ranked data)
PALETTE_SEQUENTIAL = plt.cm.viridis

# Discrete palette for datasets
DATASET_COLORS = {
    'Primary': '#3574B0',
    'Sun2023': '#E8820E',
    'PanSci': '#7B2D8E',
    'CellHint': '#00897B',
}

# ---------------------------------------------------------------------------
# Lineage assignment for 35 cell types
# ---------------------------------------------------------------------------
LINEAGE_MAP = {
    'B cell': 'Immune',
    'CD4-positive, alpha-beta T cell': 'Immune',
    'CD8-positive, alpha-beta T cell': 'Immune',
    'T cell': 'Immune',
    'classical monocyte': 'Immune',
    'granulocyte': 'Immune',
    'hematopoietic precursor cell': 'Immune',
    'hematopoietic stem cell': 'Immune',
    'intermediate monocyte': 'Immune',
    'macrophage': 'Immune',
    'mature NK T cell': 'Immune',
    'monocyte': 'Immune',
    'myeloid dendritic cell': 'Immune',
    'myeloid leukocyte': 'Immune',
    'natural killer cell': 'Immune',
    'neutrophil': 'Immune',
    'non-classical monocyte': 'Immune',
    'plasma cell': 'Immune',
    'basal cell': 'Epithelial',
    'bladder urothelial cell': 'Epithelial',
    'enterocyte of epithelium of large intestine': 'Epithelial',
    'epithelial cell': 'Epithelial',
    'large intestine goblet cell': 'Epithelial',
    'luminal epithelial cell of mammary gland': 'Epithelial',
    'pancreatic acinar cell': 'Epithelial',
    'pancreatic ductal cell': 'Epithelial',
    'adventitial cell': 'Stromal',
    'fibroblast': 'Stromal',
    'fibroblast of cardiac tissue': 'Stromal',
    'mesenchymal stem cell': 'Stromal',
    'mesenchymal stem cell of adipose tissue': 'Stromal',
    'smooth muscle cell': 'Stromal',
    'stromal cell': 'Stromal',
    'endothelial cell': 'Endothelial',
    'hepatocyte': 'Metabolic',
}

# Short display names for bar charts
SHORT_NAMES = {
    'CD4-positive, alpha-beta T cell': 'CD4+ T cell',
    'CD8-positive, alpha-beta T cell': 'CD8+ T cell',
    'enterocyte of epithelium of large intestine': 'Enterocyte (LI)',
    'large intestine goblet cell': 'Goblet cell (LI)',
    'luminal epithelial cell of mammary gland': 'Luminal epi. (mammary)',
    'mesenchymal stem cell of adipose tissue': 'MSC (adipose)',
    'mesenchymal stem cell': 'MSC',
    'fibroblast of cardiac tissue': 'Fibroblast (cardiac)',
    'hematopoietic precursor cell': 'HPC',
    'hematopoietic stem cell': 'HSC',
    'mature NK T cell': 'NKT cell',
    'myeloid dendritic cell': 'Myeloid DC',
    'myeloid leukocyte': 'Myeloid leukocyte',
    'non-classical monocyte': 'NC monocyte',
    'intermediate monocyte': 'Int. monocyte',
    'classical monocyte': 'Classical monocyte',
    'natural killer cell': 'NK cell',
    'bladder urothelial cell': 'Urothelial cell',
    'pancreatic acinar cell': 'Acinar cell',
    'pancreatic ductal cell': 'Ductal cell',
    'adventitial cell': 'Adventitial cell',
    'smooth muscle cell': 'Smooth muscle',
    'plasma cell': 'Plasma cell',
    'stromal cell': 'Stromal cell',
    'epithelial cell': 'Epithelial cell',
    'basal cell': 'Basal cell',
    'endothelial cell': 'Endothelial cell',
    'hepatocyte': 'Hepatocyte',
    'macrophage': 'Macrophage',
    'monocyte': 'Monocyte',
    'fibroblast': 'Fibroblast',
    'B cell': 'B cell',
    'T cell': 'T cell',
    'granulocyte': 'Granulocyte',
    'neutrophil': 'Neutrophil',
}


# Optional long-form override applied on top of short_name() for bar/heatmap
# panels where the SHORT_NAMES abbreviations would be ambiguous to a reader
# encountering them without the manuscript text. Keys are the post-short_name
# strings (so callers do EXPAND.get(short_name(ct), short_name(ct))). Used by
# Fig 1E, Fig 3A, and Fig 4A so abbreviations stay consistent across panels.
LABEL_EXPAND = {
    'HPC': 'Hematopoietic precursor cell',
    'HSC': 'Hematopoietic stem cell',
    'MSC': 'Mesenchymal stem cell',
    'MSC (adipose)': 'Mesenchymal stem cell (adipose)',
    'NKT cell': 'NK T cell',
    'NC monocyte': 'Non-classical monocyte',
    'Int. monocyte': 'Intermediate monocyte',
    'Myeloid DC': 'Myeloid dendritic cell',
    'Enterocyte (LI)': 'Enterocyte (large intestine)',
    'Goblet cell (LI)': 'Goblet cell (large intestine)',
    'Luminal epi. (mammary)': 'Luminal epithelial cell (mammary)',
}


def short_name(ct):
    """Return a shortened display name for a cell type."""
    return SHORT_NAMES.get(ct, ct)


def lineage_color(ct):
    """Return the lineage color for a cell type."""
    lin = LINEAGE_MAP.get(ct, 'Immune')
    return LINEAGE_COLORS.get(lin, C_GRAY)


# ---------------------------------------------------------------------------
# Apply global style
# ---------------------------------------------------------------------------
def apply_style():
    """Set matplotlib rcParams for Cell Systems formatting."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': [FONT_FAMILY, 'Helvetica', 'DejaVu Sans'],
        'font.size': FONT_SIZE_TICK,
        'axes.labelsize': FONT_SIZE_LABEL,
        'axes.titlesize': FONT_SIZE_TITLE,
        'xtick.labelsize': FONT_SIZE_TICK,
        'ytick.labelsize': FONT_SIZE_TICK,
        'legend.fontsize': FONT_SIZE_LEGEND,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.minor.size': 1.5,
        'ytick.minor.size': 1.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': DPI,
        'savefig.dpi': DPI,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'pdf.fonttype': 42,     # TrueType (editable text in PDF)
        'ps.fonttype': 42,
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.5,
    })


def add_panel_label(ax, label, x=-0.12, y=1.08):
    """Add uppercase panel label (A, B, C...) at top-left of axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=FONT_SIZE_PANEL, fontweight='bold',
            va='top', ha='left', fontfamily=FONT_FAMILY)


def save_figure(fig, path_stem, tight=True):
    """Save figure as both PDF (vector) and PNG (300 dpi raster).

    Args:
        fig: matplotlib Figure
        path_stem: path without extension (e.g. 'figures/main/fig1')
        tight: use tight_layout before saving
    """
    path = Path(path_stem)
    path.parent.mkdir(parents=True, exist_ok=True)
    if tight:
        try:
            fig.tight_layout()
        except Exception:
            pass
    # pad_inches 0.05 was too tight to survive a rebuild. matplotlib computes the
    # tight bbox from *measured* text extents; when the FreeType a figure is
    # rebuilt under measures a rotated axis label slightly narrower than it draws
    # it, the label's tail falls outside the box and is cropped -- fig3b_pre_post's
    # y-label lost its closing parenthesis exactly this way. A larger pad absorbs
    # the mismatch. Changing this rewrites nothing on its own; it changes what the
    # next rebuild of a panel produces, and no gate rebuilds panels.
    # Suppress the embedded creation timestamp so repeated runs in one environment
    # are byte-identical. This is the third and broadest of the three writers that
    # had the defect -- build_main_figures.py's save() already did it, and
    # build_fig2c_bg.py and make_figure7.py were fixed alongside this -- and it is
    # the one that governs every panel, so a panel PDF used to move on every run
    # with nothing about the drawing changed. Like the pad above, changing it
    # rewrites nothing on its own: it changes what the next rebuild of a panel
    # produces, and no gate rebuilds panels.
    fig.savefig(str(path) + '.pdf', format='pdf', dpi=DPI,
                bbox_inches='tight', pad_inches=0.12,
                metadata={'CreationDate': None})
    fig.savefig(str(path) + '.png', format='png', dpi=DPI,
                bbox_inches='tight', pad_inches=0.12,
                metadata={'CreationDate': None})
    plt.close(fig)
    print(f"  Saved: {path}.pdf and {path}.png")


def format_p(p):
    """Format a p-value for display."""
    if p < 1e-6:
        return 'p < 1e-6'
    elif p < 0.001:
        return f'p = {p:.1e}'
    elif p < 0.01:
        return f'p = {p:.4f}'
    else:
        return f'p = {p:.3f}'


def add_lineage_legend(ax, loc='lower right', ncol=1, title=None):
    """Add standalone lineage color legend to a panel.

    Used in Fig 2A as the canonical legend, referenced by other panels.
    """
    from matplotlib.patches import Patch
    elements = [Patch(facecolor=c, edgecolor='white', label=l)
                for l, c in LINEAGE_COLORS.items()]
    leg = ax.legend(handles=elements, fontsize=FONT_SIZE_LEGEND, frameon=True,
                    loc=loc, ncol=ncol, handlelength=1, handletextpad=0.4,
                    title=title, title_fontsize=FONT_SIZE_LEGEND,
                    edgecolor=C_LIGHTGRAY, fancybox=False)
    leg.get_frame().set_linewidth(0.5)
    return leg


def add_lineage_ref(ax, x=0.03, y=0.03):
    """Add 'Colors as in Fig. 3A' reference note to a panel."""
    ax.text(x, y, 'Colors as in Fig. 3A', transform=ax.transAxes,
            fontsize=5.5, color=C_GRAY, fontstyle='italic',
            ha='left', va='bottom')


def clean_spine(ax):
    """Remove top and right spines, lighten remaining ones."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(C_DARKGRAY)
    ax.spines['bottom'].set_color(C_DARKGRAY)
    ax.tick_params(colors=C_DARKGRAY)


def place_inset_in_clear_corner(ax, text, data_xy=None, **bbox_kwargs):
    """Place a text annotation in the lowest-data-density corner of `ax`.

    Rationale: matplotlib's default upper-right placement collides with
    data density in our figures (null distributions tail right; ranking
    plots accumulate in the upper-right; bar charts peak top-center). We
    pick among the four corners by counting plotted points within a
    normalized bounding box and choose the corner with the fewest points.

    Args:
        ax: matplotlib Axes with data already plotted.
        text: annotation text (string).
        data_xy: optional (xs, ys) tuple. If None, the function inspects
            ax.collections (scatter) and ax.lines (plot) for point data.
        bbox_kwargs: passed to text bbox, e.g. fc='white', ec='black', lw=0.5.

    Returns:
        The matplotlib.text.Text handle created.
    """
    import numpy as _np
    # Gather plotted data points
    if data_xy is not None:
        xs, ys = data_xy
        xs = _np.asarray(xs, dtype=float)
        ys = _np.asarray(ys, dtype=float)
    else:
        xs_list, ys_list = [], []
        for coll in ax.collections:
            off = coll.get_offsets()
            if off.size > 0:
                xs_list.append(off[:, 0]); ys_list.append(off[:, 1])
        for line in ax.lines:
            xd, yd = line.get_xdata(), line.get_ydata()
            if len(xd):
                xs_list.append(_np.asarray(xd, dtype=float))
                ys_list.append(_np.asarray(yd, dtype=float))
        if xs_list:
            xs = _np.concatenate(xs_list); ys = _np.concatenate(ys_list)
        else:
            xs = _np.array([]); ys = _np.array([])

    # Normalize to 0..1 within the current axes view
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # Guard against empty axes
    if xs.size == 0:
        best = 'upper left'
    else:
        xn = (xs - xmin) / max(xmax - xmin, 1e-9)
        yn = (ys - ymin) / max(ymax - ymin, 1e-9)
        # Corners: each is a 0.35-wide band in its corner
        corners = {
            'upper left':  ((0.0, 0.35), (0.65, 1.0)),
            'upper right': ((0.65, 1.0), (0.65, 1.0)),
            'lower left':  ((0.0, 0.35), (0.0, 0.35)),
            'lower right': ((0.65, 1.0), (0.0, 0.35)),
        }
        counts = {}
        for name, ((x0, x1), (y0, y1)) in corners.items():
            mask = (xn >= x0) & (xn <= x1) & (yn >= y0) & (yn <= y1)
            counts[name] = int(mask.sum())
        best = min(counts, key=counts.get)

    # Map to axes-fraction coordinates + alignment
    xy_map = {
        'upper left':  (0.03, 0.97, 'top',    'left'),
        'upper right': (0.97, 0.97, 'top',    'right'),
        'lower left':  (0.03, 0.03, 'bottom', 'left'),
        'lower right': (0.97, 0.03, 'bottom', 'right'),
    }
    fx, fy, va, ha = xy_map[best]
    bbox = dict(boxstyle='round,pad=0.3', fc='white',
                ec=C_DARKGRAY, lw=0.5)
    bbox.update({k: v for k, v in bbox_kwargs.items() if k in
                 ('fc', 'ec', 'lw', 'pad', 'boxstyle', 'alpha')})
    return ax.text(fx, fy, text, transform=ax.transAxes,
                   ha=ha, va=va, fontsize=7, bbox=bbox)
