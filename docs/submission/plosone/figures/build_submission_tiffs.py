#!/usr/bin/env python3
"""
Render the five main figures to PLOS-spec TIFFs.

PLOS accepts TIFF or EPS only, at 300 to 600 ppi, RGB (8 bit/channel) or
grayscale, under 10 MB, and requires the file name to match the in-text
citation, so a citation of Fig 1 needs a file named Fig1.tif. The deposited
figures exist as PDF and PNG only.

This rasterises the existing PDFs. It does not regenerate them: Fig 5D is not
bit-reproducible from tracked data, so running the producers could change a
deposited figure. The PDFs are the vector origin and at 300 dpi they land on
the pixel dimensions tabulated in FIGURES, so nothing is resampled. Figs 1 to 4
land there exactly; Fig 5's page height is 1643.798 px at 300 dpi and rounds to
the tabulated 1644, which is the whole of what its aspect tolerance covers.

The source-to-output mapping is an explicit table, never a directory glob, and
every structural assumption is a require() that aborts the build. The written
files are re-opened and checked against both the PDF they came from and the
sibling PNG, which is an independently rasterised copy of the same figure: if
alpha were dropped instead of composited, the mean RGB and the black-pixel
fraction would diverge from the PNG and the build would stop. That independence
survives Fig 5 becoming a copy of make_figure7.py's native vector -- its sibling
PNG is still matplotlib's own Agg rendering of the same figure, not a
rasterisation of the PDF, so the two paths remain genuinely separate.

Usage:  python build_submission_tiffs.py [--output-dir DIR]

Requires pymupdf and pillow. pymupdf is declared in pyproject.toml's [reproduce]
and pinned in [lock]; pillow arrives with matplotlib.
"""

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import pymupdf
    from PIL import Image
except ImportError as exc:  # pragma: no cover - dependencies are declared
    sys.exit("ERROR: pymupdf, pillow and numpy are required: %s" % exc)


# --- Aspect tolerances, per figure ---------------------------------------------
#
# The aspect check exists to catch a TIFF that is not the whole page. Four of the
# five figures can satisfy it exactly and still do; Fig 5 cannot, and the reason is
# arithmetic rather than anything wrong with the figure.
#
# ASPECT_EXACT is the original 1e-6 and applies to Figs 1 to 4. Those come from
# build_main_figures.py, which calls savefig with an explicit figsize and NO
# bbox_inches, so each page is a whole number of 300 dpi pixels in both axes
# (7.4 x 4.16 in -> 2220 x 1248 px) and the rendered aspect reproduces the page
# aspect to float precision.
#
# ASPECT_ROUNDED_PAGE applies to Fig 5 alone. Fig 5 is now the native vector PDF
# from make_figure7.py, which uses bbox_inches="tight": the page is whatever the
# trimmed drawing measures, 457.2 x 394.5117 pt, and 394.5117 pt is 1643.798 px at
# 300 dpi. A TIFF has integer dimensions, so it must round to 1644, and the
# rendered aspect then differs from the page aspect by 1.42e-04. No cropping is
# involved and none could be hidden here: a real crop moves the aspect by percent,
# three orders above this bound.
#
# The same tight bbox is why Fig 5's PNG is pinned at 1902 x 1641 while its TIFF is
# 1905 x 1644. That is a SEPARATE consequence: matplotlib measures the tight bbox
# with each backend's own renderer, so the PDF and Agg backends trim the canvas
# differently, by a few pixels whose sign varies per axis. Both pins are exact.
#
# What this trade gives up, and what still holds. The aspect guard no longer
# distinguishes a whole page from one cropped by less than about 0.02%, which is
# not a failure any producer here can produce. Everything that catches a botched
# render keeps its current strictness: the TIFF and PNG dimension pins are exact
# equality, and MEAN_RGB_TOLERANCE and BLACK_FRACTION_TOLERANCE are unchanged at
# 1.0 and 0.01. Those two are what catch the failure this script was written for --
# alpha dropped instead of composited, which moves the mean by tens of levels and
# the black fraction by tenths. Measured across Fig 5's 3 px size difference they
# move by 0.11 and 0.0014, so they retain their full margin.
#
# To tighten this back to ASPECT_EXACT: drop bbox_inches="tight" from
# make_figure7.py for an explicit figsize landing on whole pixels at 300 dpi, as
# the other four already do. Deferred because make_figure7.py writes a deposited
# figure that is also a packet canonical, and re-tuning its layout risks the
# clipping regression figure_style.save_figure's pad comment records having already
# happened once.
ASPECT_EXACT = 1e-6
ASPECT_ROUNDED_PAGE = 2e-4


# --- Paths and the explicit figure table -------------------------------------

HERE = Path(__file__).resolve().parent

# (output name, source stem, TIFF width px, TIFF height px, sibling PNG width px,
#  sibling PNG height px, aspect tolerance).
# A citation of "Fig N" in the manuscript must resolve to the file "FigN.tif".
#
# The PNG dimensions are pinned SEPARATELY from the TIFF's rather than required to
# equal them. For Figs 1 to 4 they are equal and the pin says so; for Fig 5 they are
# not, and the reason is in the aspect-tolerance block above. Both are exact: the
# guard's value is catching a wrong artifact paired with a figure, and a range would
# give that up.
FIGURES = (
    ("Fig1.tif", "Fig1_configuration_conserved", 2220, 1248, 2220, 1248, ASPECT_EXACT),
    ("Fig2.tif", "Fig2_two_layers_bg", 2220, 1860, 2220, 1860, ASPECT_EXACT),
    ("Fig3.tif", "Fig3_configuration_robust", 2220, 900, 2220, 900, ASPECT_EXACT),
    ("Fig4.tif", "Fig4_pertype_not_resolvable", 2220, 960, 2220, 960, ASPECT_EXACT),
    ("Fig5.tif", "Fig5_conserved_identity_genes", 1905, 1644, 1902, 1641, ASPECT_ROUNDED_PAGE),
)

# Fig2C_bg_replication.* is panel C of Figure 2, not a whole figure. It is left
# alone here and must not be uploaded as Fig 2; Fig2.tif is the full
# Fig2_two_layers_bg composite carrying panels A, B and C.
NOT_A_SUBMISSION_FIGURE = "Fig2C_bg_replication"


# --- PLOS limits and rendering settings --------------------------------------

DPI = 300
PLOS_MIN_WIDTH_PX = 789
PLOS_MAX_WIDTH_PX = 2250
PLOS_MAX_HEIGHT_PX = 2625
PLOS_MAX_BYTES = 10 * 1024 * 1024
# Figures sitting on a bound are likelier to need production intervention, so
# the report names the margin and flags anything this close to one.
NARROW_MARGIN_PX = 50

# The written TIFF is compared with the sibling PNG, rasterised by a different
# engine, so the tolerances allow for antialiasing but not for a botched
# flatten: dropping alpha on a transparent ground moves the mean by tens of
# levels and the black fraction by tenths, both orders above these.
MEAN_RGB_TOLERANCE = 1.0
BLACK_FRACTION_TOLERANCE = 0.01


WHITE = (255, 255, 255, 255)

# TIFF tags read back from the written file.
TAG_BITS_PER_SAMPLE = 258
TAG_COMPRESSION = 259
TAG_SAMPLES_PER_PIXEL = 277
TAG_X_RESOLUTION = 282
TAG_Y_RESOLUTION = 283
TAG_RESOLUTION_UNIT = 296
TAG_EXTRA_SAMPLES = 338
COMPRESSION_LZW = 5
RESOLUTION_UNIT_INCH = 2


# --- Failure ------------------------------------------------------------------

class BuildError(Exception):
    """Raised for any violated assumption. Always fatal."""


def require(condition, message):
    """assert, but not stripped by python -O."""
    if not condition:
        raise BuildError(message)


# --- Rendering ----------------------------------------------------------------

def composite_on_white(image):
    """Flatten an RGBA image onto white and return 8-bit RGB.

    Compositing, not dropping: an unpainted region carries alpha 0 with
    undefined colour underneath, and discarding the channel renders it black.
    """
    require(image.mode == "RGBA",
            "expected an RGBA image to composite, got %s" % image.mode)
    flattened = Image.alpha_composite(Image.new("RGBA", image.size, WHITE), image)
    return flattened.convert("RGB")


def render_pdf(path):
    """Render page 1 at DPI, composited onto white, as 8-bit RGB."""
    require(path.is_file(), "source PDF missing: %s" % path)
    document = pymupdf.open(path)
    try:
        require(document.page_count == 1,
                "%s has %d pages, expected 1" % (path.name, document.page_count))
        page = document[0]
        box = page.rect
        pixmap = page.get_pixmap(dpi=DPI, alpha=True)
        rendered = Image.frombytes("RGBA", (pixmap.width, pixmap.height), pixmap.samples)
        return composite_on_white(rendered), box.width, box.height
    finally:
        document.close()


def load_png_on_white(path):
    """The sibling PNG, composited the same way, as an independent rendering."""
    require(path.is_file(), "source PNG missing: %s" % path)
    with Image.open(path) as png:
        return composite_on_white(png.convert("RGBA"))


def statistics(image):
    """Mean RGB and the fraction of pure-black pixels."""
    pixels = np.asarray(image).astype(np.int64)
    return float(pixels.mean()), float((pixels.sum(axis=2) == 0).mean())


# --- Verification of the file actually written --------------------------------

def verify_tiff(path, width, height, source_aspect, aspect_tolerance, reference):
    """Re-open the written TIFF and check it against PLOS's rules and the PNG."""
    require(path.is_file(), "TIFF was not written: %s" % path)
    size = path.stat().st_size

    with Image.open(path) as tiff:
        require(tiff.format == "TIFF", "%s is %s, not TIFF" % (path.name, tiff.format))
        require(tiff.mode == "RGB",
                "%s is mode %s, expected RGB (PLOS takes RGB or grayscale, and an "
                "alpha channel is not accepted)" % (path.name, tiff.mode))
        require(tiff.size == (width, height),
                "%s is %dx%d, expected %dx%d"
                % (path.name, tiff.size[0], tiff.size[1], width, height))

        tags = tiff.tag_v2
        bits = tags.get(TAG_BITS_PER_SAMPLE)
        require(tuple(bits) == (8, 8, 8),
                "%s has BitsPerSample %s, expected 8 per channel" % (path.name, bits))
        samples = tags.get(TAG_SAMPLES_PER_PIXEL)
        require(samples == 3,
                "%s has SamplesPerPixel %s, expected 3" % (path.name, samples))
        require(TAG_EXTRA_SAMPLES not in tags,
                "%s carries an ExtraSamples tag (%s), so it still declares an alpha "
                "channel" % (path.name, tags.get(TAG_EXTRA_SAMPLES)))
        compression = tags.get(TAG_COMPRESSION)
        require(compression == COMPRESSION_LZW,
                "%s has Compression %s, expected %d (LZW)"
                % (path.name, compression, COMPRESSION_LZW))
        require(tags.get(TAG_RESOLUTION_UNIT) == RESOLUTION_UNIT_INCH,
                "%s has ResolutionUnit %s, expected %d (inch)"
                % (path.name, tags.get(TAG_RESOLUTION_UNIT), RESOLUTION_UNIT_INCH))
        x_resolution = float(tags.get(TAG_X_RESOLUTION))
        y_resolution = float(tags.get(TAG_Y_RESOLUTION))
        require((x_resolution, y_resolution) == (float(DPI), float(DPI)),
                "%s records %sx%s dpi in its tags, expected %dx%d"
                % (path.name, x_resolution, y_resolution, DPI, DPI))
        info_dpi = tuple(float(v) for v in tiff.info.get("dpi", (0, 0)))
        require(info_dpi == (float(DPI), float(DPI)),
                "%s reports dpi %s, expected (%d, %d)" % (path.name, info_dpi, DPI, DPI))
        mean, black = statistics(tiff)

    require(PLOS_MIN_WIDTH_PX <= width <= PLOS_MAX_WIDTH_PX,
            "%s is %d px wide, outside PLOS's %d to %d band"
            % (path.name, width, PLOS_MIN_WIDTH_PX, PLOS_MAX_WIDTH_PX))
    require(height <= PLOS_MAX_HEIGHT_PX,
            "%s is %d px tall, over PLOS's maximum of %d"
            % (path.name, height, PLOS_MAX_HEIGHT_PX))
    require(size < PLOS_MAX_BYTES,
            "%s is %d bytes, over PLOS's %d limit" % (path.name, size, PLOS_MAX_BYTES))

    aspect = width / height
    require(abs(aspect - source_aspect) < aspect_tolerance,
            "%s has aspect %.8f but its source page is %.8f (difference %.2e, "
            "tolerance %.0e), so it is not the whole figure"
            % (path.name, aspect, source_aspect,
               abs(aspect - source_aspect), aspect_tolerance))

    reference_mean, reference_black = statistics(reference)
    require(abs(mean - reference_mean) < MEAN_RGB_TOLERANCE,
            "%s has mean RGB %.4f but the PNG composited on white is %.4f. A "
            "difference this size means the flatten went wrong, most likely alpha "
            "dropped rather than composited." % (path.name, mean, reference_mean))
    require(abs(black - reference_black) < BLACK_FRACTION_TOLERANCE,
            "%s is %.6f pure black but the PNG composited on white is %.6f, which "
            "is the signature of a background that came out black"
            % (path.name, black, reference_black))

    return {
        "name": path.name,
        "bytes": size,
        "width": width,
        "height": height,
        "dpi": info_dpi,
        "mode": "RGB",
        "bits": tuple(bits),
        "compression": "LZW (%d)" % compression,
        "extra_samples": False,
        "aspect": aspect,
        "source_aspect": source_aspect,
        "mean": mean,
        "reference_mean": reference_mean,
        "black": black,
        "reference_black": reference_black,
        "width_margin": PLOS_MAX_WIDTH_PX - width,
        "height_margin": PLOS_MAX_HEIGHT_PX - height,
    }


# --- Build --------------------------------------------------------------------

def build(output_dir):
    require(HERE.is_dir(), "figure directory missing: %s" % HERE)
    for entry in FIGURES:
        stem = entry[1]
        require((HERE / (stem + ".pdf")).is_file(),
                "source PDF missing: %s" % (HERE / (stem + ".pdf")))
        require((HERE / (stem + ".png")).is_file(),
                "source PNG missing: %s" % (HERE / (stem + ".png")))
    names = [entry[0] for entry in FIGURES]
    require(len(set(names)) == len(names), "duplicate output names: %s" % names)

    output_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    for name, stem, width, height, png_width, png_height, aspect_tolerance in FIGURES:
        rendered, box_width, box_height = render_pdf(HERE / (stem + ".pdf"))
        require(rendered.size == (width, height),
                "%s rendered at %dx%d, expected %dx%d at %d dpi"
                % (stem, rendered.size[0], rendered.size[1], width, height, DPI))
        destination = output_dir / name
        rendered.save(destination, format="TIFF", dpi=(DPI, DPI),
                      compression="tiff_lzw")
        reference = load_png_on_white(HERE / (stem + ".png"))
        # Against its own pin, not against the TIFF's. They are equal for four of the
        # five and cannot be for Fig 5; see the aspect-tolerance block. Exact either
        # way, so a wrong PNG paired with a figure still fails here.
        require(reference.size == (png_width, png_height),
                "%s PNG is %s, expected %dx%d"
                % (stem, reference.size, png_width, png_height))
        report = verify_tiff(destination, width, height,
                             box_width / box_height, aspect_tolerance, reference)
        report["source"] = stem
        reports.append(report)
    return reports


def print_report(output_dir, reports):
    print("SOURCE  %s" % HERE)
    print("  renderer                pymupdf %s at %d dpi, alpha composited on white"
          % (pymupdf.__version__, DPI))
    print("  pillow                  %s" % Image.__version__)
    print("  figures                 %d" % len(reports))
    print()
    print("OUTPUT  %s" % output_dir)
    header = ("  %-9s %-6s %-11s %-11s %-5s %-4s %-9s %-8s" %
              ("file", "px w", "px h", "dpi", "mode", "bit", "compress", "bytes"))
    print(header)
    for r in reports:
        print("  %-9s %-6d %-11d %-11s %-5s %-4d %-9s %-8d"
              % (r["name"], r["width"], r["height"],
                 "%gx%g" % r["dpi"], r["mode"], r["bits"][0],
                 r["compression"], r["bytes"]))
    print()
    print("PLOS LIMITS")
    for r in reports:
        print("  %-9s width %d in %d-%d (margin %d)   height %d <= %d (margin %d)   "
              "%.2f MB < 10 MB   alpha=%s"
              % (r["name"], r["width"], PLOS_MIN_WIDTH_PX, PLOS_MAX_WIDTH_PX,
                 r["width_margin"], r["height"], PLOS_MAX_HEIGHT_PX,
                 r["height_margin"], r["bytes"] / 1048576.0, r["extra_samples"]))
    tight = [r["name"] for r in reports
             if r["width_margin"] < NARROW_MARGIN_PX or r["height_margin"] < NARROW_MARGIN_PX]
    if tight:
        print("  NOTE: within %d px of a bound, so likelier to need production "
              "intervention at NAAS: %s" % (NARROW_MARGIN_PX, ", ".join(tight)))
    print()
    print("VISUAL SANITY  TIFF versus the sibling PNG composited on white")
    print("  %-9s %-10s %-10s %-8s %-10s %-10s %-8s"
          % ("file", "mean TIFF", "mean PNG", "delta", "black TIFF", "black PNG", "delta"))
    for r in reports:
        print("  %-9s %-10.4f %-10.4f %-8.4f %-10.6f %-10.6f %-8.6f"
              % (r["name"], r["mean"], r["reference_mean"],
                 abs(r["mean"] - r["reference_mean"]), r["black"], r["reference_black"],
                 abs(r["black"] - r["reference_black"])))
    print()
    print("CONTENT IDENTITY  aspect ratio against the source PDF mediabox")
    for r in reports:
        print("  %-9s %-32s tiff %.8f   pdf %.8f   delta %.2e"
              % (r["name"], r["source"], r["aspect"], r["source_aspect"],
                 abs(r["aspect"] - r["source_aspect"])))
    print()
    print("  %s is panel C of Figure 2 and is not a submission file; it must not "
          "be uploaded as Fig 2." % NOT_A_SUBMISSION_FIGURE)
    print()
    print("OK  %d TIFFs in %s" % (len(reports), output_dir))


# --- Entry point ---------------------------------------------------------------

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Render the five main figures to PLOS-spec TIFFs.")
    parser.add_argument("--output-dir", type=Path, default=HERE,
                        help="where to write FigN.tif (default: %(default)s)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        reports = build(args.output_dir)
    except BuildError as error:
        sys.exit("ERROR: %s" % error)
    print_report(args.output_dir, reports)
    return 0


if __name__ == "__main__":
    sys.exit(main())
