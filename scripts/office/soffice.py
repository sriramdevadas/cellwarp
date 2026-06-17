#!/usr/bin/env python3
"""
Convert .docx to PDF using LibreOffice (soffice --headless).

Usage:
    python soffice.py --headless --convert-to pdf manuscript.docx
    python soffice.py --headless --convert-to pdf /path/to/manuscript.docx

Falls back to alternative methods if LibreOffice is not available.
"""

import subprocess
import sys
import shutil
from pathlib import Path


def find_soffice() -> str | None:
    """Find the soffice binary."""
    # Check PATH first
    soffice = shutil.which("soffice")
    if soffice:
        return soffice
    # Common macOS locations
    mac_paths = [
        "/Applications/LibreOffice.app/Contents/MacOS/soffice",
        "/opt/homebrew/bin/soffice",
    ]
    for p in mac_paths:
        if Path(p).exists():
            return p
    return None


def convert_docx_to_pdf(docx_path: str, output_dir: str | None = None) -> str | None:
    """Convert docx to PDF. Returns output PDF path or None on failure."""
    docx = Path(docx_path).resolve()
    if not docx.exists():
        print(f"ERROR: File not found: {docx_path}")
        return None

    out_dir = Path(output_dir) if output_dir else docx.parent
    pdf_path = out_dir / docx.with_suffix(".pdf").name

    soffice = find_soffice()
    if soffice:
        print(f"Using LibreOffice: {soffice}")
        cmd = [
            soffice, "--headless", "--convert-to", "pdf",
            "--outdir", str(out_dir), str(docx),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and pdf_path.exists():
            print(f"PDF created: {pdf_path}")
            return str(pdf_path)
        else:
            print(f"LibreOffice conversion failed (rc={result.returncode})")
            if result.stderr:
                print(f"  stderr: {result.stderr.strip()}")

    # Fallback: try docx2pdf (pip install docx2pdf)
    try:
        from docx2pdf import convert
        print("Using docx2pdf fallback...")
        convert(str(docx), str(pdf_path))
        if pdf_path.exists():
            print(f"PDF created: {pdf_path}")
            return str(pdf_path)
    except ImportError:
        pass
    except Exception as e:
        print(f"docx2pdf failed: {e}")

    print("ERROR: No PDF converter available.")
    print("  Install LibreOffice: brew install --cask libreoffice")
    print("  Or install docx2pdf: pip install docx2pdf")
    return None


def main():
    # Parse args mimicking soffice CLI: --headless --convert-to pdf file.docx
    args = sys.argv[1:]
    docx_file = None
    output_format = "pdf"

    i = 0
    while i < len(args):
        if args[i] == "--headless":
            i += 1
        elif args[i] == "--convert-to" and i + 1 < len(args):
            output_format = args[i + 1]
            i += 2
        elif args[i] == "--outdir" and i + 1 < len(args):
            i += 2  # handled by convert_docx_to_pdf
        else:
            docx_file = args[i]
            i += 1

    if not docx_file:
        print("Usage: python soffice.py --headless --convert-to pdf <file.docx>")
        sys.exit(1)

    # Resolve relative to docs/submission/ if not absolute
    p = Path(docx_file)
    if not p.is_absolute() and not p.exists():
        alt = Path(__file__).resolve().parent.parent.parent / "docs" / "submission" / p.name
        if alt.exists():
            p = alt
    docx_file = str(p)

    result = convert_docx_to_pdf(docx_file)
    if result:
        # Report file size
        size_mb = Path(result).stat().st_size / (1024 * 1024)
        print(f"PDF size: {size_mb:.2f} MB")
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
