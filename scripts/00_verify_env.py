#!/usr/bin/env python3
"""
CellWarp — Environment Verification Script

Imports every dependency listed in requirements.txt and prints its version.
Run this first to confirm the environment is correctly set up before any
data downloads or analysis.

Usage:
    python scripts/00_verify_env.py

Expected output:
    A table of library names and versions, followed by a PASS/FAIL summary.
    Any library that fails to import will be reported with the error message.
"""

import sys
import warnings
import importlib
from typing import NamedTuple


class Dependency(NamedTuple):
    """A required Python package with its import name and pip name."""
    pip_name: str       # Name used in requirements.txt / pip install
    import_name: str    # Name used in Python import statements
    purpose: str        # Why this project needs it


# All dependencies from requirements.txt, in order.
# import_name differs from pip_name for several packages.
DEPENDENCIES: list[Dependency] = [
    Dependency("numpy",           "numpy",           "Core numerical arrays"),
    Dependency("scipy",           "scipy",           "Scientific computing (Procrustes SVD, stats)"),
    Dependency("scikit-learn",    "sklearn",         "PCA, machine learning utilities"),
    Dependency("pandas",          "pandas",          "Tabular data manipulation"),
    Dependency("matplotlib",      "matplotlib",      "Plotting"),
    Dependency("seaborn",         "seaborn",         "Statistical visualization"),
    Dependency("scanpy",          "scanpy",          "Single-cell RNA-seq analysis"),
    Dependency("anndata",         "anndata",         "Annotated data matrices (.h5ad)"),
    Dependency("cellxgene-census","cellxgene_census", "CZ CELLxGENE data access"),
    Dependency("pybiomart",       "pybiomart",       "Ensembl BioMart ortholog queries"),
    Dependency("gseapy",          "gseapy",          "GO / gene set enrichment analysis"),
    Dependency("samap",           "samap",           "Cross-species atlas mapping (validation only)"),
    Dependency("tqdm",            "tqdm",            "Progress bars for long loops"),
    Dependency("jupyter",         "jupyter",         "Notebook environment"),
]


def get_version(module_name: str) -> str:
    """
    Return the version string for an imported module.

    Tries module.__version__ first (standard), then falls back to
    importlib.metadata for packages that don't expose __version__
    directly (e.g., jupyter).
    """
    try:
        mod = importlib.import_module(module_name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            if hasattr(mod, "__version__"):
                return mod.__version__
    except ImportError:
        pass

    # Fallback: use importlib.metadata (works for installed packages)
    try:
        from importlib.metadata import version
        # Map import name back to distribution name for metadata lookup
        dist_name = module_name.replace("_", "-")
        return version(dist_name)
    except Exception:
        pass

    return "unknown"


def main() -> None:
    print("=" * 65)
    print("CellWarp — Environment Verification")
    print(f"Python {sys.version}")
    print(f"Executable: {sys.executable}")
    print("=" * 65)
    print()

    passed: list[str] = []
    failed: list[tuple[str, str]] = []

    print(f"{'Package':<22} {'Import Name':<20} {'Version':<15} {'Status'}")
    print("-" * 65)

    for dep in DEPENDENCIES:
        try:
            importlib.import_module(dep.import_name)
            ver = get_version(dep.import_name)
            status = "OK"
            passed.append(dep.pip_name)
            print(f"{dep.pip_name:<22} {dep.import_name:<20} {ver:<15} {status}")
        except ImportError as e:
            status = "FAIL"
            failed.append((dep.pip_name, str(e)))
            print(f"{dep.pip_name:<22} {dep.import_name:<20} {'---':<15} {status}")

    print("-" * 65)
    print()

    # Summary
    total = len(DEPENDENCIES)
    n_pass = len(passed)
    n_fail = len(failed)

    if n_fail == 0:
        print(f"RESULT: ALL {total} DEPENDENCIES IMPORTED SUCCESSFULLY")
        print()
        print("Environment is ready. You may proceed to Phase 1 data download.")
    else:
        print(f"RESULT: {n_pass}/{total} passed, {n_fail} FAILED")
        print()
        print("Failed packages:")
        for name, err in failed:
            print(f"  - {name}: {err}")
        print()
        print("To fix, run:")
        print(f"  pip install {' '.join(name for name, _ in failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
