"""Smoke test for the no-download fast-path reproduction wrapper.

Verifies that reproduce/fast_path.py imports, the deposited centroids are
present, and the deterministic headline Procrustes distance reproduces -- WITHOUT
running the full 1,000,000-permutation demo (that ~2-min run is exercised
manually, not in the unit-test suite).
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))


def _load_fast_path():
    fp_path = REPO_ROOT / "reproduce" / "fast_path.py"
    spec = importlib.util.spec_from_file_location("fast_path", fp_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_fast_path_imports_and_inputs_present():
    fp = _load_fast_path()
    assert fp.CENTROIDS.exists(), f"deposited centroids missing: {fp.CENTROIDS}"
    assert fp.PUBLISHED_OBS_NULL == 0.522


def test_fast_path_headline_distance_reproduces():
    fp = _load_fast_path()
    data = np.load(fp.CENTROIDS)
    observed = fp._procrustes_distance(data["human"], data["mouse"])
    # numerical tolerance (NOT exact float equality) -> robust across Mac/Linux BLAS
    assert abs(observed - 61.153) < 1e-2, (
        f"observed Procrustes distance {observed:.4f} differs from 61.153 by > 1e-2"
    )
