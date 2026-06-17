"""
Central configuration for CellWarp paper reproduction.
All paths relative to repository root.
"""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
OUTPUT_DIR = REPO_ROOT / "output"
FIGURES_DIR = REPO_ROOT / "figures"

RANDOM_SEED = 42
N_COMPONENTS = 33
N_PERMUTATIONS = 10_000
N_PERMUTATIONS_1M = 1_000_000
MAX_CELLS_PER_TYPE = 2_000
