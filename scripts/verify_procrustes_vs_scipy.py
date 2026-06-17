#!/usr/bin/env python3
"""
Verify CellWarp custom SVD Procrustes against scipy.spatial.procrustes.

Loads the primary human-mouse centroid matrices (35 types × 33 PCs) and
compares outputs from both implementations. Reports max absolute differences
for rotation matrix, aligned coordinates, and Procrustes distance.

Key difference between the two implementations:
  - scipy.spatial.procrustes normalizes both point sets to unit Frobenius
    norm BEFORE finding the rotation. It returns "disparity" = sum of squared
    differences between the two unit-norm-standardized, aligned matrices.
  - CellWarp's procrustes_align finds optimal rotation AND scaling without
    pre-normalization. It returns the Frobenius distance in the original
    (centered) coordinate system.

Because of this normalization difference, direct comparison of aligned
coordinates and distance requires rescaling. However, the ROTATION MATRIX
must be identical (SVD of a scalar multiple of a matrix yields the same
U, V matrices).

To make a fair comparison, we also re-implement scipy's normalization
on our side and verify the disparity matches.
"""

import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.spatial import procrustes as scipy_procrustes

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from cellwarp.procrustes import procrustes_align


def main():
    print("=" * 72)
    print("VERIFICATION: CellWarp Procrustes vs scipy.spatial.procrustes")
    print("=" * 72)

    # ------------------------------------------------------------------
    # 1. Load centroid matrices
    # ------------------------------------------------------------------
    data_path = Path(__file__).resolve().parent.parent / "output" / "phase2" / "scaled_35types" / "pca_centroids_35.npz"
    data = np.load(data_path, allow_pickle=True)
    human_pca = data["human"].astype(np.float64)  # (35, 33)
    mouse_pca = data["mouse"].astype(np.float64)  # (35, 33)
    cell_types = data["cell_types"]

    print(f"\nInput data: {data_path.name}")
    print(f"  Human centroids: {human_pca.shape}")
    print(f"  Mouse centroids: {mouse_pca.shape}")
    print(f"  Cell types: {len(cell_types)}")
    print(f"  dtype: {human_pca.dtype}")

    # ------------------------------------------------------------------
    # 2a. Run CellWarp Procrustes
    # ------------------------------------------------------------------
    print("\n" + "-" * 72)
    print("2a. CellWarp custom Procrustes (mouse → human)")
    print("-" * 72)
    cw_result = procrustes_align(human_pca, mouse_pca, allow_reflection=False)

    # ------------------------------------------------------------------
    # 2b. Run scipy.spatial.procrustes
    # ------------------------------------------------------------------
    print("\n" + "-" * 72)
    print("2b. scipy.spatial.procrustes")
    print("-" * 72)
    # scipy convention: procrustes(data1, data2) aligns data2 onto data1
    sp_std_human, sp_std_mouse_aligned, sp_disparity = scipy_procrustes(
        human_pca, mouse_pca
    )
    print(f"  scipy disparity (sum of squared differences): {sp_disparity:.15e}")

    # ------------------------------------------------------------------
    # 3. Extract scipy's internal rotation for comparison
    # ------------------------------------------------------------------
    # scipy internally does:
    #   1. Center both
    #   2. Normalize to unit Frobenius norm
    #   3. SVD of M = X_norm^T @ Y_norm
    #   4. R = V @ U^T (with reflection correction)
    #   5. disparity = sum((X_norm - Y_norm @ R)^2)
    #
    # We replicate this to extract the rotation matrix.

    print("\n" + "-" * 72)
    print("3. Extracting scipy rotation matrix (replicating internal steps)")
    print("-" * 72)

    # Center
    mu_h = human_pca.mean(axis=0)
    mu_m = mouse_pca.mean(axis=0)
    X_c = human_pca - mu_h
    Y_c = mouse_pca - mu_m

    # Normalize to unit Frobenius norm (scipy's standardization)
    norm_X = np.sqrt(np.sum(X_c ** 2))
    norm_Y = np.sqrt(np.sum(Y_c ** 2))
    X_norm = X_c / norm_X
    Y_norm = Y_c / norm_Y

    # SVD of cross-covariance in normalized space
    M_norm = X_norm.T @ Y_norm
    U_n, sigma_n, Vt_n = np.linalg.svd(M_norm)
    V_n = Vt_n.T

    # Reflection correction (same as scipy)
    d_n = np.linalg.det(V_n @ U_n.T)
    D_n = np.eye(X_c.shape[1])
    D_n[-1, -1] = np.sign(d_n)
    R_scipy = V_n @ D_n @ U_n.T

    # Verify: apply scipy's rotation to normalized mouse
    Y_norm_aligned = Y_norm @ R_scipy
    disparity_check = np.sum((X_norm - Y_norm_aligned) ** 2)
    print(f"  Replicated disparity: {disparity_check:.15e}")
    print(f"  scipy disparity:      {sp_disparity:.15e}")
    print(f"  Difference:           {abs(disparity_check - sp_disparity):.2e}")

    # ------------------------------------------------------------------
    # 4. Extract CellWarp's rotation matrix
    # ------------------------------------------------------------------
    # CellWarp uses SVD of X_c^T @ Y_c (NOT normalized)
    # Since X_c^T @ Y_c = (norm_X * norm_Y) * (X_norm^T @ Y_norm),
    # the SVD gives the SAME U, V (rotation), just scaled singular values.
    R_cellwarp = cw_result.rotation

    # ------------------------------------------------------------------
    # 5. Compare outputs
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("COMPARISON RESULTS")
    print("=" * 72)

    TOL = 1e-10

    # --- Rotation matrix ---
    diff_R = np.abs(R_cellwarp - R_scipy)
    max_diff_R = np.max(diff_R)
    frob_diff_R = np.sqrt(np.sum(diff_R ** 2))
    pass_R = max_diff_R < TOL

    print(f"\n  ROTATION MATRIX (33×33):")
    print(f"    Max absolute difference:  {max_diff_R:.2e}")
    print(f"    Frobenius norm of diff:   {frob_diff_R:.2e}")
    print(f"    CellWarp det(R):          {np.linalg.det(R_cellwarp):+.15f}")
    print(f"    scipy   det(R):           {np.linalg.det(R_scipy):+.15f}")
    print(f"    PASS (tol={TOL:.0e}):      {'YES ✓' if pass_R else 'NO ✗'}")

    # --- Aligned coordinates (after accounting for normalization) ---
    # CellWarp: Y_aligned = s * Y_c @ R  (in centered, original-scale coords)
    # scipy:    Y_aligned = Y_norm @ R    (in centered, unit-norm coords)
    #
    # To compare, rescale CellWarp's output to scipy's normalization:
    #   CellWarp_in_scipy_space = cw_result.aligned_target / norm_X
    #
    # Derivation:
    #   cw_aligned = s * Y_c @ R = s * (norm_Y * Y_norm) @ R = s * norm_Y * (Y_norm @ R)
    #   scipy_aligned = Y_norm @ R
    #   ratio = s * norm_Y
    #   So: cw_aligned / (s * norm_Y) = Y_norm @ R = scipy_aligned
    #
    # But scipy returns sp_std_mouse_aligned which is X_norm-space, and
    # sp_std_human = X_c / norm_X. So scipy's aligned mouse is in X_norm space.
    #
    # CellWarp's aligned_target is in X_c space (centered, original scale).
    # To compare: cw_result.aligned_target / norm_X should equal sp_std_mouse_aligned.

    cw_aligned_normalized = cw_result.aligned_target / norm_X
    diff_aligned = np.abs(cw_aligned_normalized - sp_std_mouse_aligned)
    max_diff_aligned = np.max(diff_aligned)
    pass_aligned = max_diff_aligned < TOL

    print(f"\n  ALIGNED COORDINATES (35×33, after normalization to scipy space):")
    print(f"    Max absolute difference:  {max_diff_aligned:.2e}")
    print(f"    PASS (tol={TOL:.0e}):      {'YES ✓' if pass_aligned else 'NO ✗'}")

    # Also verify the reference (human) coordinates in scipy space
    cw_ref_normalized = cw_result.centered_reference / norm_X
    diff_ref = np.abs(cw_ref_normalized - sp_std_human)
    max_diff_ref = np.max(diff_ref)
    pass_ref = max_diff_ref < TOL

    print(f"\n  REFERENCE COORDINATES (35×33, after normalization to scipy space):")
    print(f"    Max absolute difference:  {max_diff_ref:.2e}")
    print(f"    PASS (tol={TOL:.0e}):      {'YES ✓' if pass_ref else 'NO ✗'}")

    # --- Procrustes distance / disparity ---
    # scipy disparity = ||X_norm - Y_norm @ R||^2
    # CellWarp d^2 = ||X_c - s * Y_c @ R||^2
    #
    # Relationship: scipy_disparity = CellWarp_d^2 / norm_X^2
    # (since X_norm = X_c / norm_X and Y_norm_aligned = cw_aligned / norm_X)

    cw_disparity_in_scipy_space = cw_result.distance_squared / (norm_X ** 2)
    diff_disparity = abs(cw_disparity_in_scipy_space - sp_disparity)
    pass_disparity = diff_disparity < TOL

    print(f"\n  PROCRUSTES DISPARITY:")
    print(f"    scipy disparity:              {sp_disparity:.15e}")
    print(f"    CellWarp d²/‖X_c‖²_F:        {cw_disparity_in_scipy_space:.15e}")
    print(f"    Absolute difference:          {diff_disparity:.2e}")
    print(f"    PASS (tol={TOL:.0e}):          {'YES ✓' if pass_disparity else 'NO ✗'}")

    print(f"\n  CellWarp raw distance:          {cw_result.distance:.15e}")
    print(f"    CellWarp d²:                  {cw_result.distance_squared:.15e}")
    print(f"    CellWarp scaling s:           {cw_result.scaling:.15e}")
    print(f"    norm_X (‖X_c‖_F):            {norm_X:.15e}")
    print(f"    norm_Y (‖Y_c‖_F):            {norm_Y:.15e}")

    # --- Overall verdict ---
    all_pass = pass_R and pass_aligned and pass_ref and pass_disparity
    max_of_all = max(max_diff_R, max_diff_aligned, max_diff_ref, diff_disparity)

    print(f"\n{'=' * 72}")
    print(f"OVERALL VERDICT")
    print(f"{'=' * 72}")
    print(f"  Tolerance:                {TOL:.0e}")
    print(f"  Max difference (any):     {max_of_all:.2e}")
    print(f"  Rotation matrix:          {'PASS ✓' if pass_R else 'FAIL ✗'}")
    print(f"  Aligned coordinates:      {'PASS ✓' if pass_aligned else 'FAIL ✗'}")
    print(f"  Reference coordinates:    {'PASS ✓' if pass_ref else 'FAIL ✗'}")
    print(f"  Disparity:                {'PASS ✓' if pass_disparity else 'FAIL ✗'}")
    print(f"  ──────────────────────────")
    print(f"  ALL CHECKS:               {'PASS ✓' if all_pass else 'FAIL ✗'}")

    if all_pass:
        print(f"\n  ► STAR Methods sentence:")
        print(f'    "The custom Procrustes implementation was verified against')
        print(f"     scipy.spatial.procrustes (SciPy v{scipy.__version__});")
        print(f"     rotation matrices, aligned coordinates, and disparity")
        print(f"     agreed to machine precision (max |Δ| = {max_of_all:.1e},")
        print(f'     tolerance 1e-10)."')
    else:
        print(f"\n  ⚠ DISCREPANCY DETECTED — investigate before publishing.")
        if not pass_R:
            print(f"    Rotation: max diff = {max_diff_R:.2e}")
        if not pass_aligned:
            print(f"    Aligned coords: max diff = {max_diff_aligned:.2e}")
        if not pass_ref:
            print(f"    Reference coords: max diff = {max_diff_ref:.2e}")
        if not pass_disparity:
            print(f"    Disparity: diff = {diff_disparity:.2e}")

    print()
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
