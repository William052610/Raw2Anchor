"""matcher.py

Feature-to-anchor matching for Raw2Anchor.

Matching definition
-------------------
For each feature (m/z, RT), we search for anchors within tolerances
(ppm_tol, rt_tol_min). If multiple anchors are eligible, we select the anchor
that minimizes a robust, MAD-scaled cost combining ppm and RT residuals.

Sign convention (important)
---------------------------
Residuals are defined as **(anchor − feature)**:

- d_ppm = 1e6 * (mz_anchor − mz_feature) / mz_feature
- d_rt  = rt_anchor − rt_feature

This sign convention is used consistently in `models.py` when learning and
applying correction models.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree
from models import mad

def _to_float_1d(a):
    return np.asarray(a, dtype=float).reshape(-1)

def _to_float_2d_col(a):
    a = np.asarray(a, dtype=float)
    return a.reshape(-1, 1) if a.ndim == 1 else a

def kd_match_pass(
    features_df: pd.DataFrame,
    anchor_df: pd.DataFrame,
    ppm_tol: float,
    rt_tol_min: float,
) -> pd.DataFrame:
    """Match features to anchors using a KD-Tree prefilter in the m/z domain.

    Parameters
    ----------
    features_df : pandas.DataFrame
        Required columns: ['feature_id', 'mz', 'rt'].
        - mz: m/z value (float)
        - rt: retention time in minutes (float)
    anchor_df : pandas.DataFrame
        Required columns: ['anchor_id', 'mz_ref', 'rt_ref'].
    ppm_tol : float
        m/z tolerance in ppm.
    rt_tol_min : float
        RT tolerance in minutes.

    Returns
    -------
    pandas.DataFrame
        Matched rows with columns:
        ['feature_id', 'anchor_id', 'mz_ref', 'rt_ref', 'd_ppm', 'd_rt_min']
        where d_* are defined as (anchor - feature).
    """
    # --- 1) Load anchor data ---
    A_mz = _to_float_2d_col(anchor_df['mz_ref'].values)    # (N,1)
    A_rt = _to_float_1d(anchor_df['rt_ref'].values)        # (N,)
    A_id = anchor_df['anchor_id'].astype(str).values       # (N,)

    F_id = features_df["feature_id"].to_numpy()
    # Defensive: if mz/rt columns are duplicated, pandas may return a 2D array
    mz_arr = np.asarray(features_df["mz"].to_numpy())
    rt_arr = np.asarray(features_df["rt"].to_numpy())
    if mz_arr.ndim > 1:
        mz_arr = mz_arr[:, 0]  # use the first column (or raise explicitly if you prefer strictness)
    if rt_arr.ndim > 1:
        rt_arr = rt_arr[:, 0]

    F_mz = np.asarray(mz_arr, dtype=float).reshape(-1)
    F_rt = np.asarray(rt_arr, dtype=float).reshape(-1)

    # Validate length consistency
    if not (len(F_id) == len(F_mz) == len(F_rt)):
        raise ValueError(
            f"[kd_match_pass] length mismatch: feature_id={len(F_id)}, mz={len(F_mz)}, rt={len(F_rt)}. "
            f"Ensure input has exactly one 'mz' and one 'rt' column."
        )

    # Drop NaNs
    valid = np.isfinite(F_mz) & np.isfinite(F_rt)
    valid = np.asarray(valid).reshape(-1)
    if valid.shape[0] != F_id.shape[0]:
        raise ValueError(
            f"[kd_match_pass] boolean mask length {valid.shape[0]} != data length {F_id.shape[0]}. "
            "Likely due to duplicate 'mz'/'rt' columns. Fix caller to pass only ['feature_id','mz','rt']."
        )

    F_id = F_id[valid]
    F_mz = F_mz[valid]
    F_rt = F_rt[valid]

    # --- 3) Build KD-Tree in m/z domain ---
    tree = KDTree(A_mz, metric='euclidean')

    rows = []
    rt_scale = max(mad(A_rt), 1e-6)

    # --- 4) Match each feature to an anchor ---
    for fid, mz, rt in zip(F_id, F_mz, F_rt):
        tol_da = (ppm_tol / 1e6) * mz
        q = np.array([[float(mz)]], dtype=float)
        idxs = tree.query_radius(q, r=float(tol_da))[0]

        if idxs.size == 0:
            continue

        cand_rt = A_rt[idxs]
        cand_mz = A_mz[idxs, 0]  # (k,)

        # RT filter
        mrt = np.abs(cand_rt - rt) <= rt_tol_min
        if not np.any(mrt):
            continue

        cand_rt = cand_rt[mrt]
        cand_mz = cand_mz[mrt]
        cand_id = A_id[idxs][mrt]  # prevents NameError (defensive)

        # Cost calculation
        d_ppm = 1e6 * (cand_mz - mz) / mz
        d_rt  = cand_rt - rt
        s_ppm = max(mad(d_ppm), 1.0)
        s_rt  = max(mad(d_rt), rt_scale)
        cost  = (d_ppm / s_ppm) ** 2 + (d_rt / s_rt) ** 2
        j = int(np.argmin(cost))

        rows.append({
            "feature_id": fid,
            "anchor_id": str(cand_id[j]),
            "mz_ref": float(cand_mz[j]),
            "rt_ref": float(cand_rt[j]),
            "d_ppm": float(d_ppm[j]),
            "d_rt_min": float(d_rt[j])
        })

    return pd.DataFrame(rows)
