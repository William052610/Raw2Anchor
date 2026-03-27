"""models.py

Learn and apply robust correction models using pass-1 matching residuals.

Why correction?
--------------
Different software and smoothing/peak-picking settings can induce systematic
shifts in reported (m/z, RT) even for the same compound. We use a two-pass
strategy:

1) Pass-1: match with generous tolerances to observe residuals (anchor − feature)
2) Learn: model residuals robustly (RANSAC + Huber for ppm; LOESS for RT)
3) Apply: correct feature coordinates (mz_corr, rt_corr)
4) Pass-2: re-match with tighter tolerances for final assignments

Sign convention (critical)
--------------------------
`matcher.kd_match_pass()` returns residuals as **(anchor − feature)**.
To move a feature toward its anchor, apply:

- mz_corr = mz * (1 + Δppm/1e6)
- rt_corr = rt + ΔRT

Using the opposite sign (e.g., 1 − Δppm/1e6) can degrade matching.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, HuberRegressor, RANSACRegressor
from statsmodels.nonparametric.smoothers_lowess import lowess
import inspect

# --- sklearn-version-compatible RANSAC factory ---
def make_ransac(estimator, residual_threshold, random_state: int = 42):
    """Create a RANSACRegressor with sklearn-version compatibility.

    scikit-learn 1.4+ renamed `base_estimator` -> `estimator`. This helper keeps the
    rest of the code identical across versions.
    """
    sig = inspect.signature(RANSACRegressor)
    if "estimator" in sig.parameters:  # sklearn >= 1.4
        return RANSACRegressor(estimator=estimator,
                                residual_threshold=residual_threshold,
                                random_state=random_state)
    else:  # sklearn <= 1.3
        return RANSACRegressor(base_estimator=estimator,
                                residual_threshold=residual_threshold,
                                random_state=random_state)

def mad(x) -> float:
    """Median absolute deviation (MAD), with NaN-safe median."""
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    return float(np.nanmedian(np.abs(x - med)) + 1e-12)

def _X_col(x1d):
    x = np.asarray(x1d, dtype=float).reshape(-1, 1)   # (n,1)
    return x

def _y_vec(y1d):
    y = np.asarray(y1d, dtype=float).reshape(-1)      # (n,)
    return y

def learn_models(pass1_df: pd.DataFrame, features_df: pd.DataFrame, label: str) -> dict:
    """Learn m/z and RT correction models from Pass-1 matches.

    Parameters
    ----------
    pass1_df : pandas.DataFrame
        Output of matcher.kd_match_pass(...). Must contain:
        ['feature_id', 'd_ppm', 'd_rt_min'].
    features_df : pandas.DataFrame
        Feature table for the same method. Must contain:
        ['feature_id', 'mz', 'rt'].
    label : str
        Label used for error messages/logging (e.g., "MSDIAL/SG").

    Returns
    -------
    dict
        JSON-serializable model:
        - ppm_model: robust linear regression (Huber on RANSAC inliers)
        - rt_model: LOESS curve fitted on RANSAC inliers
        - inliers: summary counts
    """
    merged = features_df.merge(pass1_df, on='feature_id', how='inner')
    if merged.empty:
        raise RuntimeError(f"[{label}] 0 pass-1 matches")

    # ---------- m/z -> Δppm (linear) ----------
    X_mz = np.asarray(merged['mz'].values, dtype=float).reshape(-1, 1)   # (n,1)
    y_ppm = np.asarray(merged['d_ppm'].values, dtype=float).reshape(-1)  # (n,)

    # drop NaN/inf rows (keeps lengths consistent)
    mask_mz = np.isfinite(X_mz).all(axis=1) & np.isfinite(y_ppm)
    X_mz = X_mz[mask_mz]
    y_ppm = y_ppm[mask_mz]
    if X_mz.shape[0] == 0:
        raise RuntimeError(f"[{label}] m/z 0 valid training samples")

    ransac_ppm = make_ransac(LinearRegression(),
                             residual_threshold=max(3.0, 2.5*mad(y_ppm)))
    ransac_ppm.fit(X_mz, y_ppm)
    in_ppm = np.asarray(ransac_ppm.inlier_mask_, dtype=bool).reshape(-1)
    if in_ppm.shape[0] != X_mz.shape[0]:
        raise RuntimeError(f"[{label}] RANSAC ppm inlier length mismatch: {in_ppm.shape[0]} vs {X_mz.shape[0]}")

    huber = HuberRegressor()
    huber.fit(X_mz[in_ppm], y_ppm[in_ppm])

    # ---------- RT -> ΔRT (nonlinear) ----------
    X_rt = np.asarray(merged['rt'].values, dtype=float).reshape(-1)       # (n,)
    y_rt = np.asarray(merged['d_rt_min'].values, dtype=float).reshape(-1) # (n,)

    mask_rt = np.isfinite(X_rt) & np.isfinite(y_rt)
    X_rt = X_rt[mask_rt]
    y_rt = y_rt[mask_rt]
    if X_rt.shape[0] == 0:
        raise RuntimeError(f"[{label}] RT 0 valid training samples")

    ransac_rt = make_ransac(LinearRegression(),
                            residual_threshold=max(0.06, 2.5*mad(y_rt)))
    ransac_rt.fit(X_rt.reshape(-1, 1), y_rt)
    in_rt = np.asarray(ransac_rt.inlier_mask_, dtype=bool).reshape(-1)
    if in_rt.shape[0] != X_rt.shape[0]:
        raise RuntimeError(f"[{label}] RANSAC RT inlier length mismatch: {in_rt.shape[0]} vs {X_rt.shape[0]}")

    lo = lowess(endog=y_rt[in_rt], exog=X_rt[in_rt], frac=0.3, it=2, return_sorted=True)
    lo_x, lo_y = lo[:, 0], lo[:, 1]

    return {
        "ppm_model": {"type": "HuberLinear",
                      "coef": float(huber.coef_[0]),
                      "intercept": float(huber.intercept_)},
        "rt_model":  {"type": "LOESS",
                      "x": lo_x.tolist(),
                      "y": lo_y.tolist(),
                      "frac": 0.3, "it": 2},
        "inliers":   {"ppm_inliers": int(in_ppm.sum()),
                      "rt_inliers": int(in_rt.sum()),
                      "total": int(merged.shape[0])}
    }


def apply_models(features_df: pd.DataFrame, model: dict) -> pd.DataFrame:
    """Apply learned correction models to a feature table.

    Adds two columns:
    - mz_corr: ppm-corrected m/z
    - rt_corr: LOESS-corrected RT (minutes)

    Notes
    -----
    The sign convention assumes d_ppm = (anchor - feature) / feature * 1e6.
    Hence, feature_mz is corrected by multiplying (1 + d_ppm/1e6).
    """
    out = features_df.copy()
    a = model['ppm_model']['coef']; b = model['ppm_model']['intercept']
    delta_ppm = a * out['mz'].to_numpy(dtype=float) + b
    # d_ppm is learned as (anchor_mz - feature_mz) / feature_mz * 1e6,
    # so to move the feature toward the anchor, apply the correction in the +delta_ppm direction.
    out['mz_corr'] = out['mz'].to_numpy(dtype=float) * (1.0 + (delta_ppm / 1e6))

    lo_x = np.asarray(model['rt_model']['x'], dtype=float)
    lo_y = np.asarray(model['rt_model']['y'], dtype=float)
    rt_raw = out['rt'].to_numpy(dtype=float)
    rt_delta = np.interp(rt_raw, lo_x, lo_y, left=lo_y[0], right=lo_y[-1])
    out['rt_corr'] = rt_raw + rt_delta
    return out
