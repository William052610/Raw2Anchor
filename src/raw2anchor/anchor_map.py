"""anchor_map.py

Build an anchor map (reference grid) for Raw2Anchor.

Key idea
--------
Traditional workflows align feature tables against each other.
Here we invert the problem: we construct a common set of reference coordinates
("anchors") directly from RAW MS1 centroid points, and then map each
software/method feature list onto those anchors. Feature lists produced with
different software or smoothing methods can then be compared on the same
raw-data reference.

Anchor construction (overview)
------------------------------
1. Read MS1 centroid points from mzML/mzXML.
2. Bucket points by m/z into fixed-width bins (default 0.001 Da).
3. Within each bin, sort by RT and split into segments when the RT gap exceeds
   `rt_gap_min`.
4. Within each segment, cluster points by `ppm_seed` tolerance.
5. For clusters with at least `min_points_per_anchor` points, define the anchor
   coordinate as the intensity-weighted median m/z and RT.
6. Give one ID to anchor pairs that one ion produced on both sides of a bin
   boundary (`merge_split_anchors`).

Notes
-----
This is *raw point clustering*, not full peak picking. Anchor granularity can
change with centroiding, thresholds, and acquisition density. We recommend
reporting `mz_mad_ppm`, `rt_mad_min`, and `n_files` as anchor quality metrics.
"""
import os, re
from collections import defaultdict, namedtuple
import numpy as np
import pandas as pd

try:
    import pymzml
except Exception:
    pymzml = None
try:
    from pyteomics import mzxml
except Exception:
    mzxml = None


def _rt_iso_to_min(iso: str) -> float:
    """Convert mzXML ISO8601 duration string to minutes.

    mzXML retentionTime commonly appears as an ISO 8601 duration like:
      - PT30.5S
      - PT5.2M
      - PT1M30S

    Previous implementations sometimes handled only 'PT...S'; if 'PT...M' is parsed as 0, anchor RT values can collapse and matching will fail.
    """
    s = iso or ""
    # Capture optional minutes and seconds.
    m = re.match(r"^PT(?:(\d+(?:\.\d+)?)M)?(?:(\d+(?:\.\d+)?)S)?$", s)
    if not m:
        return 0.0
    minutes = float(m.group(1)) if m.group(1) else 0.0
    seconds = float(m.group(2)) if m.group(2) else 0.0
    return minutes + seconds / 60.0

def weighted_median(values, weights) -> float:
    """Intensity-weighted median.

    Parameters
    ----------
    values : array-like
        Values to aggregate (e.g., m/z or RT).
    weights : array-like
        Non-negative weights (typically intensity).

    Returns
    -------
    float
        Weighted median.
    """
    v = np.asarray(values)
    w = np.asarray(weights)
    idx = np.argsort(v)
    v = v[idx]; w = w[idx]
    c = np.cumsum(w) / (w.sum() if w.sum() else 1.0)
    k = np.searchsorted(c, 0.5)
    return float(v[min(k, len(v)-1)])

def mad(x) -> float:
    """Median absolute deviation (MAD).

    Robust scale estimator used for reporting anchor dispersion.
    A tiny epsilon is added to avoid division-by-zero downstream.
    """
    x = np.asarray(x)
    m = np.median(x)
    return float(np.median(np.abs(x - m)) + 1e-12)

def read_raw_points(
    raw_files,
    mz_bin_da: float = 0.001,
    rt_gap_min: float = 0.08,
    ppm_seed: float = 10.0,
    min_points_per_anchor: int = 3,
    ms1_only: bool = True,
    merge_bin_splits: bool = True,
    logger=print,
) -> pd.DataFrame:
    """Build an anchor map from RAW MS1 centroid points.

    Parameters
    ----------
    raw_files : list[str]
        Input mzML/mzXML files (centroided MS1 assumed).
    mz_bin_da : float, default 0.001
        m/z bin width (Da) for candidate bucketing.
    rt_gap_min : float, default 0.08
        Within-bin RT segmentation threshold (minutes). Larger gaps split segments.
    ppm_seed : float, default 10.0
        m/z clustering threshold (ppm) *within* a segment.
    min_points_per_anchor : int, default 3
        Minimum number of centroid points required to emit an anchor.
    ms1_only : bool, default True
        If True, ignore MS2+ scans.
    merge_bin_splits : bool, default True
        If True, merge anchors that one ion produced on both sides of an m/z bin
        boundary (see `merge_split_anchors`). Uses `ppm_seed` and `rt_gap_min`.
    logger : callable, default print
        Logger callback.

    Returns
    -------
    pandas.DataFrame
        Anchor map with columns:
        - anchor_id: unique anchor identifier
        - mz_ref, rt_ref: intensity-weighted median reference coordinates
        - n_points: number of centroid points in the cluster
        - n_files: number of RAW files contributing to the anchor
        - I_sum: sum of intensities
        - mz_mad_ppm, rt_mad_min: robust dispersion metrics
        - source_anchor_id, n_merged: pairing by `merge_split_anchors`
          (paired anchors share `anchor_id`)
    """
    if pymzml is None and mzxml is None:
        raise RuntimeError("Either pymzml or pyteomics (mzxml) is required.")

    # key = int(mz / mz_bin_da)
    # value = list[(rt_min, mz, intensity, file_index)]
    buckets = defaultdict(list)
    total = 0

    for fidx, path in enumerate(raw_files):
        logger(f"[RAW] {os.path.basename(path)}")
        ext = os.path.splitext(path)[1].lower()
        if ext == ".mzml":
            if pymzml is None: raise RuntimeError("pymzml is not installed")
            run = pymzml.run.Reader(path)
            for spec in run:
                try:
                    if ms1_only and getattr(spec, "ms_level", 1) != 1:
                        continue
                except Exception:
                    pass
                try:
                    rt = spec.scan_time_in_minutes()
                    peaks = spec.peaks("centroided")
                except Exception:
                    continue
                if peaks is None or len(peaks) == 0:
                    continue
                for mz, I in peaks:
                    key = int(mz / mz_bin_da)
                    buckets[key].append((float(rt), float(mz), float(I), fidx))
                    total += 1

        elif ext == ".mzxml":
            if mzxml is None: raise RuntimeError("pyteomics is not installed")
            with mzxml.read(path) as r:
                for spec in r:
                    if ms1_only and spec.get("msLevel", 1) != 1:
                        continue
                    rt = spec.get("retentionTime")
                    rt_min = (float(rt)/60.0) if isinstance(rt,(int,float)) else _rt_iso_to_min(rt)
                    peaks = spec.get("peaks")
                    if peaks is None: continue
                    mzs = peaks["m/z array"]; ints = peaks["intensity array"]
                    for mz, I in zip(mzs, ints):
                        key = int(mz / mz_bin_da)
                        buckets[key].append((rt_min, float(mz), float(I), fidx))
                        total += 1
        else:
            logger(f"[WARN] Unsupported extension: {path}")

    logger(f"[RAW] Loaded points: {total:,} | mz bins: {len(buckets):,}")

    def ppm_ok(mz1, mz2, ppm):
        return abs((mz1 - mz2) / mz2) * 1e6 <= ppm

    Anchor = namedtuple("Anchor", "mz_ref rt_ref n_points n_files I_sum mz_mad_ppm rt_mad_min")
    anchors = []

    for pts in buckets.values():
        pts.sort(key=lambda x: x[0])  # by RT
        seg = []; last_rt = None
        for rt, mz, I, fidx in pts:
            if (last_rt is None) or (abs(rt - last_rt) <= rt_gap_min):
                seg.append((rt, mz, I, fidx))
            else:
                _emit_anchors_from_segment(seg, ppm_seed, anchors, Anchor, min_points_per_anchor)
                seg = [(rt, mz, I, fidx)]
            last_rt = rt
        if seg:
            _emit_anchors_from_segment(seg, ppm_seed, anchors, Anchor, min_points_per_anchor)

    df = pd.DataFrame([a._asdict() for a in anchors])
    if df.empty:
        raise RuntimeError("0 anchors produced — relax parameters (ppm_seed/rt_gap_min) and confirm RT units are minutes.")
    df.insert(0, "anchor_id", [f"A{idx:07d}" for idx in range(1, len(df)+1)])
    if merge_bin_splits:
        n_before = len(df)
        df = merge_split_anchors(df, mz_bin_da=mz_bin_da, ppm_tol=ppm_seed, rt_tol_min=rt_gap_min)
        logger(f"[ANCHOR] Bin-boundary pairs given one ID: {n_before:,} rows -> {df['anchor_id'].nunique():,} anchors")
    return df


def merge_split_anchors(
    anchor_df: pd.DataFrame,
    mz_bin_da: float = 0.001,
    ppm_tol: float = 10.0,
    rt_tol_min: float = 0.08,
) -> pd.DataFrame:
    """Give one anchor ID to anchor pairs that one ion produced on both sides of an m/z bin boundary.

    Anchors are built inside fixed-width m/z bins. When the centroid m/z values
    of one ion fall on both sides of a bin boundary, the ion yields two anchors
    in adjacent bins with nearly identical m/z and RT. Different feature tables
    can then be matched to different members of the pair, so one ion would be
    counted as two anchors.

    Two anchors are paired when
      1. they lie in adjacent m/z bins,
      2. their reference m/z and RT differ by no more than `ppm_tol` and
         `rt_tol_min` (the clustering tolerance and RT gap used to build
         anchors), and
      3. each is the other's closest anchor across that boundary, with
         distance (d_ppm / ppm_tol)^2 + (d_rt / rt_tol_min)^2.
    Only pairs are formed, so dense regions cannot chain into wide anchors.
    Anchors within the same bin are never paired.

    Both members keep their own m/z and RT, so feature matching is unchanged.
    The member with the smaller summed intensity takes the `anchor_id` of the
    larger member, so features matched to either member are counted as the
    same anchor. The original ID is kept in `source_anchor_id`.

    Parameters
    ----------
    anchor_df : pandas.DataFrame
        Output of the anchor construction step (columns as in `read_raw_points`).
    mz_bin_da : float, default 0.001
        m/z bin width used to build the anchors.
    ppm_tol : float, default 10.0
        Maximum m/z difference (ppm) between paired anchors.
    rt_tol_min : float, default 0.08
        Maximum RT difference (min) between paired anchors.

    Returns
    -------
    pandas.DataFrame
        Anchor map with the same rows, `anchor_id` shared within each pair,
        and the columns `source_anchor_id` and `n_merged` (1 or 2).
    """
    df = anchor_df.reset_index(drop=True).copy()
    n = len(df)
    mz = df["mz_ref"].to_numpy(float)
    rt = df["rt_ref"].to_numpy(float)
    key = (mz / mz_bin_da).astype(np.int64)

    by_key = defaultdict(list)
    for i, k in enumerate(key):
        by_key[int(k)].append(i)

    best = np.full(n, -1)
    best_d = np.full(n, np.inf)
    for k, idx_lo in by_key.items():
        idx_hi = by_key.get(k + 1)
        if not idx_hi:
            continue
        for i in idx_lo:
            for j in idx_hi:
                d_ppm = abs(mz[j] - mz[i]) / mz[i] * 1e6
                d_rt = abs(rt[j] - rt[i])
                if d_ppm <= ppm_tol and d_rt <= rt_tol_min:
                    d = (d_ppm / ppm_tol) ** 2 + (d_rt / rt_tol_min) ** 2
                    if d < best_d[i]:
                        best_d[i], best[i] = d, j
                    if d < best_d[j]:
                        best_d[j], best[j] = d, i

    ii = np.arange(n)
    jj = best
    mutual = (jj > ii) & (best[np.where(jj >= 0, jj, 0)] == ii)
    i_m, j_m = ii[mutual], jj[mutual]
    I = df["I_sum"].to_numpy(float)
    keep = np.where(I[i_m] >= I[j_m], i_m, j_m)
    minor = np.where(I[i_m] >= I[j_m], j_m, i_m)

    ids = df["anchor_id"].astype(str).to_numpy()
    new_ids = ids.copy()
    new_ids[minor] = ids[keep]
    n_merged = np.ones(n, dtype=int)
    n_merged[keep] = 2
    n_merged[minor] = 2
    df.insert(1, "source_anchor_id", ids)
    df["anchor_id"] = new_ids
    df["n_merged"] = n_merged
    return df


def _emit_anchors_from_segment(seg, ppm_seed, anchors, Anchor, min_points):
    if len(seg) < min_points:
        return
    seg.sort(key=lambda x: x[1])  # by m/z
    cluster = [seg[0]]
    def ppm_ok(m1, m2): return abs((m1 - m2)/m2) * 1e6 <= ppm_seed
    for rt, mz, I, fidx in seg[1:]:
        if ppm_ok(mz, cluster[-1][1]):
            cluster.append((rt, mz, I, fidx))
        else:
            _emit(cluster, anchors, Anchor, min_points)
            cluster = [(rt, mz, I, fidx)]
    if cluster:
        _emit(cluster, anchors, Anchor, min_points)

def _emit(cluster, anchors, Anchor, min_points):
    if len(cluster) < min_points: return
    rts = np.array([c[0] for c in cluster])
    mzs = np.array([c[1] for c in cluster])
    Is  = np.array([c[2] for c in cluster])
    files = np.array([c[3] for c in cluster])

    mz_ref = weighted_median(mzs, Is)
    rt_ref = weighted_median(rts, Is)
    anchors.append(Anchor(
        mz_ref=float(mz_ref),
        rt_ref=float(rt_ref),
        n_points=int(len(cluster)),
        n_files=int(len(np.unique(files))),
        I_sum=float(Is.sum()),
        mz_mad_ppm=float(mad(1e6*(mzs-mz_ref)/mz_ref)),
        rt_mad_min=float(mad(rts-rt_ref))
    ))
