"""anchor_map.py

Build an anchor map (reference grid) for LC–MS reverse alignment.

Key idea
--------
Traditional workflows align feature tables against each other.
Here we invert the problem: we construct a common set of reference coordinates
("anchors") directly from RAW MS1 centroid points, and then map each
software/method feature list onto those anchors. This makes cross-software and
cross-smoothing benchmarking straightforward.

Anchor construction (overview)
------------------------------
1. Read MS1 centroid points from mzML/mzXML.
2. Bucket points by m/z into fixed-width bins (default 0.001 Da).
3. Within each bin, sort by RT and split into segments when the RT gap exceeds
   `rt_gap_min`.
4. Within each segment, cluster points by `ppm_seed` tolerance.
5. For clusters with at least `min_points_per_anchor` points, define the anchor
   coordinate as the intensity-weighted median m/z and RT.

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
    ppm_seed: float = 6.0,
    min_points_per_anchor: int = 3,
    ms1_only: bool = True,
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
    ppm_seed : float, default 6.0
        m/z clustering threshold (ppm) *within* a segment.
    min_points_per_anchor : int, default 3
        Minimum number of centroid points required to emit an anchor.
    ms1_only : bool, default True
        If True, ignore MS2+ scans.
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
