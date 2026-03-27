# LC‑MS Reverse Alignment for Smoothing / Peak‑Picking Benchmarking

This repository implements a **reverse‑alignment** workflow to benchmark LC–MS
feature tables produced by different software and/or smoothing methods.

Instead of aligning feature tables against each other, we:

1) build a common reference grid (**anchor map**) from RAW MS1 centroid points, and
2) map each method’s **feature list** to those anchors using a two‑pass matching and
   correction procedure.

The outputs (union/unique anchors, overlaps, UpSet plot) can be used to quantify
how smoothing settings affect feature detection downstream.

---

## Quick summary

- Build an **anchor map** from RAW MS1 centroid points (mzML/mzXML).
- Map each software/method feature table (CSV) to anchors via two-pass matching and robust correction.
- Quantify overlaps/uniques and visualize intersection sizes using an UpSet-like plot.

---

## Repository layout

| File | Purpose |
|---|---|
| `gui.py` | Tkinter GUI (end‑to‑end pipeline) |
| `anchor_map.py` | Build anchor map from RAW MS1 centroid points |
| `matcher.py` | KD‑Tree based feature→anchor matching |
| `models.py` | Robust correction models (ppm + RT) |
| `summarize.py` | Presence matrix, pairwise overlaps, unique anchors/features |
| `plot_upset.py` | Lightweight UpSet‑like plot |
| `requirements.txt` | Python dependencies |

---

## Installation

### 1) Create a Python environment

Python 3.9+ is recommended.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install -r requirements.txt
```

### 2) Run the GUI

```bash
python gui.py
```

---

## Inputs

### 1) RAW files

- Supported: `*.mzML`, `*.mzXML`
- Assumption: MS1 centroid peaks are available (or retrievable via the reader).

### 2) Feature tables (CSV)

The GUI auto‑detects `mz` and `rt` columns using common names:

- m/z column candidates: `mz`, `m/z`, `mass`, `precursor_mz`
- RT column candidates: `rt`, `retention_time`, `ret_time`, `tr`

**Important:** RT must be in **minutes** for correct matching.
If your feature table has RT in seconds, convert to minutes before loading.

### 3) Method naming (optional)

The GUI tries to infer `(software, method)` from the filename.
For example:

- `Ecoli_msdial_SG.csv` → software=`MSDIAL`, method=`SG`
- `Ecoli_mzmine_LWMA.csv` → software=`MZmine`, method=`LWMA`

Tokens are matched using a delimiter‑aware regex to reduce false positives (e.g., `MA`).

---

## Pipeline details

### Step A — Anchor map construction (`anchor_map.py`)

Anchor map is built from RAW MS1 centroid points using:

1) **m/z bucketing** with width `mz_bin_da` (Da)
2) **RT segmentation** within each bucket using `rt_gap_min` (minutes)
3) **ppm clustering** within each segment using `ppm_seed`
4) **Intensity‑weighted medians** to define anchor coordinates:

- `mz_ref = weighted_median(mz, intensity)`
- `rt_ref = weighted_median(rt, intensity)`

Each anchor also stores dispersion metrics (`mz_mad_ppm`, `rt_mad_min`) and how many
RAW files contribute (`n_files`).

### Step B — Two‑pass matching (`matcher.py`, `models.py`)

For each method’s feature list:

**Pass‑1 (generous tolerances)**
- Find anchor candidates by KD‑Tree radius search in m/z
- Filter candidates by RT
- Choose the best anchor using a robust scaled cost:

\[ d_{ppm} = 10^6 \cdot \frac{(mz_{anchor} - mz_{feature})}{mz_{feature}} \]
\[ d_{rt} = rt_{anchor} - rt_{feature} \]
\[ cost = (d_{ppm}/MAD(d_{ppm}))^2 + (d_{rt}/MAD(d_{rt}))^2 \]

**Model learning (robust correction)**
- ppm correction: RANSAC(LinearRegression) → Huber regression on inliers
- RT correction: RANSAC(LinearRegression) → LOESS(LOWESS) on inliers

**Apply correction**
- `mz_corr = mz * (1 + d_ppm_pred/1e6)`
- `rt_corr = rt + d_rt_pred`

**Pass‑2 (tight tolerances)**
- Repeat matching with tighter ppm/RT to obtain final mapping.

> Sign convention note: `d_ppm` and `d_rt` are defined as **(anchor − feature)**.

### Step C — Union/unique anchors and UpSet (`summarize.py`, `plot_upset.py`)

- `union_presence.csv`: anchor × method presence matrix
- `pairwise_overlap.csv`: intersection/union/Jaccard for all method pairs
- `unique_anchor_*.csv`, `unique_feature_*.csv`: anchors/features unique to a single method
- `upset_plot.png`: top intersections (lightweight UpSet‑like figure)

---

## GUI parameters (recommended defaults)

| Parameter | Meaning | Typical |
|---|---|---|
| 1st ppm | Pass‑1 m/z tolerance (ppm) | 10 |
| 1st RT | Pass‑1 RT tolerance (min) | 0.20 |
| 2nd ppm | Pass‑2 m/z tolerance (ppm) | 4 |
| 2nd RT | Pass‑2 RT tolerance (min) | 0.08 |
| mz bin (Da) | RAW bucketing width | 0.001 |
| RT gap (min) | segmentation gap | 0.08 |
| min pts/anchor | minimum points to form an anchor | 3 |

---

## Outputs

All outputs are written to the selected output folder.

- `anchor_map.xlsx`
- `matches_<SW>_<METHOD>_pass1.csv`
- `matches_<SW>_<METHOD>_final.csv`
- `model_<SW>_<METHOD>.json`
- `union_presence.csv`
- `pairwise_overlap.csv`
- `unique_anchor_<SW>_<METHOD>.csv`
- `unique_feature_<SW>_<METHOD>.csv`
- `upset_plot.png`

---

## Manuscript‑ready Methods paragraph (copy/paste template)

> We constructed a reference “anchor map” directly from RAW MS1 centroid data (mzML/mzXML).
> Centroid points were bucketed by m/z (0.001 Da), segmented by retention time (RT gap 0.08 min),
> and clustered within segments using a 6‑ppm threshold. For each cluster with ≥3 points, the
> anchor coordinate was defined as the intensity‑weighted median m/z and RT. Feature tables generated
> by each software/smoothing method were mapped to the anchor map using a two‑pass matching strategy.
> In pass‑1, anchors within ±10 ppm and ±0.20 min were queried using a KD‑Tree in the m/z domain and
> the best candidate was selected by a robust cost using median absolute deviation (MAD) scaling.
> Systematic shifts were modeled from pass‑1 residuals using robust regression (RANSAC followed by Huber
> regression for ppm errors) and LOESS smoothing for RT errors, then applied as coordinate corrections
> (mz_corr = mz·(1 + Δppm/10^6), rt_corr = rt + ΔRT). Pass‑2 matching was performed with tighter tolerances
> (±4 ppm, ±0.08 min) to obtain final feature‑to‑anchor assignments. Anchor‑level presence matrices were
> used to compute method overlaps, unique anchors, and intersection sizes visualized by an UpSet‑like plot.

---

## Notes / limitations

- Anchor map is **raw point clustering**, not full peak‑picking. Very low intensity noise can inflate anchors.
  Consider adding an MS1 intensity threshold or top‑N filtering if RAW files are extremely large.
- Ensure RT units are consistent (minutes recommended throughout).

---

## Contact / attribution

If you include this repository or the anchor-based benchmarking concept in a paper, please cite the associated manuscript/code release.
