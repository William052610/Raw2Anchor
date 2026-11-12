# Raw2Anchor

Raw2Anchor is a Python GUI tool for comparing LC-MS feature tables produced
with different software or smoothing methods on a common raw-data reference.
Instead of aligning feature tables against each other, it builds a map of
**anchors** from raw MS1 centroid points and maps every feature table onto
those anchors. Features from different tables that map to the same anchor are
counted as the same signal.

Raw2Anchor was developed for the manuscript ** Smoothing Alters Feature
Composition and Low-Intensity Signal Recovery in Liquid Chromatography–Mass
Spectrometry (LC–MS) Untargeted Metabolomics**.

## Workflow

| Step | GUI button | Module | Result |
| --- | --- | --- | --- |
| A. Anchor map | **Build anchor map** | `anchor_map.py` | Anchors from raw MS1 centroid points |
| B. Two-pass matching | **Per-method matching + model + Pass-2 (KD-Tree)** | `matcher.py`, `models.py` | Feature-to-anchor assignments after offset correction |
| C. Summary | **Union summary** | `summarize.py` | Presence matrix, pairwise overlap, unique features, and an UpSet intersection table |

Smoothing and peak detection can shift the reported m/z and retention time
(RT) of the same signal, so feature tables are hard to match directly. Mapping
every table to one raw-data reference avoids choosing one software or method
as the reference.

## Installation on Windows

Raw2Anchor supports Python 3.10-3.12. A project-local virtual environment is
recommended because it keeps Raw2Anchor separate from packages already
installed in Anaconda `base` or another Python environment.

Open PowerShell or the PyCharm Terminal and move to the extracted Raw2Anchor
folder:

```powershell
cd "C:\path\to\Raw2Anchor"
```

Confirm that the selected Python version is supported:

```bash
python --version
```

Create the environment and install the dependencies:

```bash
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

On macOS or Linux, use `.venv/bin/python` in place of
`.\.venv\Scripts\python.exe`.

## GUI

Start Raw2Anchor with the Python executable inside the project environment:

```bash
.\.venv\Scripts\python.exe src\raw2anchor\gui.py
```

Activation is not required. In PyCharm, the same environment can be selected
as the project interpreter at `.venv\Scripts\python.exe`; after that, the Run
button can launch `src\raw2anchor\gui.py` directly.

1. On the **Input & Parameters** tab, click **Select RAW (mzML/mzXML)** and
   choose the raw files of one dataset.
2. Click **Add feature CSV** for each software and smoothing method.
3. Check that the software and method shown for each CSV are correct.
4. Choose an **Output folder** and check the parameters.
5. On the **Run** tab, click the three buttons in order.
6. Progress and errors are shown on the **Log** tab.

## Inputs

### Raw files

`*.mzML` or `*.mzXML` files with centroided MS1 spectra. MS2 and higher scans
are ignored. Vendor raw files can be converted with, for example,
ProteoWizard MSConvert with the `peakPicking` filter enabled for MS level 1.
Profile-mode spectra are not checked for, but they produce far more points
and are not the intended input.

### Feature tables (CSV)

The m/z and RT columns are detected by name (case-insensitive).

| Value | Accepted column names |
| --- | --- |
| m/z | `mz`, `m/z`, `mass`, `precursor_mz` |
| RT | `rt`, `retention_time`, `ret_time`, `tr` |

**RT must be in minutes.** Convert RT in seconds to minutes before loading.
Rows with a missing or non-numeric m/z or RT are dropped.

Minimal example (`Ecoli_msdial_SG.csv`):

```csv
mz,rt
132.10191,1.284
166.08626,3.512
205.09715,4.876
298.09684,6.143
```

Additional columns, such as feature IDs, intensities, or annotations, may be
kept; only the m/z and RT columns are read. One row is one feature.

Column names in tables exported from MS-DIAL or MZmine differ by software
and version (for example, `Average Mz` or `row m/z`) and are not detected
automatically. Rename the m/z and RT columns to `mz` and `rt` before loading,
or the run stops with "Could not detect mz/rt columns".

### Software and method names

The software and smoothing method are read from the CSV file name.

| File name | Software | Method |
| --- | --- | --- |
| `Ecoli_msdial_SG.csv` | `MSDIAL` | `SG` |
| `Ecoli_mzmine_LS.csv` | `MZmine` | `LS` |

- Software: the file name must contain `msdial` or `mzmine`; otherwise it is
  shown as `Unknown`.
- Method: recognized tokens are `TBLWMA`, `LWMA`, `SMA`, `BF`, `SG`, `LS`, and
  `MA`. A token is matched only between delimiters (`_`, `-`, space, or `.`);
  otherwise the method is shown as `UNK`.
- Each CSV must have a unique software–method combination. Tables with the
  same combination share one column in the presence matrix and overwrite each
  other's output files.

## Example data

The `examples/` folder contains the feature tables of the five datasets
analyzed in the manuscript. Each dataset folder holds seven tables, one per
software and smoothing method, with two columns: `m/z` and `rt` (minutes).

**These tables contain true-positive (TP) features only**, not the full
feature lists exported from MS-DIAL or MZmine. The criteria used to select TP
features are described in the Methods section of the manuscript.

The file names already follow the naming rule, so each table is assigned the
correct software and method when it is added in the GUI.

The raw mzML/mzXML files are not included because of their size, so Step A
cannot be run from this folder alone. To reproduce a comparison, select the
raw files of one dataset with **Select RAW (mzML/mzXML)** and add all seven
tables from the matching folder with **Add feature CSV**.

## How it works

### Step A. Anchor map

1. Raw MS1 centroid points from all selected files are grouped into fixed m/z
   bins (`mz bin`).
2. Within each bin, points are sorted by RT and split where the gap between
   consecutive points exceeds `RT gap`.
3. Within each segment, points are sorted by m/z and joined while adjacent
   points are within `Anchor ppm`.
4. Each cluster with at least `min pts/anchor` points becomes an anchor at its
   intensity-weighted median m/z and RT.
5. **Split anchors share one ID.** When the centroids of one ion fall on both
   sides of a bin boundary, the ion yields two anchors in adjacent bins. Two
   such anchors are paired when their m/z and RT differ by no more than
   `Anchor ppm` and `RT gap` and each is the other's closest anchor across
   that boundary. Only pairs are formed. Both members keep their own m/z and
   RT, so matching is unchanged, but they share the ID of the member with the
   larger summed intensity. `source_anchor_id` keeps the original ID.

Each anchor also stores its dispersion (`mz_mad_ppm`, `rt_mad_min`), number of
centroid points (`n_points`), summed intensity (`I_sum`), and number of
contributing raw files (`n_files`).

### Step B. Two-pass matching

Residuals are defined as anchor − feature:

```text
d_ppm = 1e6 × (mz_anchor − mz_feature) / mz_feature
d_rt  = rt_anchor − rt_feature
```

**Pass 1** uses wide tolerances. Anchors within the m/z tolerance are found
with a KD-tree and filtered by RT. When several anchors remain, the one with
the lowest MAD-scaled cost is chosen:

```text
cost = (d_ppm / s_ppm)² + (d_rt / s_rt)²
s_ppm = max(MAD(d_ppm), 1 ppm)
s_rt  = max(MAD(d_rt), MAD of all anchor RTs)
```

**Offset correction.** Systematic shifts are modeled from the pass-1
residuals of each feature table.

| Coordinate | Model | Applied as |
| --- | --- | --- |
| m/z | RANSAC (linear, threshold max(3 ppm, 2.5 × MAD)), then Huber regression of Δppm on m/z using the inliers | `mz_corr = mz × (1 + Δppm / 1e6)` |
| RT | RANSAC (linear, threshold max(0.06 min, 2.5 × MAD)), then LOESS of ΔRT on RT using the inliers (`statsmodels` lowess, frac 0.3, 2 iterations) | `rt_corr = rt + ΔRT` |

ΔRT outside the fitted RT range is held at the nearest end value. A feature
table with no pass-1 matches is skipped and reported in the log.

**Pass 2** repeats the matching with tight tolerances on the corrected
coordinates and gives the final feature-to-anchor assignments.

### Step C. Summary

Each table is labeled `<Software>_<Method>`. An anchor is present for a table
when at least one of its features is assigned to that anchor in pass 2.

Raw2Anchor does not draw plots. Instead, `upset_intersections.xlsx` contains
the numbers needed for an UpSet plot, already in plotting order.

**`intersections` sheet** — one row per combination of methods. Each anchor is
counted once, in the exact combination of methods that found it, so the
intersection sizes add up to the total number of anchors. Example (numbers
are illustrative):

| rank | combination | n_methods | intersection_size | MSDIAL_SG | MZmine_LS | MSDIAL_BF |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | MSDIAL_SG & MZmine_LS | 2 | 1520 | 1 | 1 | 0 |
| 2 | MSDIAL_SG & MZmine_LS & MSDIAL_BF | 3 | 1184 | 1 | 1 | 1 |
| 3 | MZmine_LS | 1 | 412 | 0 | 1 | 0 |

- Rows are sorted by `intersection_size` (largest first); ties are ordered by
  `n_methods`, then by `combination`. Use `rank` as the x-axis order.
- `intersection_size` is the bar height.
- The 0/1 method columns are the dot matrix below the bars. They follow the
  order in which the CSV files were added in the GUI.
- All combinations are listed. For a top-*N* plot, keep the first *N* rows.

**`set_sizes` sheet** — one row per method, in the same order: `set_size` is
the number of anchors found by the method (the horizontal set-size bars), and
`unique_size` is the number found only by that method.

## Parameters and defaults

| GUI field | Meaning | Default |
| --- | --- | --- |
| mz bin (Da) | m/z bin width for raw points | 0.001 |
| RT gap (min) | RT gap that splits a segment | 0.08 |
| Anchor ppm | m/z clustering tolerance within a segment, also used to pair split anchors | 10 |
| min pts/anchor | Minimum centroid points per anchor | 3 |
| Pass-1 ppm | Pass-1 m/z tolerance | 10 |
| Pass-1 RT (min) | Pass-1 RT tolerance | 0.20 |
| Pass-2 ppm | Pass-2 m/z tolerance | 4 |
| Pass-2 RT (min) | Pass-2 RT tolerance | 0.08 |

These defaults were used in the manuscript.

## Outputs

All files are written to the selected output folder.

| File | Step | Content |
| --- | --- | --- |
| `anchor_map.xlsx` | A | Anchor ID, reference m/z and RT, and anchor quality columns |
| `matches_<SW>_<METHOD>_pass1.csv` | B | Pass-1 feature-to-anchor matches |
| `matches_<SW>_<METHOD>_final.csv` | B | All features with their pass-2 anchor (blank if unmatched) |
| `model_<SW>_<METHOD>.json` | B | m/z and RT correction models and inlier counts |
| `union_presence.csv` | C | Anchor × method presence matrix (0/1) |
| `pairwise_overlap.csv` | C | Intersection, union, and Jaccard index for each pair of methods |
| `unique_anchor_<SW>_<METHOD>.csv` | C | Anchors found by only this method |
| `unique_feature_<SW>_<METHOD>.csv` | C | Features of this method assigned to those unique anchors |
| `upset_intersections.xlsx` | C | UpSet plot data: intersection sizes in plotting order and set sizes |

## Tests

```bash
.\.venv\Scripts\python.exe -m pip install pytest
.\.venv\Scripts\python.exe -m pytest
```

The tests check the pairing of split anchors, the two-pass matching, and the
UpSet intersection counts on synthetic data.

## Scope and limitations

- Raw2Anchor compares feature tables; it is not a feature detector or a
  metabolite identification program.
- Anchors come from clustering raw centroid points, not from peak picking.
  Very low-intensity noise can add anchors.
- Only pairs of anchors share an ID. A signal spread over more than two m/z
  bins can still yield more than one anchor.
- Keep RT in minutes throughout.

## Citation

Citation metadata are provided in `CITATION.cff`.

Author: Gun Woo Lee, Konkuk University  
Contact: [dlrjsdn0526@konkuk.ac.kr](mailto:dlrjsdn0526@konkuk.ac.kr)  
Repository: <https://github.com/William052610/Raw2Anchor>

## License

Raw2Anchor is distributed under the [MIT License](LICENSE).

Copyright (c) 2026 Gun Woo Lee.
