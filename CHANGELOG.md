# Changelog

## Unreleased

- Replaced the UpSet-style plot (`plot_upset.py`, `upset_plot.png`) with
  `upset_intersections.xlsx`, which lists every method combination with its
  intersection size in plotting order, plus per-method set sizes. The plot can
  be drawn in any external tool. Matplotlib is no longer required.
- The third GUI button is now **Union summary**.
- Added `examples/` with the true-positive (TP) feature tables of the five
  datasets analyzed in the manuscript.

## 1.1.0 - 2026-11-08

- Anchor pairs that one ion produced on both sides of an m/z bin boundary now
  share one anchor ID (`merge_split_anchors`, applied by default). Both keep
  their coordinates, so matching is unchanged.
- Added a separate **Anchor ppm** field in the GUI. Earlier versions used the
  pass-1 m/z tolerance as the anchor clustering tolerance; the default (10 ppm)
  keeps that behavior.
- Renamed the package folder to `src/raw2anchor` and fixed the run command in
  the README.
- Added tests and continuous integration.

## 1.0.0

- First version: anchor map construction, two-pass matching with offset
  correction, presence matrix, unique features, and UpSet-style plot.
