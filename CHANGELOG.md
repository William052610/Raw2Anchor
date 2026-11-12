# Changelog

## 1.0 - 2026-11-12

First public release.

- Anchor map construction from raw MS1 centroid points (mzML/mzXML).
- Anchor pairs that one ion produced on both sides of an m/z bin boundary
  share one anchor ID (`merge_split_anchors`, applied by default).
- Two-pass feature-to-anchor matching with m/z and RT offset correction.
- Presence matrix, pairwise overlap, and features unique to one method.
- `upset_intersections.xlsx` with intersection sizes in UpSet plotting order
  and per-method set sizes.
- Example true-positive (TP) feature tables of the five datasets analyzed in
  the manuscript (`examples/`).
- Tests and continuous integration.
