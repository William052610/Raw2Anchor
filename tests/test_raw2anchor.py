import numpy as np
import pandas as pd

from anchor_map import merge_split_anchors
from matcher import kd_match_pass
from models import apply_models, learn_models
from summarize import build_presence, upset_tables


def _anchor(aid, mz, rt, n, I):
    return dict(anchor_id=aid, mz_ref=mz, rt_ref=rt, n_points=n, n_files=1,
                I_sum=I, mz_mad_ppm=0.3, rt_mad_min=0.01)


def test_split_pair_is_merged():
    # one ion straddling the 275.038 bin boundary, plus an unrelated anchor
    A = pd.DataFrame([
        _anchor("A1", 275.03796, 0.538, 16, 1.9e7),
        _anchor("A2", 275.03827, 0.538, 104, 1.2e8),
        _anchor("A3", 300.10000, 2.000, 20, 5.0e6),
    ])
    M = merge_split_anchors(A)
    assert len(M) == 3                      # both coordinates are kept
    assert M.anchor_id.nunique() == 2       # but share one ID
    m = M.set_index("source_anchor_id")
    assert m.loc["A1", "anchor_id"] == "A2" and m.loc["A2", "anchor_id"] == "A2"
    assert m.loc["A1", "mz_ref"] == 275.03796 and m.loc["A1", "n_merged"] == 2
    assert m.loc["A3", "anchor_id"] == "A3" and m.loc["A3", "n_merged"] == 1


def test_same_bin_and_distant_anchors_are_not_merged():
    A = pd.DataFrame([
        _anchor("A1", 200.0001, 1.00, 10, 1e6),   # same bin as A2
        _anchor("A2", 200.0004, 1.00, 10, 1e6),
        _anchor("A3", 200.0012, 1.50, 10, 1e6),   # adjacent bin, RT too far
    ])
    assert merge_split_anchors(A).anchor_id.nunique() == 3


def test_only_pairs_are_merged():
    # a chain across four adjacent bins at the same RT must not collapse
    A = pd.DataFrame([_anchor(f"A{i}", 300.0005 + 0.001 * i, 1.0, 10, 1e6) for i in range(4)])
    M = merge_split_anchors(A)
    assert M.anchor_id.value_counts().max() == 2
    assert M.anchor_id.nunique() >= 2


def test_two_pass_matching_recovers_shifted_features():
    rng = np.random.default_rng(0)
    mz = np.sort(rng.uniform(100, 900, 300))
    rt = rng.uniform(0.5, 9.5, 300)
    A = pd.DataFrame({"anchor_id": [f"A{i}" for i in range(300)], "mz_ref": mz, "rt_ref": rt})
    F = pd.DataFrame({"feature_id": np.arange(300), "mz": mz * (1 - 3e-6), "rt": rt - 0.05})
    p1 = kd_match_pass(F, A, 10, 0.20)
    model = learn_models(p1, F, "test")
    Fc = apply_models(F, model)
    p2 = kd_match_pass(Fc[["feature_id", "mz_corr", "rt_corr"]].rename(columns={"mz_corr": "mz", "rt_corr": "rt"}), A, 4, 0.08)
    assert len(p2) >= 295
    hit = p2.merge(A, on="anchor_id", suffixes=("", "_a"))
    assert (hit.feature_id.map(lambda i: f"A{i}") == hit.anchor_id).mean() > 0.98


def test_presence_matrix():
    a = pd.DataFrame({"anchor_id": ["A1", "A2"]})
    b = pd.DataFrame({"anchor_id": ["A2", "A3"]})
    P = build_presence([("MSDIAL", "SG", a), ("MSDIAL", "SMA", b)]).set_index("anchor_id")
    assert P.loc["A2"].tolist() == [1, 1]
    assert P.loc["A1", "MSDIAL_SMA"] == 0


def test_features_on_either_member_share_an_anchor():
    A = merge_split_anchors(pd.DataFrame([
        _anchor("A1", 275.03796, 0.538, 16, 1.9e7),
        _anchor("A2", 275.03827, 0.538, 104, 1.2e8),
    ]))
    F = pd.DataFrame({"feature_id": [0, 1], "mz": [275.03795, 275.03828], "rt": [0.538, 0.538]})
    m = kd_match_pass(F, A, 4, 0.08)
    assert set(m.anchor_id) == {"A2"}


def test_upset_tables_counts_exclusive_intersections():
    pres = pd.DataFrame({
        "anchor_id": ["A1", "A2", "A3", "A4", "A5", "A6"],
        "MSDIAL_SG": [1, 1, 1, 1, 0, 0],
        "MZmine_LS": [1, 1, 1, 0, 1, 1],
        "MSDIAL_BF": [1, 0, 0, 0, 0, 1],
    })
    t = upset_tables(pres, method_order=["MSDIAL_BF", "MSDIAL_SG", "MZmine_LS"])
    inter = t["intersections"]
    assert list(inter.columns[4:]) == ["MSDIAL_BF", "MSDIAL_SG", "MZmine_LS"]
    assert list(inter["combination"]) == [
        "MSDIAL_SG & MZmine_LS",
        "MSDIAL_SG",
        "MZmine_LS",
        "MSDIAL_BF & MZmine_LS",
        "MSDIAL_BF & MSDIAL_SG & MZmine_LS",
    ]
    assert list(inter["intersection_size"]) == [2, 1, 1, 1, 1]
    assert list(inter["rank"]) == [1, 2, 3, 4, 5]
    assert inter["intersection_size"].sum() == len(pres)
    sizes = t["set_sizes"].set_index("method")
    assert sizes.loc["MSDIAL_SG", "set_size"] == 4
    assert sizes.loc["MZmine_LS", "unique_size"] == 1
    assert sizes.loc["MSDIAL_BF", "unique_size"] == 0
