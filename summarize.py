"""summarize.py

Post-processing utilities for the reverse-alignment benchmarking pipeline.

Given per-method final match tables (feature -> anchor), this module derives:

* Presence matrix (anchor × method, 0/1)
* Pairwise overlap metrics (intersection, union, Jaccard)
* Unique anchors and the corresponding per-method features
"""
import pandas as pd

def build_presence(matches_final_list) -> pd.DataFrame:
    """Build an anchor-level presence matrix.

    Parameters
    ----------
    matches_final_list : list[tuple[str, str, pandas.DataFrame]]
        Each element is (software, method, df_final).
        df_final must include an 'anchor_id' column.

    Returns
    -------
    pandas.DataFrame
        Columns: ['anchor_id', '{software}_{method}', ...] with 0/1 entries.
    """
    flags = {}
    for sw, m, df in matches_final_list:
        key = f"{sw}_{m}"
        for aid in df['anchor_id'].dropna().unique():
            flags.setdefault(aid, {})[key] = 1

    presence = pd.DataFrame.from_dict(flags, orient='index').fillna(0).astype(int)
    presence.index.name = "anchor_id"
    return presence.reset_index()

def pairwise_overlap(presence_df: pd.DataFrame) -> pd.DataFrame:
    """Compute pairwise overlap statistics from a presence matrix."""
    pres = presence_df.set_index('anchor_id')
    methods = list(pres.columns)
    rows=[]
    for i in range(len(methods)):
        for j in range(i+1, len(methods)):
            a = pres[methods[i]]; b = pres[methods[j]]
            inter = int(((a==1) & (b==1)).sum())
            union = int(((a==1) | (b==1)).sum())
            rows.append({"m1":methods[i], "m2":methods[j], "intersect":inter, "union":union,
                         "jaccard": (inter/union if union>0 else 0)})
    return pd.DataFrame(rows)

def get_unique_features(matches_final_list) -> dict:
    """Extract unique anchors (and their features) for each method.

    Return unique feature/anchor information per software_method.

    Return format:
    {
      "MSDIAL_LWMA": {
          "anchor_level": DataFrame,
          "feature_level": DataFrame
      },
      "MSDIAL_TBLWMA": {
          ...
      },
      ...
    }
    """
    # 1) Build presence matrix
    presence = build_presence(matches_final_list)
    pres = presence.set_index("anchor_id")
    methods = pres.columns.tolist()

    # 2) Map method key -> df_final
    method2df = {
        f"{sw}_{meth}": df for (sw, meth, df) in matches_final_list
    }

    result = {}

    for method_key in methods:
        # 3) Anchors present only in this method (unique anchors)
        mask = (pres[method_key] == 1) & (pres.drop(columns=[method_key]).sum(axis=1) == 0)
        unique_aids = pres.index[mask]

        if len(unique_aids) == 0:
            # Skip if there are no unique anchors
            continue

        df = method2df[method_key]

        # 4) Feature-level unique: features mapped to unique anchors in this method
        sub = df[df["anchor_id"].isin(unique_aids)].copy()

        # feature_level: use this to inspect raw feature mz/rt
        feature_level = sub[[
            "feature_id", "mz", "rt",
            "anchor_id", "mz_ref", "rt_ref",
            "d_ppm", "d_rt_min",
            "software", "method", "file"
        ]]

        # anchor_level: use this to inspect canonical anchor mz/rt
        # (If multiple features map to the same anchor, use median)
        anchor_level = (
            sub.groupby("anchor_id")
               .agg(
                   mz_ref=("mz_ref", "median"),
                   rt_ref=("rt_ref", "median"),
                   n_features=("feature_id", "nunique")
               )
               .reset_index()
        )

        result[method_key] = {
            "anchor_level": anchor_level,
            "feature_level": feature_level,
        }

    return result