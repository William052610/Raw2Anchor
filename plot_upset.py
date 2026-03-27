"""plot_upset.py

Simplified UpSet-like plot for intersection sizes across methods.

This implementation intentionally avoids external dependencies (e.g., upsetplot)
to keep the repository lightweight for reviewers. For manuscript figures,
we recommend exporting the `union_presence.csv` matrix and plotting in your
preferred environment (R/ggplot2, Python/upsetplot, etc.) if additional styling
is required.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_upset_top(
    presence_df: pd.DataFrame,
    top_k: int = 25,
    outfile: str = "upset_plot.png",
    title: str = "Upset plot",
    method_order=None,
):
    """
    presence_df: columns = ['anchor_id', m1, m2, ...] (0/1)
    method_order: method order on the y-axis (e.g., ["MSDIAL_LWMA", "MSDIAL_TBLWMA", ...])
    """
    # Use anchor_id as index; remaining columns as int
    pres = presence_df.set_index("anchor_id").astype(int)

    # ---- Column (method) ordering ----
    cols = list(pres.columns)
    if method_order is not None:
        # Prioritize columns listed in method_order; append the rest
        methods = [m for m in method_order if m in cols]
        methods += [c for c in cols if c not in methods]
        pres = pres[methods]
    else:
        methods = cols

    # ---- Count combinations per anchor ----
    combos = {}
    for aid, row in pres.iterrows():
        active = tuple(m for m in methods if row[m] == 1)
        if not active:
            continue
        combos[active] = combos.get(active, 0) + 1

    # If there are no intersections
    if not combos:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No intersections", ha="center", va="center")
        ax.axis("off")
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    # Sort by intersection size (desc) and keep top_k
    items = sorted(combos.items(), key=lambda x: x[1], reverse=True)[:top_k]
    combo_list = [c for c, _ in items]
    counts = [n for _, n in items]

    # x coordinates (shared by both panels)
    x = np.arange(len(counts))

    # ---- Plot (shared x-axis) ----
    fig, (ax_bar, ax_mat) = plt.subplots(
        2, 1,
        figsize=(12, 6),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,           # IMPORTANT: share the same x-axis for both panels
    )

    # Top bar plot
    ax_bar.bar(x, counts)
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_title(title)
    ax_bar.set_xticks([])     # show x ticks only on the bottom panel

    max_count = max(counts)
    for xi, yi in zip(x, counts):
        ax_bar.text(
            xi, yi + max_count * 0.01, f"{yi}",
            ha="center", va="bottom", fontsize=9
        )

    # Bottom dot plot
    ax_mat.set_ylim(-0.5, len(methods) - 0.5)
    ax_mat.set_yticks(range(len(methods)))
    ax_mat.set_yticklabels(methods)

    # With sharex=True, xlim is inherited from the top axis
    # (Optionally set explicitly if desired)
    ax_mat.set_xlim(-0.5, len(combo_list) - 0.5)

    # For each combination: vertical line + dots
    for xi, cols_for_combo in enumerate(combo_list):
        ys = [methods.index(m) for m in cols_for_combo]
        ys.sort()
        ax_mat.plot([xi] * len(ys), ys, marker="o", linestyle="-", linewidth=1.2)

    ax_mat.set_xlabel("")  # bottom title

    # Layout adjustments
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.12, hspace=0.05)

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
