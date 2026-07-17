"""
Volcano plots (log2 fold change vs. -log10 p-value) for every tissue x
comparison pair, with the 11 metabolites named in Soto et al. 2018
highlighted. Produces a static PNG and an interactive HTML plot per pair.

Requires: st000885_fold_change_results.csv (from analyze_fold_change.py)
Requires: pip install pandas numpy matplotlib plotly
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

SEED_METABOLITES = {
    "tryptophan", "GABA", "serotonin", "glutamate", "dopamine",
    "C14 carnitine", "C5-DC carnitine", "C10:2 carnitine",
    "guanidinoacetic acid", "SDMA", "hydroxyproline",
}


def plot_volcano(group, tissue, comparison):
    group = group.copy()
    group["neg_log10_p"] = -np.log10(group["p_value"].replace(0, np.nan))
    group["is_seed"] = group["metabolite"].isin(SEED_METABOLITES)
    group["sig_uncorrected"] = group["p_value"] < 0.05

    safe_tissue = tissue.replace(" ", "_")
    safe_comparison = comparison.replace("+", "plus")

    fig, ax = plt.subplots(figsize=(7, 6))
    non_seed = group[~group["is_seed"]]
    seed = group[group["is_seed"]]

    ax.scatter(non_seed["log2_fold_change"], non_seed["neg_log10_p"],
               c=np.where(non_seed["sig_uncorrected"], "#d97757", "#999999"),
               s=25, alpha=0.6, edgecolors="none")
    ax.scatter(seed["log2_fold_change"], seed["neg_log10_p"],
               c="#2b6cb0", s=70, edgecolors="black", linewidths=0.8,
               label="Paper-seed metabolite")

    for _, row in seed.iterrows():
        ax.annotate(row["metabolite"], (row["log2_fold_change"], row["neg_log10_p"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")

    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8,
               label="p = 0.05 (uncorrected)")
    ax.axvline(0, color="gray", linestyle="-", linewidth=0.5)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{tissue} — {comparison}", fontsize=9)
    ax.legend(fontsize=7, loc="upper left")
    plt.tight_layout()
    plt.savefig(f"volcano_{safe_tissue}_{safe_comparison}.png", dpi=150)
    plt.close(fig)

    group["marker_group"] = np.where(
        group["is_seed"], "Paper-seed metabolite",
        np.where(group["sig_uncorrected"], "p < 0.05 (uncorrected)", "not significant")
    )
    fig_html = px.scatter(
        group, x="log2_fold_change", y="neg_log10_p",
        color="marker_group", hover_name="metabolite",
        hover_data={"p_value": ":.4f", "p_adj_fdr": ":.4f",
                    "log2_fold_change": ":.2f", "marker_group": False},
        color_discrete_map={
            "Paper-seed metabolite": "#2b6cb0",
            "p < 0.05 (uncorrected)": "#d97757",
            "not significant": "#999999",
        },
        title=f"{tissue} — {comparison}",
    )
    fig_html.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
    fig_html.update_layout(xaxis_title="log2 fold change", yaxis_title="-log10(p-value)")
    fig_html.write_html(f"volcano_{safe_tissue}_{safe_comparison}.html")


def main():
    results_df = pd.read_csv("st000885_fold_change_results.csv")
    for (tissue, comparison), group in results_df.groupby(["tissue", "comparison"]):
        plot_volcano(group, tissue, comparison)
        print(f"Saved volcano plot for {tissue} / {comparison}")


if __name__ == "__main__":
    main()
