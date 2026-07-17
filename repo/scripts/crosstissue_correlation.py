"""
Cross-tissue metabolite correlation across hypothalamus, nucleus
accumbens, and plasma.

For each pair of tissues, averages each metabolite's level within each
of the 12 treatment groups, then correlates the group-mean profiles
between tissues on a log2 scale. Produces a static PNG and an
interactive HTML scatter plot per tissue pair, with the 11 metabolites
named in Soto et al. 2018 highlighted.

Requires: st000885_metabolite_data_combined.csv (from fetch_data.py)
Requires: pip install pandas numpy scipy matplotlib plotly mwtab
"""

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import plotly.express as px
import mwtab

SEED_METABOLITES = {
    "tryptophan", "GABA", "serotonin", "glutamate", "dopamine",
    "C14 carnitine", "C5-DC carnitine", "C10:2 carnitine",
    "guanidinoacetic acid", "SDMA", "hydroxyproline",
}


def build_sample_metadata(analysis_id="AN001442"):
    mwfile = next(mwtab.read_files(analysis_id))
    ssf = mwfile["SUBJECT_SAMPLE_FACTORS"]
    return {
        entry["Sample ID"]: {
            "diet": entry["Factors"]["Diet"],
            "treatment": entry["Factors"]["Treatment"],
            "tissue": entry["Factors"]["Tissue"],
        }
        for entry in ssf
    }


def group_label(sample_id, sample_meta):
    m = sample_meta.get(sample_id, {})
    diet, treatment = m.get("diet"), m.get("treatment")
    if diet is None:
        return None
    return diet if treatment == "none" else f"{diet}+{treatment}"


def build_group_means(csv_path, sample_meta):
    df = pd.read_csv(csv_path)
    name_col = next(c for c in df.columns if "metabolite" in c.lower())
    sample_cols = [c for c in df.columns if c not in {name_col, "analysis_id"}]

    long_df = df.melt(id_vars=[name_col], value_vars=sample_cols,
                       var_name="sample_id", value_name="value")
    long_df["value"] = pd.to_numeric(long_df["value"], errors="coerce")
    long_df["tissue"] = long_df["sample_id"].map(
        lambda s: sample_meta.get(s, {}).get("tissue")
    )
    long_df["group"] = long_df["sample_id"].map(lambda s: group_label(s, sample_meta))
    long_df = long_df.dropna(subset=["tissue", "group"])

    group_means = (
        long_df.groupby([name_col, "tissue", "group"])["value"].mean().reset_index()
    )
    pivot = group_means.pivot_table(index=name_col, columns=["tissue", "group"],
                                      values="value")
    return pivot, sorted(long_df["tissue"].unique()), sorted(long_df["group"].unique())


def correlate_tissue_pair(pivot, tissue_a, tissue_b, groups):
    rows = []
    for metabolite in pivot.index:
        for group in groups:
            try:
                val_a = pivot.loc[metabolite, (tissue_a, group)]
                val_b = pivot.loc[metabolite, (tissue_b, group)]
            except KeyError:
                continue
            if pd.notna(val_a) and pd.notna(val_b):
                rows.append({
                    "metabolite": metabolite, "group": group,
                    f"{tissue_a}_value": val_a, f"{tissue_b}_value": val_b,
                    "is_seed": metabolite in SEED_METABOLITES,
                })
    return pd.DataFrame(rows)


def plot_pair(pair_df, tissue_a, tissue_b, r, p):
    a_col, b_col = f"{tissue_a}_value", f"{tissue_b}_value"
    pair_df = pair_df.copy()
    pair_df["log_a"] = np.log2(pair_df[a_col].where(pair_df[a_col] > 0))
    pair_df["log_b"] = np.log2(pair_df[b_col].where(pair_df[b_col] > 0))

    safe_name = f"{tissue_a}_vs_{tissue_b}".replace(" ", "_")

    fig, ax = plt.subplots(figsize=(7, 7))
    non_seed = pair_df[~pair_df["is_seed"]]
    seed = pair_df[pair_df["is_seed"]]

    ax.scatter(non_seed["log_a"], non_seed["log_b"], c="#999999", s=20, alpha=0.5)
    ax.scatter(seed["log_a"], seed["log_b"], c="#2b6cb0", s=60,
               edgecolors="black", label="Paper-seed metabolite")

    ax.set_xlabel(f"log2({tissue_a}), group mean")
    ax.set_ylabel(f"log2({tissue_b}), group mean")
    ax.set_title(f"{tissue_a} vs {tissue_b}\nSpearman r={r:.3f} (p={p:.2e})", fontsize=9)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"crosstissue_{safe_name}.png", dpi=150)
    plt.close(fig)

    pair_df["marker"] = np.where(pair_df["is_seed"], "Paper-seed metabolite", "other")
    fig_html = px.scatter(
        pair_df, x="log_a", y="log_b", color="marker",
        hover_name="metabolite", hover_data={"group": True, "marker": False},
        color_discrete_map={"Paper-seed metabolite": "#2b6cb0", "other": "#999999"},
        title=f"{tissue_a} vs {tissue_b} (Spearman r={r:.3f}, p={p:.2e})",
        labels={"log_a": f"log2({tissue_a})", "log_b": f"log2({tissue_b})"},
    )
    fig_html.write_html(f"crosstissue_{safe_name}.html")


def main():
    sample_meta = build_sample_metadata()
    pivot, tissues, groups = build_group_means(
        "st000885_metabolite_data_combined.csv", sample_meta
    )

    summary = []
    for tissue_a, tissue_b in combinations(tissues, 2):
        pair_df = correlate_tissue_pair(pivot, tissue_a, tissue_b, groups)
        if pair_df.empty:
            continue

        a_col, b_col = f"{tissue_a}_value", f"{tissue_b}_value"
        valid = (pair_df[a_col] > 0) & (pair_df[b_col] > 0)
        log_a = np.log2(pair_df.loc[valid, a_col])
        log_b = np.log2(pair_df.loc[valid, b_col])
        r, p = stats.spearmanr(log_a, log_b)

        print(f"{tissue_a} vs {tissue_b}: r={r:.3f} p={p:.2e} n={valid.sum()}")
        summary.append({"tissue_pair": f"{tissue_a}_vs_{tissue_b}",
                         "spearman_r": r, "p_value": p, "n_points": int(valid.sum())})

        plot_pair(pair_df, tissue_a, tissue_b, r, p)

    pd.DataFrame(summary).to_csv("crosstissue_correlation_summary.csv", index=False)
    print("Saved crosstissue_correlation_summary.csv")


if __name__ == "__main__":
    main()
