"""
PCA of metabolite profiles per tissue, colored by treatment group.

This uses tissue and treatment group as the integration axes. There is
no cytokine, behavioral, or insulin-signaling data in this deposit, so
this is not a multi-omic integration in the strict sense -- it's a
single-omic view across tissues and treatment groups.

Requires: st000885_metabolite_data_combined.csv (from fetch_data.py)
Requires: pip install pandas numpy matplotlib plotly scikit-learn mwtab
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import mwtab

GROUP_COLORS = {
    "Chow": "#2b6cb0", "HFD": "#d97757",
    "HFD+metronidazole": "#38a169", "HFD+vancomycin": "#805ad5",
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
        return "unknown"
    return diet if treatment == "none" else f"{diet}+{treatment}"


def run_pca_for_tissue(df, name_col, sample_meta, tissue):
    tissue_samples = [sid for sid, m in sample_meta.items()
                       if m["tissue"] == tissue and sid in df.columns]
    wide = df.set_index(name_col)[tissue_samples].T

    missing_frac = wide.isna().mean(axis=0)
    wide = wide.loc[:, missing_frac <= 0.80]
    wide = wide.apply(lambda col: col.fillna(
        col.min(skipna=True) / 2 if pd.notna(col.min(skipna=True)) else np.nan
    ), axis=0)
    wide_log2 = np.log2(wide.where(wide > 0)).dropna(axis=1)

    if wide_log2.shape[1] < 2 or wide_log2.shape[0] < 3:
        return None

    scaled = StandardScaler().fit_transform(wide_log2.values)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(scaled)

    records = []
    for i, sid in enumerate(wide_log2.index):
        records.append({
            "sample_id": sid,
            "tissue": tissue,
            "group": group_label(sid, sample_meta),
            "PC1": coords[i, 0],
            "PC2": coords[i, 1],
            "PC1_var_explained": pca.explained_variance_ratio_[0],
            "PC2_var_explained": pca.explained_variance_ratio_[1],
        })
    return records


def main():
    sample_meta = build_sample_metadata()
    df = pd.read_csv("st000885_metabolite_data_combined.csv")
    name_col = next(c for c in df.columns if "metabolite" in c.lower())

    all_records = []
    for tissue in sorted(set(m["tissue"] for m in sample_meta.values())):
        records = run_pca_for_tissue(df, name_col, sample_meta, tissue)
        if records:
            all_records.extend(records)
        else:
            print(f"Skipping PCA for {tissue}: insufficient data")

    pca_df = pd.DataFrame(all_records)
    pca_df.to_csv("st000885_pca_coordinates.csv", index=False)

    tissues = sorted(pca_df["tissue"].unique())
    fig, axes = plt.subplots(1, len(tissues), figsize=(6 * len(tissues), 5))
    if len(tissues) == 1:
        axes = [axes]

    for ax, tissue in zip(axes, tissues):
        sub = pca_df[pca_df["tissue"] == tissue]
        for grp, color in GROUP_COLORS.items():
            grp_sub = sub[sub["group"] == grp]
            ax.scatter(grp_sub["PC1"], grp_sub["PC2"], label=grp, c=color, s=60,
                       edgecolors="black", linewidths=0.5, alpha=0.85)
        var1 = sub["PC1_var_explained"].iloc[0] * 100 if len(sub) else 0
        var2 = sub["PC2_var_explained"].iloc[0] * 100 if len(sub) else 0
        ax.set_xlabel(f"PC1 ({var1:.1f}% var)")
        ax.set_ylabel(f"PC2 ({var2:.1f}% var)")
        ax.set_title(tissue)
        ax.legend(fontsize=8)

    plt.suptitle("PCA of metabolite profiles by tissue and treatment group", fontsize=10)
    plt.tight_layout()
    plt.savefig("pca_multiomic_space.png", dpi=150)
    plt.close(fig)

    fig_pca = px.scatter(
        pca_df, x="PC1", y="PC2", color="group", facet_col="tissue",
        hover_name="sample_id", color_discrete_map=GROUP_COLORS,
        title="PCA of metabolite profiles, faceted by tissue",
    )
    fig_pca.write_html("pca_multiomic_space.html")

    print("Saved pca_multiomic_space.png, pca_multiomic_space.html, "
          "st000885_pca_coordinates.csv")


if __name__ == "__main__":
    main()
