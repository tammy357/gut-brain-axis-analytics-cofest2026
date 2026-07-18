"""
Per-metabolite fold-change and significance testing across ST000885's
donor cohort (Chow / HFD / HFD+metronidazole / HFD+vancomycin), split
by tissue (hypothalamus / nucleus accumbens / plasma).

Comparisons mirror the source publication (Soto et al. 2018, Mol
Psychiatry, DOI: 10.1038/s41380-018-0086-5):
  - Diet effect: Chow vs HFD
  - Antibiotic effect: HFD vs HFD+metronidazole, HFD vs HFD+vancomycin

Preprocessing follows the paper's stated approach: metabolites missing
in more than 80% of samples are dropped, remaining missing values are
imputed with half the per-metabolite minimum, and values are log2
transformed. Significance testing here uses Welch's t-test with
Benjamini-Hochberg FDR correction, which is not identical to the
paper's limma-based moderated t-test but follows the same logic.

Requires: st000885_metabolite_data_combined.csv (from fetch_data.py)
Requires: pip install mwtab pandas scipy numpy
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
import mwtab

COMPARISONS = [
    ("Chow_vs_HFD",              ("Chow", "none"),  ("HFD", "none")),
    ("HFD_vs_HFD+Metronidazole", ("HFD", "none"),   ("HFD", "metronidazole")),
    ("HFD_vs_HFD+Vancomycin",    ("HFD", "none"),   ("HFD", "vancomycin")),
]

SEED_METABOLITES = {
    "tryptophan", "GABA", "serotonin", "glutamate", "dopamine",
    "C14 carnitine", "C5-DC carnitine", "C10:2 carnitine",
    "guanidinoacetic acid", "SDMA", "hydroxyproline",
}


def build_sample_metadata(analysis_id="AN001442"):
    mwfile = next(mwtab.read_files(analysis_id))
    ssf = mwfile["SUBJECT_SAMPLE_FACTORS"]
    meta = {}
    for entry in ssf:
        factors = entry["Factors"]
        meta[entry["Sample ID"]] = {
            "diet": factors["Diet"],
            "treatment": factors["Treatment"],
            "tissue": factors["Tissue"],
        }
    return meta


def load_long_format(csv_path, sample_meta):
    df = pd.read_csv(csv_path)
    name_col = next(c for c in df.columns if "metabolite" in c.lower())
    sample_cols = [c for c in df.columns if c not in {name_col, "analysis_id"}]

    long_df = df.melt(id_vars=[name_col], value_vars=sample_cols,
                       var_name="sample_id", value_name="value")

    for field in ("diet", "treatment", "tissue"):
        long_df[field] = long_df["sample_id"].map(
            lambda s: sample_meta.get(s, {}).get(field)
        )

    long_df = long_df.dropna(subset=["diet"])
    long_df["value"] = pd.to_numeric(long_df["value"], errors="coerce")
    return long_df, name_col


def preprocess_tissue(tissue_df, metabolite_col):
    wide = tissue_df.pivot_table(index=metabolite_col, columns="sample_id",
                                   values="value", aggfunc="first")

    missing_frac = wide.isna().mean(axis=1)
    wide = wide.loc[missing_frac[missing_frac <= 0.80].index]

    for m in wide.index:
        row = wide.loc[m]
        if row.isna().any():
            min_val = row.min(skipna=True)
            wide.loc[m] = row.fillna(min_val / 2 if pd.notna(min_val) else np.nan)

    return np.log2(wide.where(wide > 0))


def sample_ids_for_group(tissue_df, diet, treatment):
    mask = (tissue_df["diet"] == diet) & (tissue_df["treatment"] == treatment)
    return tissue_df.loc[mask, "sample_id"].unique().tolist()


def run_comparisons(long_df, name_col):
    results = []
    for tissue in sorted(long_df["tissue"].dropna().unique()):
        tissue_df = long_df[long_df["tissue"] == tissue]
        wide_log2 = preprocess_tissue(tissue_df, name_col)

        for comp_name, (diet_a, treat_a), (diet_b, treat_b) in COMPARISONS:
            group_a = [s for s in sample_ids_for_group(tissue_df, diet_a, treat_a)
                       if s in wide_log2.columns]
            group_b = [s for s in sample_ids_for_group(tissue_df, diet_b, treat_b)
                       if s in wide_log2.columns]

            if len(group_a) < 2 or len(group_b) < 2:
                continue

            for metabolite in wide_log2.index:
                vals_a = wide_log2.loc[metabolite, group_a].dropna()
                vals_b = wide_log2.loc[metabolite, group_b].dropna()
                if len(vals_a) < 2 or len(vals_b) < 2:
                    continue

                t_stat, p_val = stats.ttest_ind(vals_a, vals_b, equal_var=False)
                results.append({
                    "tissue": tissue,
                    "comparison": comp_name,
                    "metabolite": metabolite,
                    "is_paper_seed": metabolite in SEED_METABOLITES,
                    "n_group_a": len(vals_a),
                    "n_group_b": len(vals_b),
                    "log2_fold_change": vals_b.mean() - vals_a.mean(),
                    "p_value": p_val,
                })
    return pd.DataFrame(results)


def add_fdr(results_df):
    chunks = []
    for _, group in results_df.groupby(["tissue", "comparison"]):
        group = group.copy()
        group["p_adj_fdr"] = false_discovery_control(group["p_value"])
        chunks.append(group)
    return pd.concat(chunks, ignore_index=True).sort_values(
        ["tissue", "comparison", "p_value"]
    )


def main():
    sample_meta = build_sample_metadata()
    long_df, name_col = load_long_format(
        "st000885_metabolite_data_combined.csv", sample_meta
    )

    results_df = run_comparisons(long_df, name_col)
    results_df = add_fdr(results_df)
    results_df.to_csv("st000885_fold_change_results.csv", index=False)
    print(f"Saved {results_df.shape[0]} rows to st000885_fold_change_results.csv")

    seed_results = results_df[results_df["is_paper_seed"]].sort_values(
        ["metabolite", "tissue", "comparison"]
    )
    for _, row in seed_results.iterrows():
        flag = "*" if row["p_value"] < 0.05 else " "
        print(f"{flag} {row['metabolite']:25} {row['tissue']:18} "
              f"{row['comparison']:28} log2FC={row['log2_fold_change']:+.2f} "
              f"p={row['p_value']:.4f} FDR={row['p_adj_fdr']:.4f}")


if __name__ == "__main__":
    main()
