"""
Rank the top non-seed metabolites by significance in each tissue x
antibiotic-comparison pair, and list every result that survives FDR
correction across the whole dataset.

Requires: st000885_fold_change_results.csv (from analyze_fold_change.py)
"""

import pandas as pd

SEED_METABOLITES = {
    "tryptophan", "GABA", "serotonin", "glutamate", "dopamine",
    "C14 carnitine", "C5-DC carnitine", "C10:2 carnitine",
    "guanidinoacetic acid", "SDMA", "hydroxyproline",
}

ANTIBIOTIC_COMPARISONS = {"HFD_vs_HFD+Metronidazole", "HFD_vs_HFD+Vancomycin"}


def top_non_seed_hits(results_df, n=5):
    non_seed = results_df[
        (~results_df["metabolite"].isin(SEED_METABOLITES))
        & (results_df["comparison"].isin(ANTIBIOTIC_COMPARISONS))
    ]
    output = {}
    for (tissue, comparison), group in non_seed.groupby(["tissue", "comparison"]):
        output[(tissue, comparison)] = group.sort_values("p_value").head(n)
    return output


def fdr_significant(results_df, threshold=0.05):
    return results_df[results_df["p_adj_fdr"] < threshold].sort_values("p_adj_fdr")


def main():
    results_df = pd.read_csv("st000885_fold_change_results.csv")

    for (tissue, comparison), top5 in top_non_seed_hits(results_df).items():
        print(f"{tissue} | {comparison}")
        for _, row in top5.iterrows():
            print(f"  {row['metabolite']:28} log2FC={row['log2_fold_change']:+.2f} "
                  f"p={row['p_value']:.2e} FDR={row['p_adj_fdr']:.4f}")
        print()

    fdr_sig = fdr_significant(results_df)
    print("FDR-significant results across the whole dataset:")
    if fdr_sig.empty:
        print("  none")
    else:
        for _, row in fdr_sig.iterrows():
            seed_flag = " [seed]" if row["metabolite"] in SEED_METABOLITES else ""
            print(f"  {row['tissue']:18} {row['comparison']:28} "
                  f"{row['metabolite']:25}{seed_flag} log2FC={row['log2_fold_change']:+.2f} "
                  f"p={row['p_value']:.2e} FDR={row['p_adj_fdr']:.4f}")


if __name__ == "__main__":
    main()
