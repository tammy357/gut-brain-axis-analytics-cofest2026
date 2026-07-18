"""
Merge targeted-extraction results with MTBLS218 sample metadata (visit
number and time), then run paired baseline-vs-1-year comparisons for
each compound.

Requires:
  multi_metabolite_extraction_results.csv (from multi_metabolite_extraction.py)
  s_MTBLS218.txt (ISA-Tab sample file, from `mtbls public download MTBLS218`)

Requires: pip install pandas numpy scipy matplotlib
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

COMPOUNDS = ["TMAO", "Hydroxyproline", "Guanidinoacetic_acid", "SDMA", "Carnitine_free"]


def classify_file(fname):
    if "_Bl.mzXML" in fname:
        return "blank"
    if "pool_QC" in fname:
        return "QC"
    return "sample"


def extract_sample_name(fname):
    stem = fname.replace(".mzXML", "")
    parts = stem.split("_", 3)
    return parts[3] if len(parts) == 4 else None


def merge_with_metadata(extraction_csv, sample_meta_txt):
    extraction = pd.read_csv(extraction_csv)
    meta = pd.read_csv(sample_meta_txt, sep="\t")
    meta_lookup = meta.set_index("Sample Name")[["Factor Value[Visit]", "Factor Value[Time]"]]

    extraction["file_type"] = extraction["file"].apply(classify_file)
    extraction["sample_name"] = extraction["file"].apply(extract_sample_name)

    merged = extraction.merge(meta_lookup, left_on="sample_name", right_index=True, how="left")
    merged = merged.rename(columns={
        "Factor Value[Visit]": "visit_number",
        "Factor Value[Time]": "time",
    })
    return merged


def sanity_check_blanks(merged, compounds=COMPOUNDS):
    print("=== Blank vs sample intensity, per compound (sanity check) ===")
    for compound in compounds:
        col = f"{compound}_intensity"
        blanks = merged[merged["file_type"] == "blank"][col]
        samples = merged[merged["file_type"] == "sample"][col]
        ratio = samples.median() / blanks.median() if blanks.median() > 0 else float("inf")
        print(f"{compound}: sample/blank ratio = {ratio:.1f}x")


def paired_comparison(merged, compounds=COMPOUNDS):
    samples = merged[merged["file_type"] == "sample"].copy()
    samples["subject"] = samples["sample_name"].str.split("_").str[0]

    results = []
    for compound in compounds:
        col = f"{compound}_intensity"
        pivot = samples.pivot_table(index="subject", columns="visit_number", values=col)
        complete = pivot.dropna()

        if len(complete) < 2 or 2.0 not in complete.columns or 6.0 not in complete.columns:
            continue

        t, p = stats.ttest_rel(complete[2.0], complete[6.0])
        fold_change = complete[6.0].mean() / complete[2.0].mean()
        results.append({
            "compound": compound,
            "n_paired": len(complete),
            "baseline_mean": complete[2.0].mean(),
            "one_year_mean": complete[6.0].mean(),
            "fold_change": fold_change,
            "p_value": p,
        })
    return pd.DataFrame(results)


def plot_summary(summary_df, out_path="human_metabolite_summary.png"):
    summary_df = summary_df.sort_values("fold_change", ascending=False)
    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ["#2b6cb0" if p < 0.05 else "#999999" for p in summary_df["p_value"]]
    ax.barh(summary_df["compound"].str.replace("_", " "), summary_df["fold_change"] - 1,
            left=1, color=colors, edgecolor="black", linewidth=0.8)
    ax.axvline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Fold-change, baseline to 1 year post-surgery")
    ax.set_title("Human plasma metabolite changes after bariatric surgery\nMTBLS218, n=44 paired subjects", fontsize=10)

    for i, (fc, p) in enumerate(zip(summary_df["fold_change"], summary_df["p_value"])):
        label = f"{fc:.2f}x" + ("  *" if p < 0.05 else "  ns")
        ax.text(fc + 0.05, i, label, va="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def main():
    merged = merge_with_metadata(
        "multi_metabolite_extraction_results.csv", "s_MTBLS218.txt"
    )
    merged.to_csv("multi_metabolite_merged.csv", index=False)

    sanity_check_blanks(merged)

    summary = paired_comparison(merged)
    print(summary.to_string(index=False))
    summary.to_csv("human_metabolite_summary.csv", index=False)

    plot_summary(summary)


if __name__ == "__main__":
    main()
