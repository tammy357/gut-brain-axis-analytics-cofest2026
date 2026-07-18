"""
Metabolite correlation network restricted to the 11 metabolites named
in Soto et al. 2018 (tryptophan, GABA, serotonin, glutamate, dopamine,
C14/C5-DC/C10:2 carnitines, guanidinoacetic acid, SDMA, hydroxyproline),
one network per tissue.

Nodes are colored by whether the paper reports the metabolite as
reversed by antibiotic treatment. Edge labels show the Spearman
correlation. Given n=5-6 per group, treat this as exploratory rather
than a validated biological network.

Requires: pip install mwtab pandas numpy scipy matplotlib networkx
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import networkx as nx
import mwtab

SEED_METABOLITES = [
    "tryptophan", "GABA", "serotonin", "glutamate", "dopamine",
    "C14 carnitine", "C5-DC carnitine", "C10:2 carnitine",
    "guanidinoacetic acid", "SDMA", "hydroxyproline",
]

# Direction/significance context from the paper, used to color nodes
# (blue = reversed by antibiotics per the paper; gray = not reversed)
REVERSED_BY_ABX = {
    "tryptophan": False, "GABA": True, "serotonin": False, "glutamate": False,
    "dopamine": False, "C14 carnitine": True, "C5-DC carnitine": False,
    "C10:2 carnitine": False, "guanidinoacetic acid": True, "SDMA": True,
    "hydroxyproline": True,
}

CORRELATION_THRESHOLD = 0.6  # lower than the original 0.7, since with only
                              # 11 nodes we want to see most real relationships,
                              # not just the very strongest few

# ---------------------------------------------------------------------------
# 1. Pull sample metadata and metabolite data
# ---------------------------------------------------------------------------

print("Rebuilding sample metadata mapping from AN001442...")
mwfile = next(mwtab.read_files("AN001442"))
ssf = mwfile["SUBJECT_SAMPLE_FACTORS"]
sample_meta = {
    entry["Sample ID"]: {
        "diet": entry["Factors"]["Diet"],
        "treatment": entry["Factors"]["Treatment"],
        "tissue": entry["Factors"]["Tissue"],
    }
    for entry in ssf
}

data_block = mwfile["MS_METABOLITE_DATA"]
rows = data_block["Data"]
df = pd.DataFrame(rows)
name_col = [c for c in df.columns if "metabolite" in c.lower()][0]

# Restrict to seed metabolites only, right away
df_seed = df[df[name_col].isin(SEED_METABOLITES)].copy()
print(f"Found {df_seed.shape[0]} of {len(SEED_METABOLITES)} seed metabolites in the data.")
missing = set(SEED_METABOLITES) - set(df_seed[name_col])
if missing:
    print(f"WARNING: these seed metabolites were not found in the data: {missing}")

sample_cols = [c for c in df.columns if c not in {name_col, "analysis_id"}]

# ---------------------------------------------------------------------------
# 2. Build per-tissue correlation matrices among just these 11 metabolites
# ---------------------------------------------------------------------------

tissues = sorted(set(m["tissue"] for m in sample_meta.values()))

for target_tissue in tissues:
    print(f"\n=== {target_tissue} ===")
    tissue_samples = [sid for sid, m in sample_meta.items()
                       if m["tissue"] == target_tissue and sid in sample_cols]

    wide = df_seed.set_index(name_col)[tissue_samples].apply(pd.to_numeric, errors="coerce").T

    # Same preprocessing as the fold-change analysis: impute, log2
    missing_frac = wide.isna().mean(axis=0)
    wide = wide.loc[:, missing_frac <= 0.80]
    wide = wide.apply(lambda col: col.fillna(col.min(skipna=True) / 2
                       if pd.notna(col.min(skipna=True)) else np.nan), axis=0)
    wide_log2 = np.log2(wide.where(wide > 0)).dropna(axis=1)

    if wide_log2.shape[1] < 3:
        print(f"  Skipping {target_tissue}: too few metabolites survived preprocessing.")
        continue

    corr_matrix = wide_log2.corr(method="spearman")

    G = nx.Graph()
    for m in corr_matrix.columns:
        G.add_node(m, reversed_by_abx=REVERSED_BY_ABX.get(m, False))

    edge_count = 0
    for i, m1 in enumerate(corr_matrix.columns):
        for m2 in corr_matrix.columns[i+1:]:
            r = corr_matrix.loc[m1, m2]
            if abs(r) > CORRELATION_THRESHOLD:
                G.add_edge(m1, m2, weight=r)
                edge_count += 1

    print(f"  {G.number_of_nodes()} nodes, {edge_count} edges (|r| > {CORRELATION_THRESHOLD})")

    if edge_count == 0:
        print(f"  No edges above threshold in {target_tissue} -- consider lowering "
              f"CORRELATION_THRESHOLD if you want to see this tissue's network.")
        continue

    # --- Plot: use a circular layout instead of spring_layout, since with
    #     only 11 nodes a circle is far more legible and reproducible than
    #     a force-directed layout that can look different every run ---
    pos = nx.circular_layout(G)

    fig, ax = plt.subplots(figsize=(8, 8))

    reversed_nodes = [n for n in G.nodes if G.nodes[n]["reversed_by_abx"]]
    other_nodes = [n for n in G.nodes if not G.nodes[n]["reversed_by_abx"]]

    # Edge colors: positive correlation in one color, negative in another
    pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d["weight"] > 0]
    neg_edges = [(u, v) for u, v, d in G.edges(data=True) if d["weight"] < 0]

    nx.draw_networkx_edges(G, pos, edgelist=pos_edges, ax=ax, edge_color="#2b6cb0",
                            alpha=0.5, width=2)
    nx.draw_networkx_edges(G, pos, edgelist=neg_edges, ax=ax, edge_color="#c0544a",
                            alpha=0.5, width=2, style="dashed")

    nx.draw_networkx_nodes(G, pos, nodelist=other_nodes, ax=ax,
                            node_color="#cccccc", node_size=1800, edgecolors="black")
    nx.draw_networkx_nodes(G, pos, nodelist=reversed_nodes, ax=ax,
                            node_color="#5a9e6f", node_size=1800, edgecolors="black")

    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9, font_weight="bold")

    # Edge labels showing the actual correlation value
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax, font_size=7)

    ax.set_title(
        f"Metabolite correlations — {target_tissue}\n"
        f"(the 11 metabolites named in Soto et al. 2018 only; Spearman |r| > {CORRELATION_THRESHOLD})\n"
        f"Green = reversed by antibiotics per the paper · Gray = not reversed · "
        f"Blue edge = positive correlation · Red dashed = negative\n"
        f"EXPLORATORY: n=5-6/group, not a validated biological network",
        fontsize=9,
    )
    ax.axis("off")
    plt.tight_layout()

    safe_tissue = target_tissue.replace(" ", "_")
    out_name = f"seed_correlation_network_{safe_tissue}.png"
    plt.savefig(out_name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_name}")

print("\nDone. n=5-6/group is small -- treat this as exploratory,")
print("not a validated biological network.")
