"""
Four-panel summary combining metabolite, microbiome, cytokine, and
behavioral data across the four treatment groups (Chow / HFD /
HFD+Metronidazole / HFD+Vancomycin).

Panel A (TMAO) and Panel B (microbiome Shannon entropy) use real,
sourced values -- Panel A from the fold-change reanalysis, Panel B
read directly from Supplementary Figure 1b of Soto et al. 2018.

Panel C (Nacc cytokines) and Panel D (marble burying) are schematic.
The paper describes these as "markedly increased" / "reversed toward
chow-fed levels" without giving full group-level numbers in the text
(only a single 1.3-fold value for marble burying in HFD vs. chow), so
the bar heights here illustrate the reported pattern rather than
digitized data. This is labeled directly on the figure.

Requires: pip install matplotlib
"""

import matplotlib.pyplot as plt

GROUPS = ["Chow", "HFD", "HFD+Metro", "HFD+Vanco"]
COLORS = ["#8c9c8a", "#c0544a", "#4a6fa5", "#5a9e6f"]


def plot_tmao(ax):
    bar_groups = ["HFD\n(baseline)", "HFD+Metro", "HFD+Vanco"]
    bar_colors = [COLORS[1], COLORS[2], COLORS[3]]
    bar_heights = [1.0, 0.85, 0.85]  # log2FC = -0.23 for both antibiotics
    ax.bar(bar_groups, bar_heights, color=bar_colors, edgecolor="black", linewidth=0.8)
    ax.set_ylabel("TMAO level, relative to HFD-alone")
    ax.set_title("A. Trimethylamine-N-oxide (TMAO)\nNucleus accumbens", fontsize=10, fontweight="bold")
    ax.text(0.5, 0.95, "Both antibiotics: FDR p<0.001 vs HFD\n(log2FC = -0.23)",
            transform=ax.transAxes, ha="center", va="top", fontsize=8, style="italic")
    ax.set_ylim(0, 1.2)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.6)


def plot_microbiome_diversity(ax):
    shannon = [2.5, 2.1, 1.75, 0.3]
    shannon_err = [0.25, 0.25, 0.2, 0.1]
    ax.bar(GROUPS, shannon, yerr=shannon_err, color=COLORS, edgecolor="black",
           linewidth=0.8, capsize=4)
    ax.set_ylabel("Shannon entropy\n(microbial diversity)")
    ax.set_title("B. Gut microbiome diversity\nSupplementary Fig. 1b", fontsize=10, fontweight="bold")
    ax.text(2, 2.0, "p=0.07", ha="center", fontsize=8)
    ax.text(3, 0.55, "***", ha="center", fontsize=11, fontweight="bold")
    ax.set_ylim(0, 3)


def plot_cytokines(ax):
    cytokine_pattern = [1.0, 2.8, 1.2, 1.1]  # schematic, not digitized
    ax.bar(GROUPS, cytokine_pattern, color=COLORS, edgecolor="black", linewidth=0.8)
    ax.set_ylabel("Nacc cytokine mRNA\n(schematic — pattern only)")
    ax.set_title("C. Nacc inflammation (TNFα/IL-1β/IL-6/IL-10)\nFigure 5a — schematic", fontsize=10, fontweight="bold")
    ax.text(0.5, -0.22, "Paper states 'markedly increased' by HFD,\n'returned to normal' with either antibiotic",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5, style="italic", color="#555555")


def plot_marble_burying(ax):
    marble_relative = [1.0, 1.3, 1.05, 1.05]  # 1.3-fold is the paper's stated number
    ax.bar(GROUPS, marble_relative, color=COLORS, edgecolor="black", linewidth=0.8)
    ax.set_ylabel("Marbles buried, relative to chow")
    ax.set_title("D. Marble burying (anxiety-like behavior)\nFigure 2d", fontsize=10, fontweight="bold")
    ax.text(1, 1.35, "1.3-fold\n(stated in text)", ha="center", fontsize=7.5)
    ax.text(0.5, -0.22, "HFD value (1.3-fold) is stated in the text; antibiotic-treated\nvalues shown schematically as reversed toward chow",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5, style="italic", color="#555555")
    ax.set_ylim(0, 1.6)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    plot_tmao(axes[0, 0])
    plot_microbiome_diversity(axes[0, 1])
    plot_cytokines(axes[1, 0])
    plot_marble_burying(axes[1, 1])

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="x", labelsize=9)

    fig.suptitle(
        "Integrated summary across four treatment groups: metabolites, "
        "microbiome,\ncytokines, and behavior",
        fontsize=12, fontweight="bold", y=1.00,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("multiomic_integrated_summary.png", dpi=150, bbox_inches="tight")
    print("Saved multiomic_integrated_summary.png")


if __name__ == "__main__":
    main()
