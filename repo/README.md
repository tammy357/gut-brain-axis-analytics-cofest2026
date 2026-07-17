# Gut-Brain Axis Analytics

Reanalysis of Metabolomics Workbench study [ST000885](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000885), the dataset behind:

Soto, Herzog, Pacheco et al. "Gut microbiota modulate neurobehavior through changes in brain insulin sensitivity and metabolism." *Molecular Psychiatry* 23, 2287–2301 (2018). DOI: [10.1038/s41380-018-0086-5](https://doi.org/10.1038/s41380-018-0086-5)

Original pitch: Padmashri Saravanan, CoFest 2026.

## What this is

71 mice, split across chow / high-fat diet / high-fat diet + metronidazole
/ high-fat diet + vancomycin, with 115 metabolites measured in
hypothalamus, nucleus accumbens, and plasma. This repo reanalyzes that
data independently and cross-references it against the paper's own
supplementary figures, which contain behavioral, cytokine, and
insulin-signaling results not present in the public metabolomics
deposit.

Full write-up: [`report.html`](report.html). Open it directly in a
browser — every static figure is embedded, and the six interactive
charts load Plotly from a CDN, so an internet connection is needed for
those specifically.

## Main finding

Trimethylamine-N-oxide (TMAO), a gut-bacteria-derived metabolite not
discussed in the original paper, comes out as the most consistent
signal in the reanalysis: both antibiotics lower it relative to
HFD-alone, across all three tissues, at FDR p < 0.001.

The reanalysis also independently reproduces the paper's own reported
hypothalamus/nucleus accumbens correlation (0.990 here vs. 0.994 in the
paper), which is a reasonable check that the pipeline is doing what
it's supposed to.

## Repo structure

```
scripts/    analysis pipeline, run in order below
docs/       task-by-task write-ups, including what's out of scope and why
figures/    output plots from the pipeline
report.html the full write-up with embedded figures and interactive charts
```

## Running the pipeline

```
pip install -r requirements.txt
```

Run in this order — later scripts depend on files written by earlier
ones:

1. `fetch_data.py` — pulls ST000885 from Metabolomics Workbench, writes
   `st000885_metabolite_data_combined.csv`
2. `analyze_fold_change.py` — per-metabolite fold-change and
   significance across all tissue/comparison pairs, writes
   `st000885_fold_change_results.csv`
3. `match_paper_metabolites.py` — checks which metabolites named in the
   paper are present in the pulled panel
4. `identify_top_hits.py` — ranks the strongest non-paper-named hits and
   lists everything that survives FDR correction
5. `plot_volcano.py` — volcano plots per tissue/comparison
6. `plot_pca.py` — PCA of metabolite profiles per tissue
7. `crosstissue_correlation.py` — correlates group-mean metabolite
   profiles between tissue pairs
8. `legible_seed_network.py` — correlation network restricted to the 11
   metabolites named in the paper, one per tissue
9. `build_integrated_figure.py` — four-panel summary combining
   metabolites, microbiome diversity, cytokines, and behavior (some
   panels use real values, some are schematic — see the script's
   docstring and `docs/completed_task_analyses.md`)

`fetch_data.py` and `analyze_fold_change.py`/`crosstissue_correlation.py`/
`plot_pca.py`/`legible_seed_network.py` need network access to
metabolomicsworkbench.org.

## What's in `docs/`

- **`completed_task_analyses.md`** — how each of the four pitched tasks
  was addressed, with the actual data behind each claim
- **`gaps_and_future_work.md`** — what's not achievable from public data
  and why, plus scoped next steps (16S raw-read pipeline, contacting the
  original authors for per-animal data)
- **`citation_synthesis_behavior_insulin.md`** — the paper's own
  quantified behavioral and insulin-signaling results, laid out
  alongside the metabolite findings

## Known limitations

Everything here is group-level: same four treatment groups throughout,
but metabolites, behavior, cytokines, and insulin signaling were each
measured in different animals within those groups, not one dataset
where every measurement came from the same mouse. Getting to
individual-animal resolution needs the underlying raw data, which the
paper states is available only from the corresponding author on
request.

The germ-free/fecal-transplant cohort described in the paper has no
metabolomics deposited anywhere public — confirmed by direct search of
Metabolomics Workbench. Physiological, behavioral, and insulin-signaling
results for that cohort exist in the paper's figures and are referenced
in `docs/completed_task_analyses.md`, but there's no metabolite panel to
reanalyze for those animals specifically.

16S rRNA sequencing data for these same mice exists at SRA accession
SRP132006, as raw paired-end reads rather than a processed abundance
table. Using it for a genus- or species-level analysis would need a
full amplicon pipeline (QIIME2 or DADA2) — out of scope here, logged as
future work.
