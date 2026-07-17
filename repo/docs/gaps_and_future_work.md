# Gut-Brain Axis Analytics — Scope, Gaps & Future Work

**Project:** Gut-Brain Axis Analytics (CoFest 2026)
**Pitched by:** Padmashri Saravanan
**Dataset:** Metabolomics Workbench ST000885 / AN001442 (Project PR000885)
**Source publication:** Soto, Herzog, Pacheco et al., "Gut microbiota modulate
neurobehavior through changes in brain insulin sensitivity and metabolism,"
*Molecular Psychiatry* 23, 2287–2301 (2018). DOI: 10.1038/s41380-018-0086-5

This document tracks what each pitched task actually needed, what the
deposited data supports, and what's explicitly missing — so gaps are
documented as scoped future work rather than silently dropped.

---

## Task 1 — Map metabolites/inflammatory markers to behavioral outcomes

**What the pitch asked for:** link metabolite and inflammatory marker changes
to anxiety/depression-like behavioral phenotypes.

**What's confirmed to exist (in the source paper, not in the deposit):**
- Four named, quantified behavioral tests: dark/light box, novelty-suppressed
  feeding, open field exploration, marble burying (n=14–18/group)
- Cytokine data: TNFα, IL-1β, IL-6, IL-10 (mRNA + protein) in the nucleus
  accumbens (n=12–14/group)

**What's missing from ST000885/AN001442:**
- No behavioral scores of any kind are in the deposited metabolomics data
- No cytokine measurements are in the deposited metabolomics data

**Status:** Not achievable from this dataset alone. The behavioral and
cytokine values exist only as group-level means ± SEM in the paper's
Figures 2 and 5.

**Scoped future work:**
- Manually digitize group means/SEM from Figs. 2 and 5 (e.g. via
  WebPlotDigitizer or careful visual estimation) to enable a
  metabolite-vs-behavior correlation at the group level (n=4 groups per
  tissue, not per-animal)
- Alternative: contact the corresponding author (Julian Avila-Pacheco,
  jravilap@broadinstitute.org) — the paper states data is "available from
  the corresponding author on reasonable request"

---

## Task 2 — Compare metabolic profiles across DIO / antibiotic / germ-free+transplant

**What the pitch asked for:** cross-group comparison to identify causal
microbiota-derived metabolites, across DIO, antibiotic-treated, and
germ-free+transplant groups.

**What's confirmed to exist in the deposit:** Chow, HFD, HFD+metronidazole,
HFD+vancomycin — 71 samples, 3 tissues (hypothalamus, nucleus accumbens,
plasma), n≈6/group. Fully supported.

**What's missing:**
- **Germ-free + fecal-transplant cohort is not in any deposited dataset.**
  Confirmed via direct REST API query: `ST000885` is the only study
  returned when searching by this project's exact title — there is no
  second study or analysis covering the GF-recipient mice. The paper's own
  GF-transplant results (Figs. 1e/f, 2e, 4) are behavioral/metabolic/
  western-blot readouts only; no metabolomics table for that cohort was
  ever deposited to a public repository, as far as this search could
  determine.

**Status:** Donor-cohort comparison (Chow/HFD/antibiotics) is complete and
delivered — see `st000885_fold_change_results.csv`. Germ-free+transplant
comparison is not achievable from public data.

**Scoped future work:**
- Same author-contact route as Task 1, if germ-free metabolomics data was
  generated but never deposited
- Otherwise, report the GF-transplant causality claim as citation-only from
  the source paper, not as an independently reanalyzed result

---

## Task 3 — Multi-omic integrated visualization

**What the pitch asked for:** visualize metabolites, cytokines, and
behavioral scores together in integrated space to find co-varying patterns.

**What's confirmed to exist in the deposit:** metabolite data only, across
3 tissues and 4 groups. No cytokine or behavioral data (see Task 1).

**What was delivered instead:** a cross-tissue metabolite correlation
analysis (`crosstissue_correlation_summary.csv` and associated plots),
which reproduces the paper's own stated observation that brain regions
correlate more strongly with each other than with plasma. Result: r=0.990
(hypothalamus vs nucleus accumbens) vs. r=0.724 / r=0.708 (either brain
region vs. plasma).

**Status:** This is a legitimate single-omic, multi-tissue integration —
**not** true multi-omic integration, since only one data type (metabolites)
exists in this deposit. Any presentation of this result should state that
distinction explicitly.

**Scoped future work — the real path to multi-omic:**
- 16S rRNA community composition data for these same mice **does exist**,
  at SRA accession **SRP132006** (referenced in the paper's Methods)
- This is **raw, unassembled paired-end amplicon sequencing data**
  (~27,000 reads/sample), not a processed abundance table — confirmed via
  the paper's own methods description
- Integrating it would require a full amplicon bioinformatics pipeline
  before any taxa-vs-metabolite correlation is possible:
  1. Quality filtering and primer trimming
  2. Denoising / OTU or ASV clustering (QIIME2 or DADA2)
  3. Taxonomic assignment against a reference database (SILVA or Greengenes)
  4. Building a sample × taxon abundance table
  5. Only then: correlate taxon abundance against metabolite levels
    (e.g. does vancomycin/metronidazole's suppression of TMAO track with
    loss of specific TMA-producing bacterial taxa?)
- Realistic effort: multi-day pipeline work, not a CoFest-sprint task.
  Flagged here as the concrete next step for true multi-omic integration,
  not attempted in this pass to avoid producing unreliable taxonomy under
  time pressure

---

## Task 4 — Critical metabolic nodes affecting brain insulin signaling

**What the pitch asked for:** identify metabolic nodes where microbiota
intervention has the largest effect on brain insulin signaling and
behavior.

**What's confirmed to exist (in the paper, not the deposit):** insulin
receptor (IR) and IRS-1 phosphorylation, measured by western blot in
hypothalamus and nucleus accumbens, following in vivo insulin injection
(Figs. 3 and 4)

**What's missing from the deposit:**
- No insulin-signaling protein/phosphorylation data of any kind
- No direct annotation linking any of the 115 deposited metabolites to
  insulin-signaling pathway involvement

**What was delivered instead:** an 11-metabolite "paper-seed" panel
(tryptophan, GABA, serotonin, glutamate, dopamine, C14/C5-DC/C10:2
carnitines, guanidinoacetic acid, SDMA, hydroxyproline) — every one
explicitly named in the source paper as significantly changed by HFD
and/or reversed by antibiotics. Fold-change reanalysis confirms hydroxyproline,
C5-DC carnitine, and C10:2 carnitine as FDR-significant in this
independent reanalysis (Chow vs HFD). Additionally, **trimethylamine-N-oxide
(TMAO)** — not named in the paper's main text — was found to be
FDR-significant across every antibiotic comparison in all three tissues,
a novel finding from this reanalysis worth highlighting.

**Status:** A real, citable, hypothesis-generating metabolite panel exists.
The link from any of these metabolites to insulin signaling specifically
is established only in the source paper's western blot data — this
reanalysis cannot independently confirm a metabolite-to-insulin-signaling
causal link, only a metabolite-to-treatment-group association.

**Scoped future work:**
- Same figure-extraction route as Task 1, to get quantified pIR/pIRS-1
  values for a direct statistical link to metabolite levels
- The correlation network built here (`correlation_network.png/html`) is
  exploratory only, given n=5–6/group — should not be presented as a
  validated biological network

---

## Stretch goal — Map mouse signatures to human dietary interventions (MetaboLights)

**Status: attempted via direct MetaboLights study search. No study allows a
direct quantitative replication of the antibiotic-depletion→TMAO finding,
but two relevant human reference studies were identified.**

**Search process:** Searched metabolights-home/study-search directly for
"high fat diet" (133 hits, all mouse/rat except one tea-compound study with
mixed human/mouse serum — not a clean diet-comparison match) and
"trimethylamine" filtered to Homo sapiens (8 hits). The second search
produced genuinely relevant candidates.

**MTBLS310** — "A red wine intervention does not modify plasma
trimethylamine N-oxide but is associated with broad shifts in the plasma
metabolome and gut microbiota composition" (Haas et al., *Am J Clin Nutr*,
DOI: 10.1093/ajcn/nqac286). Randomized crossover trial, 42 men with coronary
artery disease, plasma TMAO measured by LC-MS/MS alongside 16S microbiota
sequencing. Key finding: plasma TMAO did **not** change despite
significant gut microbiota remodeling, and TMAO showed low intraindividual
concordance over time (ICC = 0.049) — explicitly flagged by the authors as
a challenge to TMAO's reliability as an individual-level biomarker.

**MTBLS218** — "An untargeted metabolomics approach to characterize
short-term and long-term metabolic changes after bariatric surgery"
(Narath et al., *PLOS ONE*, DOI: 10.1371/journal.pone.0161425). Identified
TMAO as one of several metabolites (with alanine, phenylalanine,
indoxyl-sulfate) shifting after a major diet/metabolic status change in
humans — supports TMAO as a recognized, reproducible human biomarker of
diet-driven metabolic change generally.

**Honest conclusion:** No MetaboLights study directly tests
antibiotic-induced gut depletion and TMAO in humans, so this reanalysis's
mouse finding (antibiotics lower TMAO in an HFD-obese context) cannot be
quantitatively replicated in human data. What these two studies support is
narrower and should be stated precisely: TMAO is an established,
microbiota-dependent human biomarker that responds to major dietary/
metabolic perturbations (MTBLS218) but can also be surprisingly stable
under some microbiota-altering interventions (MTBLS310) — both useful
context for interpreting the mouse result, neither a direct cross-species
replication. Any write-up should frame this as biological plausibility
context, not a conserved-signature confirmation.

---

## Summary table

| Task | Deposited data supports it? | Delivered | Key gap |
|---|---|---|---|
| 1. Metabolites/inflammation → behavior | No | Nothing yet | Behavioral/cytokine data not deposited |
| 2. Cross-group comparison | Partial (donor cohort only) | Full fold-change analysis | Germ-free+transplant not deposited anywhere |
| 3. Multi-omic visualization | No (single-omic only) | Cross-tissue correlation | True multi-omic needs 16S pipeline (SRP132006) |
| 4. Critical nodes → insulin signaling | Partial (metabolites only) | Paper-seed panel + fold-change + exploratory network | Insulin-signaling link is citation-only |
| Stretch: MetaboLights human comparison | N/A | Not attempted | — |
