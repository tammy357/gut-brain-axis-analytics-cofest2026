# Day 2 — Mapping to Human Interventions

Direction proposed by Padmashri Saravanan: map the mouse metabolic
fingerprints from ST000885 to real-world human nutritional/metabolic
interventions, using MetaboLights, to look for conserved biomarker
signatures.

## What was checked

Two candidate human studies were identified from an earlier MetaboLights
search for TMAO (see `gaps_and_future_work.md`):

- **MTBLS310** — red wine intervention in men with coronary artery
  disease (crossover design)
- **MTBLS218** — bariatric surgery, short- and long-term metabolic
  changes

## MTBLS218 — full comparison completed

MTBLS218's public deposit turned out to include only raw instrument
files (`.mzXML`), not a processed abundance table for most compounds.
Rather than run a full untargeted processing pipeline, five compounds
from the mouse paper-seed list were extracted directly from the raw
files at their known [M+H]+ m/z, using a 10 ppm tolerance window:
TMAO, hydroxyproline, guanidinoacetic acid, SDMA, and free carnitine
(the last one as an exploratory sanity check, not the full carnitine
ester panel — see Limitations below).

224 samples were processed (blanks, pooled QC, and 132 real subject
samples across 44 subjects at three visits: baseline, 2–3 weeks
post-surgery, and 52 weeks post-surgery). Blank samples showed
intensities 7–265x lower than real samples across all five compounds,
confirming the extraction is picking up real signal rather than noise.

### Results (baseline vs. 1 year, paired, n=44)

| Compound | Fold-change | p-value | Significant |
|---|---|---|---|
| TMAO | 2.84x ↑ | <0.00001 | Yes |
| Guanidinoacetic acid | 1.27x ↑ | 0.0001 | Yes |
| Hydroxyproline | 1.09x ↑ | 0.19 | No |
| SDMA | 1.03x ↑ | 0.45 | No |
| Carnitine (free, exploratory) | 1.03x ↑ | 0.51 | No |

See `figures/human_metabolite_summary.png`,
`figures/tmao_human_bariatric_trend.png`, and
`figures/guanidinoacetic_acid_human_trend.png`.

### How this connects to the mouse findings

Both TMAO and guanidinoacetic acid are confirmed dynamic across two
completely different gut-microbiome-altering interventions in two
different species. The direction differs between mouse and human for
both compounds, but the interventions are opposite in kind: the mouse
data measures antibiotic depletion of an already-obese gut microbiome,
while this human cohort measures bariatric surgery's long-term
remodeling of gut anatomy and microbiome composition. The honest
reading: these two compounds are genuinely microbiome-sensitive
biomarkers in both species, not that either dataset predicts the
direction in the other.

Hydroxyproline, SDMA, and carnitine showed no significant change in
this human cohort at this sample size. This does not contradict the
mouse findings (different tissue -- plasma here vs. brain and plasma in
the mouse study -- different intervention, different species) but means
these three specific compounds cannot be said to replicate in humans
from this data.

## MTBLS310 — partially checked, blocked on missing metadata

TMAO and hydroxyproline are confirmed present and quantified in this
study's MAF file (both matched by ChEBI ID: TMAO = CHEBI:15724,
hydroxyproline = CHEBI:24741), with real, non-null values across all 80
sample columns. Unlike MTBLS218, this MAF file contains pre-computed
abundances directly -- no raw-file extraction pipeline was needed here.

The comparison could not be completed because **the deposited sample
metadata (`s_MTBLS310.txt`) does not include a machine-readable factor
distinguishing which samples were collected during the red-wine period
versus the abstention period** -- confirmed by checking every column in
that file for any value that could plausibly encode treatment
assignment. This is a genuine gap in the public deposit, not a
processing failure on this end.

**To complete this comparison**, the treatment-period assignment per
sample would need to come from the study's corresponding author, or
from cross-referencing the original published paper's supplementary
tables if they list which `IDCF-` sample codes correspond to which
period.

## MTBLS469 — checked, ruled out

An apple-intake trial (Ulaszewska et al. 2020, DOI 10.1007/s00394-020-02201-8):
40 mildly hypercholesterolemic volunteers eating two whole apples a day
for 8 weeks versus a sugar- and energy-matched control beverage, in a
randomised crossover design. Plasma and urine, 320 samples. On design
grounds this is the closest human analogue to the mouse HFD experiment
found anywhere in this search: a dietary intervention with a matched
control, each subject serving as their own comparison.

**TMAO cannot be measured in this study.** The assay file states a scan
range of m/z 100-1000 for both plasma and urine. TMAO's [M+H]+ ion sits
at m/z 76.0757, below the lower bound, so it was never recorded.
Targeted extraction across all 160 plasma mzXML files returned zero
signal in zero scans, which is the expected result when the mass was
never scanned rather than evidence of absence. This is a hard
instrumental limit, not a processing failure, and it rules out the
comparison regardless of how well-designed the trial is.

**Tryptophan was checked instead** (also on the mouse paper-seed list,
and within the scan range at [M+H]+ m/z 205.0972). Signal was strong
and unambiguous: 160/160 files, median intensity 4.5e7, retention time
4.51 min against the study's stated 4.9 min.

Note for anyone reusing this MAF: its `mass_to_charge` column lists
204.0898 for tryptophan, which is the neutral monoisotopic mass, while
the `modifications` column says [M+H]+. Searching 204.0898 returns a
different, much weaker peak (median 1.1e5) eluting at 11.9 min, nowhere
near the stated retention time. The observed ion is 205.0972.

Result: **no effect.** Apple versus control beverage, paired within
subject, n=40:

| | Value |
|---|---|
| Apple mean | 44,970,923 |
| Control mean | 43,823,709 |
| Fold-change | 1.026x |
| Paired t-test | t=0.398, p=0.693 |
| Wilcoxon | W=376.0, p=0.656 |

Both tests agree, so this is not an artefact of distributional
assumptions. Eating two apples a day for eight weeks does not move bulk
plasma tryptophan. This is consistent with the original paper's own
framing, which concerned microbial co-metabolism of tryptophan rather
than tryptophan concentration.

This result also has limited bearing on the mouse tryptophan finding
either way: that was tryptophan in nucleus accumbens under HFD, and the
cross-tissue analysis in this project already showed plasma is a poor
proxy for brain (r=0.71 plasma-brain versus r=0.99 brain-brain).

## Summary

| Study | Intervention | Outcome |
|---|---|---|
| MTBLS218 | Bariatric surgery | TMAO 2.84x (p<0.00001), guanidinoacetic acid 1.27x (p=0.0001); hydroxyproline, SDMA, carnitine n.s. |
| MTBLS310 | Red wine vs abstention | Blocked: no machine-readable treatment-period labels |
| MTBLS469 | Two apples/day vs control | TMAO unmeasurable (scan range excludes m/z 76); tryptophan measured, no effect (p=0.69) |

Day 2 rests on MTBLS218. The other two are documented dead ends.

## Limitations

- **Only five compounds were checked in MTBLS218**, not the full mouse
  paper-seed panel. C14, C5-DC, and C10:2 carnitine esters were left out
  because their exact molecular formulas could not be independently
  confirmed against a citable source in the time available -- the
  "Carnitine (free)" mass checked instead is the unesterified backbone
  only, a lower-confidence stand-in, not the real target compounds.
- **Targeted extraction gives relative peak intensity, not absolute
  concentration.** Comparisons are valid within this study (same
  instrument, same samples) but are not directly comparable in absolute
  terms to the mouse study's LC-MS measurements.
- **MTBLS310's group comparison remains incomplete** for the reason
  above.
- **Only one of three studies yielded a usable comparison.** This is one
  human data point, not a survey. A broader claim about conserved
  biomarker signatures across human interventions would need more
  studies that clear both bars: measuring the compound in question, and
  depositing enough metadata to compare groups.
