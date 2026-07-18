# Connecting the Reanalysis to the Paper's Own Behavioral & Insulin-Signaling Results

**Purpose of this document:** Tasks 1 and 4 of the original pitch (metabolites
→ behavior; metabolic nodes → insulin signaling) cannot be completed as
independent statistical analyses, because behavioral scores, cytokine
values, and insulin-signaling measurements are not in the public
Metabolomics Workbench deposit (ST000885) — confirmed directly from the
paper's own Data Availability statement: *"All data and original codes
used in this study are available from the corresponding author on
reasonable request. The complete data set of metabolites in brain and
plasma... is available at Metabolomics Workbench."*

What follows instead is a **citation-level synthesis**: the paper's own
published, quantified effects (not raw data, not digitized figures — the
specific numbers stated in its Results and Discussion text) laid alongside
our independently-reanalyzed metabolite findings. This is explicitly
**not** a new statistical result — every number below is directly
attributable to Soto et al. 2018, cited by figure/section. Where we
connect it to our own findings (TMAO, hydroxyproline, etc.), that
connection is our own interpretation, clearly marked as such.

---

## Quantified behavioral effects (from the paper's text, Figs. 2, 5)

| Test | HFD effect vs. chow | Reversed by antibiotics? | Source |
|---|---|---|---|
| Dark/light box | 26% less time in light compartment | Yes, both vancomycin and metronidazole | Fig. 2a |
| Novelty-suppressed feeding | 1.6-fold increase in latency to feed | Yes | Fig. 2b |
| Open field exploration | Decreased center entries/time; 23% less distance traveled | Yes | Fig. 2c, Suppl. Fig. 3a |
| Marble burying | 1.3-fold increase in marbles buried | Yes, toward chow-fed levels | Fig. 2d |
| Nacc cytokines (TNFα, IL-1β, IL-6, IL-10 mRNA) | Markedly increased | Yes, returned to normal | Fig. 5a |
| Nacc BDNF (pro + mature protein) | 2-fold increase | Yes | Fig. 6e |

All group sizes: n=14–18/group for behavioral tests, n=12–14/group for
cytokine mRNA (stated directly in figure legends).

## Quantified insulin-signaling effects (from the paper's text, Figs. 3, 4)

| Measurement | Chow-fed response | HFD response | Antibiotic-treated response | Source |
|---|---|---|---|---|
| Hypothalamus pIR/pIRS-1 (after insulin injection) | ~2.5-fold increase within 10 min | Little to no increase (central insulin resistance) | Returned to near-normal (both antibiotics) | Fig. 3a–c |
| Nucleus accumbens pIR/pIRS-1 | Similar pattern to hypothalamus | Similar decrease | Similar rescue | Fig. 3d,e |
| GF-mice hypothalamus pIR (via fecal transplant from chow donors) | 44-fold stimulation | Almost completely lost (HFD donor microbiota) | Robust stimulation restored (antibiotic-donor microbiota) | Fig. 4b,c |

n=5–8/group (donor mice, Fig. 3), n=4/group (GF recipient mice, Fig. 4).

---

## Where our reanalysis findings connect — and the limits of that connection

**TMAO and inflammation/insulin signaling.** Our reanalysis found TMAO
dropping with both antibiotics across all three tissues (FDR p<0.01 in
most comparisons). The paper's own Discussion (citing external literature,
not their own TMAO data — they never measured TMAO) notes that HFD
induces brain inflammation "with increased levels of cytokines and Cd11b
positive macrophages in the Nacc and VTA" and that antibiotics reduce this
in parallel with restored insulin signaling. TMAO is independently
documented in the broader literature (see gaps_and_future_work.md,
stretch-goal section — MTBLS310, and the CIRI/POCD papers found in the
TMAO literature search) as capable of activating NF-κB and promoting
neuroinflammation. **The connection we can make:** TMAO's antibiotic-driven
decline is temporally and directionally consistent with the paper's own
reported decline in Nacc cytokines under the same antibiotic treatments —
both move the same way, in the same groups, in the same direction. **What
we cannot claim:** that TMAO drives, mediates, or is mechanistically linked
to the cytokine or insulin-signaling changes. The paper measured cytokines
and insulin signaling; we measured metabolites. No single dataset ties
these together at the level of individual animals, so this is a
plausibility argument, not a demonstrated pathway.

**Hydroxyproline, carnitines, and behavior.** Our reanalysis confirms
(FDR-significant) the paper's own reported hydroxyproline and carnitine
changes with HFD. The paper's Discussion explicitly links hydroxyproline
to anxiety-like behavior via external literature (Funatsu et al., cited as
ref. 61: "Rats given fish protein dietary supplements have been shown to
have increased hydroxyproline in the brain, and this was associated with
decreased anxiety-like behaviors"). **The connection we can make:** this
gives a literature-supported, plausible mechanistic link between a
metabolite we've independently confirmed and the behavioral phenotype the
paper measured (marble burying, dark/light box) — but it is a link drawn
from a *different* study in rats given dietary supplements, not from
Soto et al.'s own mice. **What we cannot claim:** a direct, within-study
statistical relationship between hydroxyproline levels and behavioral
scores in these specific animals.

**Guanidinoacetic acid and insulin.** The paper's Discussion notes (citing
external literature, refs. 57–58) that guanidinoacetic acid "has been
shown to decrease plasma glucose, increase insulin secretion, and improve
insulin sensitivity" in the periphery. Our reanalysis confirms
guanidinoacetic acid decreases with HFD and reverses with metronidazole in
the Nacc. **The connection:** directionally consistent with an
insulin-sensitizing role, same caveat as above — external literature, not
this study's own insulin-signaling measurement tied statistically to this
metabolite.

---

## What this synthesis is — and is not

**Is:** A citation-backed, directionally-consistent narrative linking our
independently reanalyzed metabolite findings to the paper's own
behavioral and insulin-signaling results, using the specific quantified
effects the paper itself reports.

**Is not:** A statistical test, a correlation, or a causal claim. No
p-value in this document reflects a relationship between metabolites and
behavior/insulin-signaling computed from data — every number here is
either (a) directly quoted from the paper's own reported statistics on
its own behavioral/insulin data, or (b) our own reanalysis of metabolite
data alone. The two were never statistically joined, because the
underlying per-animal data needed to do so is not public.

**To make this a real statistical claim rather than a narrative**, the
per-animal behavioral scores and insulin-signaling values would need to
come from the corresponding author (see gaps_and_future_work.md) — at
which point a genuine metabolite-by-behavior and metabolite-by-insulin-
signaling correlation could be computed directly.
