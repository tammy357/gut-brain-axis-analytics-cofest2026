"""
Match the metabolites named in Soto et al. 2018 as significantly
changed by HFD and/or reversed by antibiotics against the metabolite
panel pulled from ST000885/AN001442.

Requires: st000885_metabolite_data_combined.csv (from fetch_data.py)
"""

import pandas as pd

PAPER_CANDIDATES = {
    "Tryptophan":             {"tissue": "Nacc",              "hfd_direction": "up",   "reversed_by": "vancomycin"},
    "GABA":                   {"tissue": "Nacc",              "hfd_direction": "up",   "reversed_by": "metronidazole"},
    "Serotonin":              {"tissue": "Nacc/plasma",       "hfd_direction": "down (brain) / up (plasma)", "reversed_by": "vancomycin (brain, trend only)"},
    "Glutamate":              {"tissue": "Nacc/plasma",       "hfd_direction": "down", "reversed_by": "vancomycin (trend only)"},
    "Dopamine":               {"tissue": "Nacc",              "hfd_direction": "down (trend only)", "reversed_by": "vancomycin (trend only)"},
    "C14 carnitine":          {"tissue": "hypothalamus",      "hfd_direction": "down", "reversed_by": "vancomycin"},
    "C5-DC carnitine":        {"tissue": "Nacc/hypothalamus", "hfd_direction": "up",   "reversed_by": "none"},
    "C10:2 carnitine":        {"tissue": "Nacc/hypothalamus", "hfd_direction": "up",   "reversed_by": "none"},
    "Guanidinoacetic acid":   {"tissue": "Nacc",              "hfd_direction": "down", "reversed_by": "metronidazole"},
    "SDMA":                   {"tissue": "hypothalamus",      "hfd_direction": "down", "reversed_by": "one or both antibiotics"},
    "Hydroxyproline":         {"tissue": "hypothalamus",      "hfd_direction": "down", "reversed_by": "vancomycin"},
}


def load_panel_names(csv_path):
    df = pd.read_csv(csv_path)
    name_col = next(c for c in df.columns if "metabolite" in c.lower())
    return name_col, set(df[name_col].astype(str).str.strip())


def match_candidates(panel_names, candidates=PAPER_CANDIDATES):
    matched, unmatched = [], []
    for candidate, info in candidates.items():
        exact = [n for n in panel_names if n == candidate]
        ci = [n for n in panel_names if n.lower() == candidate.lower()]
        substr = [n for n in panel_names
                  if candidate.lower() in n.lower() or n.lower() in candidate.lower()]
        hits = exact or ci or substr
        (matched if hits else unmatched).append((candidate, hits, info) if hits else (candidate, info))
    return matched, unmatched


def main():
    name_col, panel_names = load_panel_names("st000885_metabolite_data_combined.csv")
    print(f"Loaded panel with column '{name_col}', {len(panel_names)} unique metabolites")

    matched, unmatched = match_candidates(panel_names)

    print(f"\nMatched ({len(matched)}):")
    for candidate, hits, info in matched:
        print(f"  {candidate} -> {hits} | tissue={info['tissue']} | "
              f"HFD direction={info['hfd_direction']} | reversed by={info['reversed_by']}")

    print(f"\nNot found ({len(unmatched)}):")
    for candidate, info in unmatched:
        print(f"  {candidate} | expected tissue={info['tissue']}")


if __name__ == "__main__":
    main()
