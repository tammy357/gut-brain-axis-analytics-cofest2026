"""
Pull metabolite data and sample metadata for Metabolomics Workbench
study ST000885, and write it out as a combined CSV plus per-analysis
JSON metadata.

Requires: pip install mwtab pandas requests
"""

import json
import requests
import mwtab
import pandas as pd

STUDY_ID = "ST000885"


def get_analysis_ids(study_id):
    url = f"https://www.metabolomicsworkbench.org/rest/study/study_id/{study_id}/analysis"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    ids = []
    if isinstance(data, dict):
        if "analysis_id" in data:
            ids = [data["analysis_id"]]
        else:
            ids = [v["analysis_id"] for v in data.values()
                   if isinstance(v, dict) and "analysis_id" in v]
    return sorted(set(ids))


def load_analyses(analysis_ids):
    mwfiles = []
    for an_id in analysis_ids:
        try:
            mwfiles.append(next(mwtab.read_files(an_id)))
        except Exception as e:
            print(f"Could not load {an_id}: {e}")
    return mwfiles


def build_inventory(mwfiles):
    rows = []
    for mwfile in mwfiles:
        ssf = mwfile.get("SUBJECT_SAMPLE_FACTORS", [])
        study = mwfile.get("STUDY", {})
        analysis = mwfile.get("ANALYSIS", {})
        rows.append({
            "analysis_id": getattr(mwfile, "analysis_id", None),
            "study_id": getattr(mwfile, "study_id", None),
            "blocks_present": list(mwfile.keys()),
            "species": study.get("SUBJECT_SPECIES", ""),
            "analysis_type": analysis.get("ANALYSIS_TYPE", ""),
            "n_samples": len(ssf) if isinstance(ssf, list) else 0,
        })
    return pd.DataFrame(rows)


def print_sample_factors(mwfiles):
    for mwfile in mwfiles:
        an_id = getattr(mwfile, "analysis_id", "unknown")
        ssf = mwfile.get("SUBJECT_SAMPLE_FACTORS", [])
        print(f"{an_id}: {len(ssf)} samples")

        factor_names = set()
        for entry in ssf:
            factors = entry.get("Factors", {})
            if isinstance(factors, dict):
                factor_names.update(factors.keys())
        print(f"  factor fields: {sorted(factor_names)}")


def extract_metabolite_data(mwfiles):
    frames = []
    for mwfile in mwfiles:
        an_id = getattr(mwfile, "analysis_id", "unknown")
        data_block = None
        for key in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"):
            if key in mwfile:
                data_block = mwfile[key]
                break
        if data_block is None:
            continue
        rows = data_block.get("Data", [])
        if not rows:
            continue
        df = pd.DataFrame(rows)
        df["analysis_id"] = an_id
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main():
    analysis_ids = get_analysis_ids(STUDY_ID)
    print(f"Found {len(analysis_ids)} analysis ID(s) for {STUDY_ID}: {analysis_ids}")

    mwfiles = load_analyses(analysis_ids)
    print(f"Loaded {len(mwfiles)} analyses")

    inventory_df = build_inventory(mwfiles)
    print(inventory_df.to_string(index=False))

    print_sample_factors(mwfiles)

    combined = extract_metabolite_data(mwfiles)
    if not combined.empty:
        combined.to_csv("st000885_metabolite_data_combined.csv", index=False)
        print(f"Saved {combined.shape[0]} rows to st000885_metabolite_data_combined.csv")

    for mwfile in mwfiles:
        an_id = getattr(mwfile, "analysis_id", "unknown")
        with open(f"{an_id}_metadata.json", "w") as f:
            json.dump(dict(mwfile), f, indent=2, default=str)
        print(f"Saved metadata for {an_id}")


if __name__ == "__main__":
    main()
